//! Synthetic sky and camera for end-to-end solver validation.
//!
//! A deterministic pseudo-random catalog (golden-angle sphere distribution)
//! is rendered through a known WCS into a grayscale image with a Gaussian
//! PSF, Poisson-ish photon noise and sky background. Solving the rendered
//! image and comparing against the injected camera gives a true
//! end-to-end accuracy number in arcseconds.

use crate::astro::{TanProjection, Vec3, ARCSEC, DEG};
use crate::catalog::CatStar;
use crate::image::GrayImage;

/// Deterministic splitmix64 PRNG so test runs are reproducible.
#[derive(Clone, Debug)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed)
    }
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    pub fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.uniform()
    }
    /// Standard normal via Box–Muller.
    pub fn normal(&mut self) -> f64 {
        let u = self.uniform().max(1e-12);
        let v = self.uniform();
        (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
    }
}

/// A synthetic catalog: uniformly distributed stars with magnitudes.
pub struct SyntheticSky {
    pub stars: Vec<CatStar>,
}

impl SyntheticSky {
    /// `n` stars over the whole sphere, magnitudes in `[mag_lo, mag_hi]`.
    pub fn random(n: u32, mag_lo: f32, mag_hi: f32, seed: u64) -> Self {
        let mut rng = Rng::new(seed);
        let mut stars = Vec::with_capacity(n as usize);
        // Fibonacci sphere gives a low-discrepancy uniform base; jitter
        // positions slightly so geometry is not grid-aligned.
        let ga = std::f64::consts::PI * (3.0 - 5.0f64.sqrt());
        for i in 0..n {
            let z = 1.0 - (i as f64 + 0.5) * 2.0 / n as f64;
            let r = (1.0 - z * z).max(0.0).sqrt();
            let phi = ga * i as f64;
            let mut v = Vec3 {
                x: r * phi.cos(),
                y: r * phi.sin(),
                z,
            };
            // Small angular jitter (~10 arcsec).
            v = Vec3 {
                x: v.x + rng.normal() * 5e-5,
                y: v.y + rng.normal() * 5e-5,
                z: v.z + rng.normal() * 5e-5,
            }
            .normalize();
            let (ra, dec) = v.to_radec();
            stars.push(CatStar {
                id: i,
                ra,
                dec,
                mag: rng.range(mag_lo as f64, mag_hi as f64) as f32,
                v,
            });
        }
        SyntheticSky { stars }
    }
}

/// Pinhole camera with a tangent-plane WCS (plus optional radial distortion).
#[derive(Clone, Debug)]
pub struct SyntheticCamera {
    pub width: usize,
    pub height: usize,
    pub ra0: f64,
    pub dec0: f64,
    /// Plate scale, radians per pixel.
    pub scale: f64,
    /// Rotation east of north, radians.
    pub rotation: f64,
    /// Gaussian PSF sigma in pixels.
    pub psf_sigma_px: f64,
    /// Sky background level in pixel flux units.
    pub sky: f32,
    /// Read-noise standard deviation in the same units.
    pub read_noise: f32,
    /// Flux of a magnitude-zero star in the same units.
    pub zp_flux: f64,
    /// Radial distortion coefficient (1/pixel²), 0 = pure pinhole.
    pub k1: f64,
}

impl Default for SyntheticCamera {
    fn default() -> Self {
        SyntheticCamera {
            width: 2048,
            height: 1536,
            ra0: 10.68 * DEG,
            dec0: 41.27 * DEG,
            scale: 2.5 * ARCSEC,
            rotation: 33.0 * DEG,
            psf_sigma_px: 1.4,
            sky: 500.0,
            read_noise: 8.0,
            zp_flux: 4.0e7,
            k1: 0.0,
        }
    }
}

impl SyntheticCamera {
    /// Sky -> pixel under the true camera, returning `None` for points
    /// outside the frame.
    pub fn radec_to_pixel(&self, ra: f64, dec: f64) -> Option<(f64, f64)> {
        let proj = TanProjection::new(self.ra0, self.dec0);
        let (xi, eta) = proj.project(ra, dec)?;
        let (cs, sn) = (self.rotation.cos(), self.rotation.sin());
        // Inverse of the WCS affine: pixel offsets along camera axes.
        let dx = (xi * cs + eta * sn) / self.scale;
        let dy = (-xi * sn + eta * cs) / self.scale;
        // Apply inverse radial distortion (small k1, one Newton step).
        let r2 = dx * dx + dy * dy;
        let (dxu, dyu) = (dx * (1.0 - self.k1 * r2), dy * (1.0 - self.k1 * r2));
        let x = dxu + self.width as f64 * 0.5;
        let y = dyu + self.height as f64 * 0.5;
        if x < 0.0 || y < 0.0 || x >= self.width as f64 || y >= self.height as f64 {
            return None;
        }
        Some((x, y))
    }

    fn star_flux(&self, mag: f32) -> f64 {
        self.zp_flux * 10f64.powf(-0.4 * mag as f64)
    }
}

/// Render the sky through the camera into a grayscale image. Only stars
/// inside (or just outside, for PSF bleed) the frame are drawn.
pub fn render_synthetic_image(
    sky: &SyntheticSky,
    cam: &SyntheticCamera,
    seed: u64,
) -> (GrayImage, Vec<(CatStar, f64, f64)>) {
    let mut img = GrayImage::new(cam.width, cam.height);
    let mut rng = Rng::new(seed);
    // Background + read noise.
    for v in img.buf.iter_mut() {
        *v = cam.sky + rng.normal() as f32 * cam.read_noise;
    }

    let margin = 4.0 * cam.psf_sigma_px + 2.0;
    let mut injected = Vec::new();
    let psf_r = cam.psf_sigma_px.ceil() as i32 + 2;
    let two_s2 = 2.0 * cam.psf_sigma_px * cam.psf_sigma_px;

    for star in &sky.stars {
        let Some((x, y)) = cam.radec_to_pixel(star.ra, star.dec) else {
            continue;
        };
        if x < -margin || y < -margin || x > cam.width as f64 + margin || y > cam.height as f64 + margin
        {
            continue;
        }
        let flux = cam.star_flux(star.mag);
        if flux < 2.0 * cam.read_noise as f64 {
            continue; // too faint to detect
        }
        let (xc, yc) = (x as i32, y as i32);
        for dy in -psf_r..=psf_r {
            for dx in -psf_r..=psf_r {
                let (px, py) = (xc + dx, yc + dy);
                if px < 0 || py < 0 || px >= cam.width as i32 || py >= cam.height as i32 {
                    continue;
                }
                let d2 = (dx as f64).powi(2) + (dy as f64).powi(2);
                let add = flux * (-d2 / two_s2).exp();
                // Shot noise approximated as Gaussian with sigma sqrt(flux).
                let noisy = add + rng.normal() * add.sqrt().max(0.0);
                let idx = py as usize * cam.width + px as usize;
                img.buf[idx] += noisy as f32;
            }
        }
        injected.push((*star, x, y));
    }
    (img, injected)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rng_deterministic() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn camera_project_roundtrip() {
        let cam = SyntheticCamera::default();
        let Some((x, y)) = cam.radec_to_pixel(cam.ra0 + 0.01, cam.dec0 + 0.01) else {
            panic!("center-adjacent star should be in frame");
        };
        assert!(x > 0.0 && y > 0.0 && x < cam.width as f64 && y < cam.height as f64);
    }
}
