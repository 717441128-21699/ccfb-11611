//! Spherical geometry and gnomonic (tangent-plane, WCS TAN) projection.
//!
//! All angles in radians unless stated. Coordinate convention:
//! right ascension in `[0, 2π)`, declination in `[-π/2, π/2]`.

#![allow(clippy::needless_range_loop)]

pub const TWO_PI: f64 = std::f64::consts::TAU;
pub const DEG: f64 = std::f64::consts::PI / 180.0;
pub const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);

#[inline]
pub fn wrap_two_pi(mut a: f64) -> f64 {
    a %= TWO_PI;
    if a < 0.0 {
        a += TWO_PI;
    }
    a
}

#[inline]
pub fn wrap_pi(a: f64) -> f64 {
    wrap_two_pi(a + std::f64::consts::PI) - std::f64::consts::PI
}

/// Unit 3-vector on the celestial sphere.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    #[inline]
    pub fn from_radec(ra: f64, dec: f64) -> Self {
        let (sd, cd) = (dec.sin(), dec.cos());
        Vec3 {
            x: cd * ra.cos(),
            y: cd * ra.sin(),
            z: sd,
        }
    }

    #[inline]
    pub fn to_radec(self) -> (f64, f64) {
        let dec = self.z.clamp(-1.0, 1.0).asin();
        let ra = wrap_two_pi(self.y.atan2(self.x));
        (ra, dec)
    }

    #[inline]
    pub fn dot(self, o: Vec3) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    #[inline]
    pub fn cross(self, o: Vec3) -> Vec3 {
        Vec3 {
            x: self.y * o.z - self.z * o.y,
            y: self.z * o.x - self.x * o.z,
            z: self.x * o.y - self.y * o.x,
        }
    }

    #[inline]
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }

    #[inline]
    pub fn normalize(self) -> Vec3 {
        let n = self.norm();
        Vec3 {
            x: self.x / n,
            y: self.y / n,
            z: self.z / n,
        }
    }

    /// Great-circle angular separation in radians, `[0, π]`.
    #[inline]
    pub fn separation(self, o: Vec3) -> f64 {
        // Stable for both small and large angles.
        self.cross(o).norm().atan2(self.dot(o)).abs()
    }
}

/// Gnomonic tangent-plane projection centered on `(ra0, dec0)`.
///
/// Projected coordinates `(xi, eta)` are radians on the tangent plane
/// (standard WCS intermediate world coordinates).
#[derive(Clone, Copy, Debug)]
pub struct TanProjection {
    pub ra0: f64,
    pub dec0: f64,
    pub sin_d0: f64,
    pub cos_d0: f64,
}

impl TanProjection {
    pub fn new(ra0: f64, dec0: f64) -> Self {
        TanProjection {
            ra0,
            dec0,
            sin_d0: dec0.sin(),
            cos_d0: dec0.cos(),
        }
    }

    /// Sky `(ra, dec)` -> tangent plane `(xi, eta)` radians.
    /// Returns `None` on the anti-meridian hemisphere (point cannot project).
    pub fn project(&self, ra: f64, dec: f64) -> Option<(f64, f64)> {
        let sd = dec.sin();
        let cd = dec.cos();
        let dra = ra - self.ra0;
        let c = self.sin_d0 * sd + self.cos_d0 * cd * dra.cos();
        if c <= 1e-12 {
            return None;
        }
        let xi = (cd * dra.sin()) / c;
        let eta = (self.cos_d0 * sd - self.sin_d0 * cd * dra.cos()) / c;
        Some((xi, eta))
    }

    /// Tangent plane `(xi, eta)` radians -> sky `(ra, dec)`.
    pub fn unproject(&self, xi: f64, eta: f64) -> (f64, f64) {
        let rho = (xi * xi + eta * eta).sqrt();
        if rho < 1e-15 {
            return (self.ra0, self.dec0);
        }
        let c = rho.atan();
        let (sc, cc) = (c.sin(), c.cos());
        let dec = (cc * self.sin_d0 + eta * sc * self.cos_d0 / rho).asin();
        let ra = self.ra0
            + (xi * sc)
                .atan2(rho * self.cos_d0 * cc - eta * self.sin_d0 * sc);
        (wrap_two_pi(ra), dec)
    }
}

/// Position angle of the direction from `(ra1,dec1)` to `(ra2,dec2)`
/// measured east of north, radians `[0, 2π)`.
pub fn position_angle(ra1: f64, dec1: f64, ra2: f64, dec2: f64) -> f64 {
    let dra = ra2 - ra1;
    let y = dra.sin() * dec2.cos();
    let x = dec1.cos() * dec2.sin() - dec1.sin() * dec2.cos() * dra.cos();
    wrap_two_pi(y.atan2(x))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn radec_roundtrip_vec3() {
        let (ra, dec) = (3.21, -0.45);
        let (ra2, dec2) = Vec3::from_radec(ra, dec).to_radec();
        assert!((ra - ra2).abs() < 1e-12);
        assert!((dec - dec2).abs() < 1e-12);
    }

    #[test]
    fn tan_project_roundtrip() {
        let proj = TanProjection::new(2.0, 0.3);
        // Points stay within one hemisphere of the projection center, as
        // required by the gnomonic projection.
        for &(ra, dec) in &[
            (2.01, 0.31),
            (1.95, 0.25),
            (2.3, 0.1),
            (2.2, 0.6),
            (1.8, -0.1),
        ] {
            let (xi, eta) = proj.project(ra, dec).unwrap();
            let (ra2, dec2) = proj.unproject(xi, eta);
            assert!((wrap_pi(ra - ra2)).abs() < 1e-11, "ra {}", ra);
            assert!((dec - dec2).abs() < 1e-11, "dec {}", ra);
        }
        // Points on the anti-meridian hemisphere cannot project.
        assert!(proj.project(2.0 + std::f64::consts::PI, 0.3).is_none());
    }

    #[test]
    fn separation_known() {
        // 1 degree on the equator.
        let a = Vec3::from_radec(0.0, 0.0);
        let b = Vec3::from_radec(DEG, 0.0);
        assert!((a.separation(b) - DEG).abs() < 1e-12);
    }

    #[test]
    fn position_angle_cardinal() {
        // Due east on equator -> PA = 90°.
        let pa = position_angle(0.0, 0.0, DEG, 0.0);
        assert!((pa - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
        // Due north -> PA = 0.
        let pa2 = position_angle(0.0, 0.0, 0.0, DEG);
        assert!(pa2.abs() < 1e-12 || (pa2 - TWO_PI).abs() < 1e-12);
    }
}
