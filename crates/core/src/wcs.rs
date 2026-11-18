//! WCS fitting: pixel → sky transform from star correspondences.
//!
//! Model (tangent plane around the solution center):
//! ```text
//! xi  = a0 + a1·x + a2·y + a3·r²·x' + a4·r²·y'
//! eta = b0 + b1·x + b2·y + b3·r²·x' + b4·r²·y'
//! ```
//! where `(x', y')` are pixel offsets from the image center and
//! `r² = x'² + y'²`. The linear 6-coefficient fit is the standard WCS SIP
//! zero-th order; the two radial terms absorb mild lens distortion and are
//! what gets the residual down toward the arcsecond regime on wide fields.
//!
//! The normal equations are solved by Cholesky decomposition (symmetric
//! positive-definite, small matrices) — no external linear algebra needed.

use crate::astro::{TanProjection, Vec3, ARCSEC};

/// A matched pair: image source index and catalog star identity.
#[derive(Clone, Copy, Debug)]
pub struct Correspondence {
    pub img_x: f64,
    pub img_y: f64,
    pub cat_ra: f64,
    pub cat_dec: f64,
}

/// Fitted world coordinate system.
#[derive(Clone, Debug)]
pub struct WcsSolution {
    /// Image-center sky coordinates.
    pub ra0: f64,
    pub dec0: f64,
    /// Linear affine coefficients (radians per pixel etc.).
    pub a: [f64; 3],
    pub b: [f64; 3],
    /// Radial distortion coefficients (radians per pixel^3 effectively).
    pub a_rad: f64,
    pub b_rad: f64,
    /// Image width/height the fit was built for (distortion reference).
    pub img_w: f64,
    pub img_h: f64,
    /// RMS residual in arcseconds over the inliers.
    pub rms_arcsec: f64,
    pub n_inliers: usize,
}

impl WcsSolution {
    /// Pixel `(x, y)` -> `(ra, dec)`.
    pub fn pixel_to_radec(&self, x: f64, y: f64) -> (f64, f64) {
        let (xi, eta) = self.pixel_to_tan(x, y);
        TanProjection::new(self.ra0, self.dec0).unproject(xi, eta)
    }

    /// Pixel -> tangent-plane radians relative to the solution center.
    pub fn pixel_to_tan(&self, x: f64, y: f64) -> (f64, f64) {
        let xc = x - self.img_w * 0.5;
        let yc = y - self.img_h * 0.5;
        let r2 = xc * xc + yc * yc;
        let xi = self.a[0] + self.a[1] * x + self.a[2] * y + self.a_rad * r2 * xc;
        let eta = self.b[0] + self.b[1] * x + self.b[2] * y + self.b_rad * r2 * yc;
        (xi, eta)
    }

    /// Field of view of the diagonal (radians).
    pub fn fov_diagonal_rad(&self) -> f64 {
        let (x1, y1) = self.pixel_to_tan(0.0, 0.0);
        let (x2, y2) = self.pixel_to_tan(self.img_w, self.img_h);
        ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt()
    }

    /// Horizontal and vertical field of view (radians). Each axis is the full
    /// tangent-plane *distance* between opposite edge midpoints (not a single
    /// coordinate difference), which makes the values independent of image
    /// rotation. The diagonal uses corner-to-corner distance analogously.
    pub fn fov_xy_rad(&self) -> (f64, f64) {
        let (w, h) = (self.img_w, self.img_h);
        let (xl, yl) = self.pixel_to_tan(0.0, h * 0.5);
        let (xr, yr) = self.pixel_to_tan(w, h * 0.5);
        let (xt, yt) = self.pixel_to_tan(w * 0.5, 0.0);
        let (xb, ybb) = self.pixel_to_tan(w * 0.5, h);
        let fov_x = ((xr - xl).powi(2) + (yr - yl).powi(2)).sqrt();
        let fov_y = ((xb - xt).powi(2) + (ybb - yt).powi(2)).sqrt();
        (fov_x, fov_y)
    }

    /// Plate scale in arcseconds per pixel, averaged over the two axes.
    pub fn plate_scale_arcsec(&self) -> f64 {
        let sx = (self.a[1].powi(2) + self.b[1].powi(2)).sqrt();
        let sy = (self.a[2].powi(2) + self.b[2].powi(2)).sqrt();
        0.5 * (sx + sy) / ARCSEC
    }

    /// Rotation angle of the image +x axis relative to celestial east,
    /// measured east of north (radians, `[0, 2π)`).
    pub fn rotation_east_of_north(&self) -> f64 {
        // Direction of increasing x on the tangent plane.
        let (dx, dy) = (self.a[1], self.b[1]);
        let theta = dy.atan2(dx);
        let mut t = theta;
        if t < 0.0 {
            t += 2.0 * std::f64::consts::PI;
        }
        t
    }

    /// Residual of one correspondence in arcseconds (catalog vs predicted).
    pub fn residual_arcsec(&self, c: &Correspondence) -> f64 {
        let (xi, eta) = self.pixel_to_tan(c.img_x, c.img_y);
        let (ra_p, dec_p) = TanProjection::new(self.ra0, self.dec0).unproject(xi, eta);
        let vp = Vec3::from_radec(ra_p, dec_p);
        let vc = Vec3::from_radec(c.cat_ra, c.cat_dec);
        vp.separation(vc) / ARCSEC
    }
}

#[derive(Clone, Copy, Debug)]
pub struct FitConfig {
    /// Include the two radial distortion terms.
    pub with_distortion: bool,
    /// RANSAC-style inlier cutoff in arcseconds.
    pub inlier_arcsec: f64,
}

impl Default for FitConfig {
    fn default() -> Self {
        FitConfig {
            with_distortion: true,
            inlier_arcsec: 8.0,
        }
    }
}

/// Cholesky solve `A·x = b` for symmetric positive-definite `A` (n×n).
fn cholesky_solve(a: &[f64], n: usize, b: &[f64]) -> Option<Vec<f64>> {
    let mut l = vec![0f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = a[i * n + j];
            for k in 0..j {
                s -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if s <= 1e-18 {
                    return None;
                }
                l[i * n + j] = s.sqrt();
            } else {
                l[i * n + j] = s / l[j * n + j];
            }
        }
    }
    // Forward substitution L·y = b.
    let mut y = vec![0f64; n];
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= l[i * n + k] * y[k];
        }
        y[i] = s / l[i * n + i];
    }
    // Back substitution Lᵀ·x = y.
    let mut x = vec![0f64; n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for k in (i + 1)..n {
            s -= l[k * n + i] * x[k];
        }
        x[i] = s / l[i * n + i];
    }
    Some(x)
}

/// Pure least-squares WCS around a fixed projection center. No rejection:
/// every projectable correspondence contributes. Returns `None` if the normal
/// equations are singular or too few points project.
fn lstsq_around(
    corrs: &[Correspondence],
    img_w: f64,
    img_h: f64,
    ra0: f64,
    dec0: f64,
    with_distortion: bool,
) -> Option<WcsSolution> {
    let proj = TanProjection::new(ra0, dec0);
    let ncoef = if with_distortion { 5 } else { 3 };

    let mut ata = vec![0f64; ncoef * ncoef];
    let mut atb_xi = vec![0f64; ncoef];
    let mut atb_eta = vec![0f64; ncoef];
    let mut used = 0usize;

    for c in corrs {
        let Some((xi, eta)) = proj.project(c.cat_ra, c.cat_dec) else {
            continue;
        };
        let xc = c.img_x - img_w * 0.5;
        let yc = c.img_y - img_h * 0.5;
        let r2 = xc * xc + yc * yc;
        let mut row = [0f64; 5];
        row[0] = 1.0;
        row[1] = c.img_x;
        row[2] = c.img_y;
        if with_distortion {
            row[3] = r2 * xc;
            row[4] = r2 * yc;
        }
        for i in 0..ncoef {
            for j in 0..ncoef {
                ata[i * ncoef + j] += row[i] * row[j];
            }
            atb_xi[i] += row[i] * xi;
            atb_eta[i] += row[i] * eta;
        }
        used += 1;
    }
    if used < ncoef + 1 {
        return None;
    }
    let x = cholesky_solve(&ata, ncoef, &atb_xi)?;
    let b = cholesky_solve(&ata, ncoef, &atb_eta)?;
    Some(WcsSolution {
        ra0,
        dec0,
        a: [x[0], x[1], x[2]],
        b: [b[0], b[1], b[2]],
        a_rad: if with_distortion { x[3] } else { 0.0 },
        b_rad: if with_distortion { b[3] } else { 0.0 },
        img_w,
        img_h,
        rms_arcsec: 0.0,
        n_inliers: used,
    })
}

/// Fit a WCS from correspondences, robust to gross outliers.
///
/// Pipeline:
/// 1. distortion-free affine least squares on the whole (contaminated) sample;
/// 2. rank every point by residual, keep those within `inlier_arcsec`, and
///    refit — the clean affine fit drives the rejection;
/// 3. full fit (with optional radial distortion) on the inliers;
/// 4. recenter on the predicted image center and refit once.
pub fn fit_wcs(
    corrs: &[Correspondence],
    img_w: f64,
    img_h: f64,
    cfg: &FitConfig,
) -> Option<WcsSolution> {
    if corrs.len() < 4 {
        return None;
    }
    // Initial center: robust median of catalog positions. The plain mean is
    // pulled toward gross outliers and biases the tangent-plane projection;
    // the median stays anchored on the good majority.
    let mut ras: Vec<f64> = corrs.iter().map(|c| c.cat_ra).collect();
    let mut decs: Vec<f64> = corrs.iter().map(|c| c.cat_dec).collect();
    ras.sort_by(|a, b| a.partial_cmp(b).unwrap());
    decs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let ra_med = ras[ras.len() / 2];
    let dec_med = decs[decs.len() / 2];
    let (ra0, dec0) = Vec3::from_radec(ra_med, dec_med).to_radec();

    // 1. RANSAC on a distortion-free affine model. Three non-collinear
    //    correspondences determine the 6 affine coefficients; sampling many
    //    triples makes a clean model overwhelmingly likely even when a large
    //    fraction of pairs are wrong (which is exactly the quad-vote regime).
    let ransac_cutoff = cfg.inlier_arcsec.max(8.0);
    let (inliers, _) = ransac_affine(corrs, img_w, img_h, ra0, dec0, ransac_cutoff)?;
    if inliers.len() < 4 {
        return None;
    }

    // 2. Full fit on the clean inlier set.
    let first = lstsq_around(&inliers, img_w, img_h, ra0, dec0, cfg.with_distortion)
        .or_else(|| lstsq_around(&inliers, img_w, img_h, ra0, dec0, false))?;

    // 3. Iteratively recenter the tangent plane on the predicted image
    //    center and refit. Each recentering reduces the gnomonic projection's
    //    nonlinear residual, converging in 2–3 passes to a stable center.
    let mut refined = first;
    let mut prev_center = refined.pixel_to_radec(img_w * 0.5, img_h * 0.5);
    for _ in 0..3 {
        let (ra_c, dec_c) = prev_center;
        let candidate =
            lstsq_around(&inliers, img_w, img_h, ra_c, dec_c, cfg.with_distortion)
                .or_else(|| lstsq_around(&inliers, img_w, img_h, ra_c, dec_c, false));
        let Some(next) = candidate else { break };
        let new_center = next.pixel_to_radec(img_w * 0.5, img_h * 0.5);
        let shift = Vec3::from_radec(new_center.0, new_center.1)
            .separation(Vec3::from_radec(prev_center.0, prev_center.1));
        refined = next;
        prev_center = new_center;
        if shift < 0.05 * ARCSEC {
            break; // converged to under a twentieth of an arcsecond
        }
    }

    // Report RMS over the inlier set.
    let resid: Vec<f64> = inliers.iter().map(|c| refined.residual_arcsec(c)).collect();
    let mut out = refined;
    out.n_inliers = inliers.len();
    out.rms_arcsec = (resid.iter().map(|r| r * r).sum::<f64>() / resid.len() as f64).sqrt();
    Some(out)
}

/// RANSAC over the 6-parameter affine tangent-plane model. Returns the
/// largest consensus set (all points within `cutoff_arcsec` of the best
/// triple's model) and the fitted solution.
fn ransac_affine(
    corrs: &[Correspondence],
    img_w: f64,
    img_h: f64,
    ra0: f64,
    dec0: f64,
    cutoff_arcsec: f64,
) -> Option<(Vec<Correspondence>, WcsSolution)> {
    let n = corrs.len();
    if n < 3 {
        return None;
    }
    // Deterministic stratified sampling: spread triples across the data
    // rather than relying on an RNG, so behavior is reproducible.
    let max_iters = 200.min(n * n.saturating_sub(1) * n.saturating_sub(2) / 6);
    let mut best_set: Vec<Correspondence> = Vec::new();
    let mut best_model: Option<WcsSolution> = None;

    let mut trial = |idx: [usize; 3]| -> Option<()> {
        let sample: Vec<Correspondence> = idx.iter().map(|&i| corrs[i]).collect();
        // An affine model needs 3 correspondences for 6 unknowns; add the
        // median-center correspondence as a soft fourth anchor by fitting on
        // the triple directly (3 points exactly determine xi/eta affine).
        let model = fit_affine_exact(&sample, img_w, img_h, ra0, dec0)?;
        let set: Vec<Correspondence> = corrs
            .iter()
            .copied()
            .filter(|c| model.residual_arcsec(c) <= cutoff_arcsec)
            .collect();
        if set.len() > best_set.len() {
            best_set = set;
            best_model = Some(model);
        }
        Some(())
    };

    // Cover the data with evenly spaced triples plus adjacent combinations.
    let mut iters = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                if iters >= max_iters {
                    break;
                }
                let _ = trial([i, j, k]);
                iters += 1;
            }
            if iters >= max_iters {
                break;
            }
        }
        if iters >= max_iters {
            break;
        }
    }

    let model = best_model?;
    if best_set.len() < 3 {
        return None;
    }
    // Final least-squares refit over the full consensus set.
    let refit = lstsq_around(&best_set, img_w, img_h, ra0, dec0, false).unwrap_or(model);
    Some((best_set, refit))
}

/// Exact affine fit through exactly three correspondences (6 equations,
/// 6 unknowns). Used as the RANSAC minimal-sample model.
fn fit_affine_exact(
    sample: &[Correspondence],
    img_w: f64,
    img_h: f64,
    ra0: f64,
    dec0: f64,
) -> Option<WcsSolution> {
    if sample.len() != 3 {
        return None;
    }
    let proj = TanProjection::new(ra0, dec0);
    // Build two 3x3 systems A·coef = rhs from the three points.
    let mut a = [[0f64; 3]; 3];
    let mut bx = [0f64; 3];
    let mut by = [0f64; 3];
    for (i, c) in sample.iter().enumerate() {
        let (xi, eta) = proj.project(c.cat_ra, c.cat_dec)?;
        a[i] = [1.0, c.img_x, c.img_y];
        bx[i] = xi;
        by[i] = eta;
    }
    let x = solve_3x3(&a, &bx)?;
    let b = solve_3x3(&a, &by)?;
    Some(WcsSolution {
        ra0,
        dec0,
        a: [x[0], x[1], x[2]],
        b: [b[0], b[1], b[2]],
        a_rad: 0.0,
        b_rad: 0.0,
        img_w,
        img_h,
        rms_arcsec: 0.0,
        n_inliers: 3,
    })
}

/// Solve a 3×3 linear system by Cramer's rule.
fn solve_3x3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
    let det = |m: &[[f64; 3]; 3]| -> f64 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(a);
    if d.abs() < 1e-18 {
        return None;
    }
    let mut out = [0f64; 3];
    for col in 0..3 {
        let mut rep = *a;
        for row in 0..3 {
            rep[row][col] = rhs[row];
        }
        out[col] = det(&rep) / d;
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::DEG;

    fn synth_wcs() -> WcsSolution {
        // 4000x3000 image, 2 arcsec/pixel, 35° rotation, center on M31-ish.
        // Constant terms are chosen so the pixel center (2000, 1500) maps to
        // tangent-plane origin, i.e. (xi, eta) = (0, 0) at image center.
        let scale = 2.0 * ARCSEC;
        let theta = 35.0f64.to_radians();
        let (w, h) = (4000.0f64, 3000.0f64);
        let (a1, a2) = (scale * theta.cos(), -scale * theta.sin());
        let (b1, b2) = (scale * theta.sin(), scale * theta.cos());
        WcsSolution {
            ra0: 10.68 * DEG,
            dec0: 41.27 * DEG,
            a: [-(a1 * w * 0.5 + a2 * h * 0.5), a1, a2],
            b: [-(b1 * w * 0.5 + b2 * h * 0.5), b1, b2],
            a_rad: 0.0,
            b_rad: 0.0,
            img_w: w,
            img_h: h,
            rms_arcsec: 0.0,
            n_inliers: 0,
        }
    }

    #[test]
    fn recovers_known_transform() {
        let truth = synth_wcs();
        let mut corrs = Vec::new();
        for i in 0..40 {
            let x = 200.0 + (i % 8) as f64 * 450.0;
            let y = 150.0 + (i / 8) as f64 * 380.0;
            let (ra, dec) = truth.pixel_to_radec(x, y);
            corrs.push(Correspondence {
                img_x: x,
                img_y: y,
                cat_ra: ra,
                cat_dec: dec,
            });
        }
        let fit = fit_wcs(
            &corrs,
            4000.0,
            3000.0,
            &FitConfig {
                with_distortion: false,
                inlier_arcsec: 1.0,
            },
        )
        .unwrap();
        assert!((fit.plate_scale_arcsec() - 2.0).abs() < 1e-6);
        let mut rot = fit.rotation_east_of_north() - truth.rotation_east_of_north();
        rot = (rot + std::f64::consts::PI) % (2.0 * std::f64::consts::PI)
            - std::f64::consts::PI;
        assert!(rot.abs() < 1e-5, "rotation error {rot}");
        // Center accuracy well under an arcsecond.
        let vc = Vec3::from_radec(fit.ra0, fit.dec0);
        let vt = Vec3::from_radec(truth.ra0, truth.dec0);
        assert!(vc.separation(vt) / ARCSEC < 0.05, "center off");
        assert!(fit.rms_arcsec < 0.01);
    }

    #[test]
    fn rejects_outliers() {
        let truth = synth_wcs();
        let mut corrs = Vec::new();
        for i in 0..30 {
            let x = 300.0 + (i % 6) as f64 * 600.0;
            let y = 200.0 + (i / 6) as f64 * 500.0;
            let (ra, dec) = truth.pixel_to_radec(x, y);
            corrs.push(Correspondence {
                img_x: x,
                img_y: y,
                cat_ra: ra,
                cat_dec: dec,
            });
        }
        // Three grossly wrong pairs, but kept inside the projection
        // hemisphere so tangent projection remains defined.
        for i in 0..3 {
            corrs.push(Correspondence {
                img_x: 500.0 + i as f64 * 100.0,
                img_y: 500.0,
                cat_ra: corrs[0].cat_ra + 0.05,
                cat_dec: corrs[0].cat_dec + 0.04,
            });
        }
        // Default 8-arcsec cutoff rejects the ~2.9° gross outliers while
        // keeping all noise-free synthetic inliers.
        let cfg = FitConfig::default();
        let fit = fit_wcs(&corrs, 4000.0, 3000.0, &cfg).unwrap();
        assert!(fit.n_inliers >= 29, "inliers = {}", fit.n_inliers);
        assert!(fit.rms_arcsec < 0.5, "rms = {}", fit.rms_arcsec);
    }
}
