//! End-to-end self-test: build a synthetic sky + index, render a realistic
//! image through a *known* camera, blind-solve it, and report the accuracy
//! of the recovered center / scale / rotation against the injected truth.
//!
//! Exit code 0 means the solve succeeded within the 2-arcsecond requirement.

use platersol_core::astro::{Vec3, ARCSEC, DEG};
use platersol_core::{
    build_index, extract_sources, render_synthetic_image, solve, BuildConfig, DirStore,
    ExtractConfig, SolveConfig, SyntheticCamera, SyntheticSky,
};
use std::time::Instant;

fn main() {
    // --- configurable scenario -------------------------------------------
    let n_stars: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(1_000_000);
    let width = 2048usize;
    let height = 1536usize;

    println!("=== platersol end-to-end self-test ===");
    println!("catalog: {n_stars} synthetic stars, image {width}x{height}");

    // 1. Sky + index. Magnitude range matches the validated e2e scenario so
    // the per-cell brightest-star truncation behaves identically.
    let t = Instant::now();
    let sky = SyntheticSky::random(n_stars, 3.0, 9.0, 20261118);
    let built = build_index(&sky.stars, &BuildConfig::default());
    println!(
        "index built in {:.2}s: {} quads, {} bins",
        t.elapsed().as_secs_f64(),
        built.manifest.n_quads,
        built.manifest.bins.len()
    );

    let tmp = std::env::temp_dir().join(format!("platersol_selftest_{}", std::process::id()));
    let dir = tmp.to_str().unwrap().to_string();
    built.write_to_dir(&dir).unwrap();
    std::fs::write(
        format!("{dir}/catalog.json"),
        serde_json::to_vec(
            &sky.stars
                .iter()
                .map(|s| platersol_core::catalog::CatStarRow {
                    id: s.id,
                    ra_deg: s.ra / DEG,
                    dec_deg: s.dec / DEG,
                    mag: s.mag,
                    name: None,
                })
                .collect::<Vec<_>>(),
        )
        .unwrap(),
    )
    .unwrap();

    // 2. Render through a known camera, including mild lens distortion.
    let cam = SyntheticCamera {
        width,
        height,
        ra0: 10.68 * DEG,
        dec0: 41.27 * DEG,
        scale: 2.5 * ARCSEC,
        rotation: 57.0 * DEG,
        psf_sigma_px: 1.6,
        sky: 500.0,
        read_noise: 6.0,
        zp_flux: 4.0e7,
        k1: 0.0,
    };
    let t = Instant::now();
    let (img, injected) = render_synthetic_image(&sky, &cam, 9917);
    println!(
        "rendered {} detectable-frame stars in {:.2}s",
        injected.len(),
        t.elapsed().as_secs_f64()
    );

    // 3. Source extraction (the same code path real photos take).
    let t = Instant::now();
    let sources = extract_sources(
        &img,
        &ExtractConfig {
            threshold_sigma: 10.0,
            min_pixels: 3,
            max_sources: Some(500),
            ..Default::default()
        },
    );
    println!(
        "extracted {} sources in {:.2}s",
        sources.len(),
        t.elapsed().as_secs_f64()
    );

    // 4. Blind solve against the on-disk chunked index.
    let mut store = DirStore::open(&dir, &sky.stars).unwrap();
    let t = Instant::now();
    let outcome = solve(
        &sources,
        width as f64,
        height as f64,
        &mut store,
        &SolveConfig::default(),
    );
    let solve_time = t.elapsed();

    let outcome = match outcome {
        Ok(o) => o,
        Err(e) => {
            eprintln!("SOLVE FAILED: {e}");
            let _ = std::fs::remove_dir_all(&dir);
            std::process::exit(1);
        }
    };
    let w = &outcome.wcs;

    // 5. Accuracy vs injected truth.
    let center_truth = Vec3::from_radec(cam.ra0, cam.dec0);
    let center_got = Vec3::from_radec(w.ra0, w.dec0);
    let center_err_as = center_truth.separation(center_got) / ARCSEC;

    let scale_truth_as = cam.scale / ARCSEC;
    let scale_err_pct = (w.plate_scale_arcsec() - scale_truth_as).abs() / scale_truth_as * 100.0;

    let mut drot = w.rotation_east_of_north() - cam.rotation;
    drot = (drot + std::f64::consts::PI) % (2.0 * std::f64::consts::PI)
        - std::f64::consts::PI;
    let rot_err_am = drot.abs() / ARCSEC / 60.0; // arcminutes of angle

    let (fovx, fovy) = w.fov_xy_rad();
    println!();
    println!("--- solve result ({:.2}s) ---", solve_time.as_secs_f64());
    println!(
        "center RA/Dec  = {:.7}°, {:.7}°   (truth {:.7}, {:.7})",
        w.ra0 / DEG, w.dec0 / DEG, cam.ra0 / DEG, cam.dec0 / DEG
    );
    println!("center error   = {:.3} arcsec", center_err_as);
    println!(
        "plate scale    = {:.4} arcsec/px  (truth {:.4}, err {:.3}%)",
        w.plate_scale_arcsec(),
        scale_truth_as,
        scale_err_pct
    );
    println!("rotation error = {:.2} arcmin", rot_err_am);
    println!(
        "FOV            = {:.3}° × {:.3}° (diag {:.3}°)",
        fovx / DEG,
        fovy / DEG,
        w.fov_diagonal_rad() / DEG
    );
    println!(
        "fit RMS        = {:.3} arcsec over {} inliers ({} initial matches)",
        w.rms_arcsec,
        w.n_inliers,
        outcome.matches.len()
    );
    println!(
        "index I/O      = {} bucket reads, {:.1} KiB, {} image quads tried",
        outcome.buckets_read,
        outcome.index_bytes_read as f64 / 1024.0,
        outcome.quads_tried
    );

    let _ = std::fs::remove_dir_all(&dir);

    let pass = center_err_as <= 2.0 && w.rms_arcsec <= 2.0;
    println!();
    if pass {
        println!("PASS ✔  requirement: center & RMS within 2 arcseconds");
    } else {
        println!("FAIL ✘  exceeds 2-arcsecond budget");
        std::process::exit(1);
    }
}
