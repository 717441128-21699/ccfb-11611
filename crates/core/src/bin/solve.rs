//! Plate solver CLI for real grayscale images.
//!
//! Usage:
//!   ps-solve <image.pgm> <index_dir> [threshold_sigma]
//!
//! Loads a binary PGM (P5, 8/16-bit), extracts sources, blind-solves against
//! a prebuilt chunked index, and prints the WCS plus an accuracy-independent
//! self-description. PGM keeps this tool codec-free; the browser path decodes
//! JPEG/PNG in TypeScript and feeds f32 pixels to the same WASM engine.

use platersol_core::astro::{DEG, Vec3};
use platersol_core::catalog::{CatStar, CatStarRow};
use platersol_core::{
    extract_sources, solve, DirStore, ExtractConfig, GrayImage, SolveConfig,
};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: ps-solve <image.pgm> <index_dir> [threshold_sigma]");
        std::process::exit(2);
    }
    let img_path = &args[1];
    let index_dir = &args[2];
    let threshold: f32 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(10.0);

    // 1. Load image.
    let t = Instant::now();
    let img = GrayImage::read_pgm(img_path).unwrap_or_else(|e| {
        eprintln!("cannot read image {img_path}: {e}");
        std::process::exit(1);
    });
    println!(
        "image {}x{} loaded in {:.2}s",
        img.width,
        img.height,
        t.elapsed().as_secs_f64()
    );

    // 2. Load the companion catalog the index was built from.
    let cat_path = format!("{index_dir}/catalog.json");
    let rows: Vec<CatStarRow> = serde_json::from_slice(
        &std::fs::read(&cat_path).unwrap_or_else(|e| {
            eprintln!("cannot read {cat_path}: {e}");
            std::process::exit(1);
        }),
    )
    .unwrap_or_else(|e| {
        eprintln!("bad catalog JSON: {e}");
        std::process::exit(1);
    });
    let catalog: Vec<CatStar> = rows
        .iter()
        .map(|r| {
            let ra = r.ra_deg * DEG;
            let dec = r.dec_deg * DEG;
            CatStar {
                id: r.id,
                ra,
                dec,
                mag: r.mag,
                v: Vec3::from_radec(ra, dec),
            }
        })
        .collect();

    // 3. Extract sources.
    let t = Instant::now();
    let sources = extract_sources(
        &img,
        &ExtractConfig {
            threshold_sigma: threshold,
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
    if sources.len() < 4 {
        eprintln!("too few sources to solve");
        std::process::exit(1);
    }

    // 4. Blind solve via the directory-backed chunk store.
    let mut store = DirStore::open(index_dir, &catalog).unwrap_or_else(|e| {
        eprintln!("cannot open index at {index_dir}: {e}");
        std::process::exit(1);
    });
    let t = Instant::now();
    let outcome = solve(
        &sources,
        img.width as f64,
        img.height as f64,
        &mut store,
        &SolveConfig::default(),
    );
    let elapsed = t.elapsed();

    match outcome {
        Ok(o) => {
            let w = &o.wcs;
            let (fovx, fovy) = w.fov_xy_rad();
            println!();
            println!("=== SOLVED in {:.2}s ===", elapsed.as_secs_f64());
            println!("center RA  = {:.7} deg", w.ra0 / DEG);
            println!("center Dec = {:.7} deg", w.dec0 / DEG);
            println!("FOV        = {:.4} x {:.4} deg", fovx / DEG, fovy / DEG);
            println!("rotation   = {:.4} deg east of north", w.rotation_east_of_north() / DEG);
            println!("plate scale= {:.4} arcsec/px", w.plate_scale_arcsec());
            println!("fit RMS    = {:.3} arcsec over {} inliers", w.rms_arcsec, w.n_inliers);
            println!(
                "index I/O  = {} bucket reads, {:.1} KiB, {} quads",
                o.buckets_read,
                o.index_bytes_read as f64 / 1024.0,
                o.quads_tried
            );
        }
        Err(e) => {
            eprintln!("SOLVE FAILED: {e}");
            std::process::exit(1);
        }
    }
}
