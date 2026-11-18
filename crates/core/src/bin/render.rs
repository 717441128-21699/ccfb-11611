//! Render a synthetic deep-sky PGM through a known camera, for offline demos
//! and for exercising the ps-solve CLI without a real photograph.
//!
//! Usage:
//!   ps-render <index_dir> <out.pgm> [ra_deg] [dec_deg] [rotation_deg]

use platersol_core::astro::{DEG, ARCSEC};
use platersol_core::catalog::{CatStar, CatStarRow};
use platersol_core::astro::Vec3;
use platersol_core::{render_synthetic_image, SyntheticCamera, SyntheticSky};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: ps-render <index_dir> <out.pgm> [ra_deg] [dec_deg] [rotation_deg]");
        std::process::exit(2);
    }
    let index_dir = &args[1];
    let out_path = &args[2];
    let ra: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(10.68);
    let dec: f64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(41.27);
    let rot: f64 = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(57.0);

    // Rebuild the identical synthetic catalog the index was built from.
    let rows: Vec<CatStarRow> =
        serde_json::from_slice(&std::fs::read(format!("{index_dir}/catalog.json")).unwrap())
            .unwrap();
    let stars: Vec<CatStar> = rows
        .iter()
        .map(|r| {
            let (ra, dec) = (r.ra_deg * DEG, r.dec_deg * DEG);
            CatStar { id: r.id, ra, dec, mag: r.mag, v: Vec3::from_radec(ra, dec) }
        })
        .collect();
    let sky = SyntheticSky { stars };

    let cam = SyntheticCamera {
        width: 2048,
        height: 1536,
        ra0: ra * DEG,
        dec0: dec * DEG,
        scale: 2.5 * ARCSEC,
        rotation: rot * DEG,
        psf_sigma_px: 1.6,
        sky: 500.0,
        read_noise: 6.0,
        zp_flux: 4.0e7,
        k1: 0.0,
    };

    let (img, injected) = render_synthetic_image(&sky, &cam, 9917);
    img.write_pgm(out_path).unwrap();
    println!(
        "wrote {out_path} ({}x{}, {} injected stars, truth center RA={} Dec={})",
        img.width, img.height, injected.len(), ra, dec
    );
}
