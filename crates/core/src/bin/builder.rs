//! Index builder CLI.
//!
//! Modes:
//!   ps-builder synth <out_dir> [n_stars] [seed]
//!       Build an index from a deterministic synthetic uniform catalog.
//!   ps-builder catalog <catalog.json> <out_dir>
//!       Build from a JSON array of CatStarRow ({id, ra_deg, dec_deg, mag}).
//!
//! The catalog used for solving must also be available at runtime, so the
//! builder writes `catalog.json` next to the index chunks.

use platersol_core::astro::{DEG, Vec3};
use platersol_core::catalog::{CatStar, CatStarRow};
use platersol_core::synthetic::SyntheticSky;
use platersol_core::{build_index, BuildConfig};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage:\n  ps-builder synth <out_dir> [n_stars] [seed]\n  ps-builder catalog <catalog.json> <out_dir>"
        );
        std::process::exit(2);
    }
    let mode = &args[1];
    let (catalog, out_dir) = match mode.as_str() {
        "synth" => {
            let out_dir = args[2].clone();
            let n: u32 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(20_000);
            let seed: u64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(20261118);
            let sky = SyntheticSky::random(n, 3.0, 10.0, seed);
            (sky.stars, out_dir)
        }
        "catalog" => {
            let cat_path = args[2].clone();
            let out_dir = args
                .get(3)
                .cloned()
                .expect("catalog mode needs <catalog.json> <out_dir>");
            let rows: Vec<CatStarRow> =
                serde_json::from_slice(&std::fs::read(&cat_path).unwrap_or_else(|e| {
                    eprintln!("cannot read {cat_path}: {e}");
                    std::process::exit(1);
                }))
                .unwrap_or_else(|e| {
                    eprintln!("bad catalog JSON: {e}");
                    std::process::exit(1);
                });
            let stars = rows
                .into_iter()
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
            (stars, out_dir)
        }
        other => {
            eprintln!("unknown mode {other}");
            std::process::exit(2);
        }
    };

    eprintln!(
        "building index: {} catalog stars, {} scale bins/cell cfg",
        catalog.len(),
        BuildConfig::default().bins_per_decade
    );
    let t0 = std::time::Instant::now();
    let cfg = BuildConfig::default();
    let mut built = build_index(&catalog, &cfg);
    built.manifest.catalog = mode.clone();
    eprintln!(
        "build done in {:.2}s: {} quads across {} bins",
        t0.elapsed().as_secs_f64(),
        built.manifest.n_quads,
        built.manifest.bins.len()
    );

    std::fs::create_dir_all(&out_dir).expect("create out dir");
    built.write_to_dir(&out_dir).expect("write index");

    // Persist the catalog rows required at solve time.
    let rows: Vec<CatStarRow> = catalog
        .iter()
        .map(|s| CatStarRow {
            id: s.id,
            ra_deg: s.ra / DEG,
            dec_deg: s.dec / DEG,
            mag: s.mag,
            name: None,
        })
        .collect();
    std::fs::write(
        format!("{out_dir}/catalog.json"),
        serde_json::to_vec(&rows).unwrap(),
    )
    .expect("write catalog");

    // Per-bin size report.
    let mut total_bytes = 0u64;
    for bin in 0..built.manifest.bins.len() {
        let bytes = built.write_bin(bin).len() as u64;
        total_bytes += bytes;
        let n: u64 = built.per_bin[bin].iter().map(|b| b.len() as u64).sum();
        if n > 0 {
            eprintln!(
                "  bin {:2}  ab={:7.3}–{:7.3} deg  quads={:7}  {:8.1} KiB",
                bin,
                built.manifest.bins[bin].lo_rad / DEG,
                built.manifest.bins[bin].hi_rad / DEG,
                n,
                bytes as f64 / 1024.0
            );
        }
    }
    eprintln!(
        "wrote index to {out_dir} ({:.1} MiB total, manifest + {} bins + catalog)",
        total_bytes as f64 / 1048576.0,
        built.manifest.bins.len()
    );
}
