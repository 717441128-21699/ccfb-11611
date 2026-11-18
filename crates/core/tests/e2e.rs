//! End-to-end pipeline test on synthetic data (in-memory index store).

use platersol_core::astro::{Vec3, ARCSEC, DEG};
use platersol_core::index::{build_index, IndexEntry, IndexManifest};
use platersol_core::solver::{solve, IndexStore, SolveConfig};
use platersol_core::{
    extract_sources, render_synthetic_image, BuildConfig, CatStar, ExtractConfig,
    SyntheticCamera, SyntheticSky,
};
use std::collections::HashMap;

/// In-memory store holding every built bucket directly.
struct MemStore {
    manifest: IndexManifest,
    catalog: HashMap<u32, CatStar>,
    buckets: HashMap<(usize, usize), Vec<u8>>,
    known: Vec<u32>,
}

impl MemStore {
    fn from_built(
        manifest: IndexManifest,
        per_bin: Vec<Vec<Vec<IndexEntry>>>,
        catalog: &[CatStar],
    ) -> Self {
        let mut buckets = HashMap::new();
        for (bin, bks) in per_bin.into_iter().enumerate() {
            for (bucket, entries) in bks.into_iter().enumerate() {
                if entries.is_empty() {
                    continue;
                }
                let mut raw = Vec::with_capacity(entries.len() * 24);
                for e in entries {
                    e.write_to(&mut raw).unwrap();
                }
                buckets.insert((bin, bucket), raw);
            }
        }
        let catalog: HashMap<u32, CatStar> =
            catalog.iter().map(|s| (s.id, *s)).collect();
        let known = catalog.keys().copied().collect();
        MemStore {
            manifest,
            catalog,
            buckets,
            known,
        }
    }
}

impl IndexStore for MemStore {
    fn manifest(&self) -> &IndexManifest {
        &self.manifest
    }
    fn load_bucket(&mut self, bin: usize, bucket: usize) -> std::io::Result<Vec<u8>> {
        self.buckets
            .get(&(bin, bucket))
            .cloned()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::NotFound, "empty bucket"))
    }
    fn star_by_id(&self, id: u32) -> Option<CatStar> {
        self.catalog.get(&id).copied()
    }
    fn known_ids(&self) -> Vec<u32> {
        self.known.clone()
    }
    fn bytes_fetched(&self) -> u64 {
        0
    }
}

#[test]
fn synthetic_blind_solve_end_to_end() {
    // Dense sky: a million catalog stars yields several dozen stars even in
    // a narrow ~1.4° deep-sky frame, matching real source counts.
    let sky = SyntheticSky::random(1_000_000, 3.0, 9.0, 20261118);
    let cam = SyntheticCamera {
        width: 2048,
        height: 1536,
        ra0: 10.68 * DEG,
        dec0: 41.27 * DEG,
        // 2.5 arcsec/px -> ~1.4° x 1.1° frame, quad A-B within index range.
        scale: 2.5 * ARCSEC,
        rotation: 57.0 * DEG,
        psf_sigma_px: 1.6,
        sky: 500.0,
        read_noise: 6.0,
        zp_flux: 4.0e7,
        k1: 0.0,
    };
    let (img, injected) = render_synthetic_image(&sky, &cam, 9917);
    eprintln!("injected stars in frame: {}", injected.len());
    assert!(injected.len() >= 12, "injected {}", injected.len());

    let sources = extract_sources(
        &img,
        &ExtractConfig {
            threshold_sigma: 10.0,
            min_pixels: 3,
            max_sources: Some(500),
            ..Default::default()
        },
    );
    // Extraction quality is judged by how many *real* injected stars are
    // recovered, not by suppressing every noise peak (the solver tolerates
    // noise sources; the brightest-60 quad selection ignores the faint ones).

    // Build the geometric index.
    let t0 = std::time::Instant::now();
    let built = build_index(&sky.stars, &BuildConfig::default());
    eprintln!("index build: {:.2}s, {} quads", t0.elapsed().as_secs_f64(), built.manifest.n_quads);
    let manifest = built.manifest.clone();
    let store = MemStore::from_built(manifest, built.per_bin, &sky.stars);

    // How many of the extracted sources actually coincide with an injected
    // star within 2 pixels?
    let mut matched_real = 0;
    for s in &sources {
        let near = injected
            .iter()
            .any(|(_, x, y)| (s.x - x).hypot(s.y - y) < 2.0);
        if near {
            matched_real += 1;
        }
    }
    eprintln!(
        "extracted={} injected={} real_recovered={}",
        sources.len(),
        injected.len(),
        matched_real
    );
    assert!(matched_real >= 20, "only {matched_real} real stars recovered");

    let mut store = store;
    let t1 = std::time::Instant::now();
    let outcome = match solve(
        &sources,
        cam.width as f64,
        cam.height as f64,
        &mut store,
        &SolveConfig::default(),
    ) {
        Ok(o) => o,
        Err(e) => panic!("blind solve failed: {e}"),
    };
    eprintln!("solve: {:.2}s", t1.elapsed().as_secs_f64());

    let w = &outcome.wcs;
    let center_err = Vec3::from_radec(w.ra0, w.dec0)
        .separation(Vec3::from_radec(cam.ra0, cam.dec0))
        / ARCSEC;
    eprintln!(
        "SOLVED center_err={:.3} arcsec rms={:.3} arcsec inliers={} matches={} quads={}",
        center_err, w.rms_arcsec, w.n_inliers, outcome.matches.len(), outcome.quads_tried
    );
    assert!(center_err <= 2.0, "center error {center_err} arcsec");
    assert!(w.rms_arcsec <= 2.0, "rms {} arcsec", w.rms_arcsec);
}
