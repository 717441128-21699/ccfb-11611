//! Verify DirStore returns byte-identical bucket content to an in-memory
//! reference, since the browser ChunkStore uses the same offset/slice logic.

use platersol_core::astro::{ARCSEC, DEG};
use platersol_core::index::{build_index, IndexEntry, IndexManifest};
use platersol_core::solver::{DirStore, IndexStore};
use platersol_core::{render_synthetic_image, BuildConfig, SyntheticCamera, SyntheticSky};
use std::collections::HashMap;

#[test]
fn dirstore_buckets_match_reference() {
    // Small catalog keeps this fast.
    let sky = SyntheticSky::random(50_000, 3.0, 9.0, 20261118);
    let built = build_index(&sky.stars, &BuildConfig::default());
    let manifest: IndexManifest = built.manifest.clone();

    // Reference: every non-empty bucket's raw bytes from the in-memory build.
    let mut reference: HashMap<(usize, usize), Vec<u8>> = HashMap::new();
    for (bin, buckets) in built.per_bin.iter().enumerate() {
        for (bucket, entries) in buckets.iter().enumerate() {
            if entries.is_empty() {
                continue;
            }
            let mut raw = Vec::new();
            for e in entries {
                e.write_to(&mut raw).unwrap();
            }
            reference.insert((bin, bucket), raw);
        }
    }

    // Write to disk and reopen through DirStore.
    let tmp = std::env::temp_dir().join(format!("platersol_dirstore_{}", std::process::id()));
    let dir = tmp.to_str().unwrap().to_string();
    built.write_to_dir(&dir).unwrap();
    let mut store = DirStore::open(&dir, &sky.stars).unwrap();

    let mut checked = 0;
    for ((bin, bucket), expect) in &reference {
        let got = store.load_bucket(*bin, *bucket).expect("bucket read");
        assert_eq!(&got, expect, "mismatch at bin {bin} bucket {bucket}");
        // Parse back to confirm entries survive.
        assert_eq!(got.len() % 24, 0);
        checked += 1;
    }
    assert!(checked > 100, "only {checked} buckets compared");

    // And an empty bucket must return an empty (not erroring) slice.
    let empties: Vec<(usize, usize)> = (0..manifest.bins.len())
        .flat_map(|b| (0..1024).map(move |k| (b, k)))
        .filter(|(b, k)| !reference.contains_key(&(*b, *k)))
        .take(5)
        .collect();
    for (b, k) in empties {
        let got = store.load_bucket(b, k).unwrap();
        assert!(got.is_empty());
    }

    let _ = std::fs::remove_dir_all(&dir);
    let _ = (ARCSEC, DEG);
}

#[test]
fn dirstore_full_solve_matches_memory() {
    // Same dense-sky scenario as the validated e2e test, but driven through
    // the on-disk DirStore (the path the browser ChunkStore mirrors).
    let sky = SyntheticSky::random(1_000_000, 3.0, 9.0, 20261118);
    let cam = SyntheticCamera {
        width: 2048,
        height: 1536,
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
    let (img, injected) = render_synthetic_image(&sky, &cam, 9917);
    assert!(injected.len() >= 8, "injected {}", injected.len());

    let built = build_index(&sky.stars, &BuildConfig::default());
    let tmp = std::env::temp_dir().join(format!("platersol_dirstore_solve_{}", std::process::id()));
    let dir = tmp.to_str().unwrap().to_string();
    built.write_to_dir(&dir).unwrap();

    let sources = platersol_core::extract_sources(
        &img,
        &platersol_core::ExtractConfig {
            threshold_sigma: 10.0,
            min_pixels: 3,
            max_sources: Some(500),
            ..Default::default()
        },
    );
    eprintln!("injected={} extracted={}", injected.len(), sources.len());

    let mut store = DirStore::open(&dir, &sky.stars).unwrap();
    let result = platersol_core::solve(
        &sources,
        cam.width as f64,
        cam.height as f64,
        &mut store,
        &platersol_core::SolveConfig::default(),
    );
    match result {
        Ok(o) => {
            let err = platersol_core::Vec3::from_radec(o.wcs.ra0, o.wcs.dec0)
                .separation(platersol_core::Vec3::from_radec(cam.ra0, cam.dec0))
                / ARCSEC;
            eprintln!("DirStore solve OK: center_err={:.2} arcsec rms={:.2}", err, o.wcs.rms_arcsec);
            assert!(err <= 2.0);
        }
        Err(e) => panic!("DirStore solve failed: {e}"),
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[allow(dead_code)]
fn _entry_type_used(_e: IndexEntry) {}
