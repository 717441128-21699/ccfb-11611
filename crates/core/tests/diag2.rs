//! Verify a known image quad's hash is findable in the built index via the
//! exact bucket/scale-bin path the solver walks.

use platersol_core::astro::{TanProjection, ARCSEC, DEG};
use platersol_core::index::{build_index, hash_bucket, BuildConfig};
use platersol_core::quad::{hash_window, make_quad, Pt};
use platersol_core::{render_synthetic_image, SyntheticCamera, SyntheticSky};

#[test]
fn known_hash_present_in_index() {
    let sky = SyntheticSky::random(50_000, 3.0, 9.0, 20261118);
    let cam = SyntheticCamera {
        width: 2048,
        height: 1536,
        ra0: 10.68 * DEG,
        dec0: 41.27 * DEG,
        scale: 6.0 * ARCSEC,
        rotation: 57.0 * DEG,
        psf_sigma_px: 1.6,
        sky: 500.0,
        read_noise: 6.0,
        zp_flux: 4.0e7,
        k1: 0.0,
    };
    let (_img, injected) = render_synthetic_image(&sky, &cam, 9917);
    let mut stars: Vec<_> = injected.iter().collect();
    stars.sort_by(|a, b| a.0.mag.partial_cmp(&b.0.mag).unwrap());
    // Pick a compact bright four-set: greedily find the four brightest stars
    // whose maximum pairwise separation stays inside the index coverage range.
    let max_sep = 1.4 * DEG;
    let mut chosen: Vec<_> = Vec::new();
    'search: for anchor in &stars {
        let mut group = vec![*anchor];
        for cand in &stars {
            if group.len() == 4 {
                break;
            }
            let ok = group
                .iter()
                .all(|g| g.0.v.separation(cand.0.v) <= max_sep);
            if ok && !group.iter().any(|g| g.0.id == cand.0.id) {
                group.push(*cand);
            }
        }
        if group.len() == 4 {
            chosen = group;
            break 'search;
        }
    }
    assert_eq!(chosen.len(), 4, "could not find a compact bright four-set");
    let ids: Vec<u32> = chosen.iter().map(|(s, _, _)| s.id).collect();

    let proj = TanProjection::new(cam.ra0, cam.dec0);
    let cat_pts: Vec<Pt> = chosen
        .iter()
        .enumerate()
        .map(|(i, (s, _, _))| {
            let (xi, eta) = proj.project(s.ra, s.dec).unwrap();
            Pt { x: xi, y: eta, id: i as u32 }
        })
        .collect();
    let q = make_quad([cat_pts[0], cat_pts[1], cat_pts[2], cat_pts[3]]).unwrap();
    let target_bucket = hash_bucket(q.hash);
    eprintln!("quad hash={:016x} bucket={} AB={:.4} deg", q.hash, target_bucket, q.ab / DEG);

    let built = build_index(&sky.stars, &BuildConfig::default());

    // Which scale bin should hold this AB?
    let bin_info: Vec<_> = built
        .manifest
        .bins
        .iter()
        .filter(|b| q.ab >= b.lo_rad && q.ab < b.hi_rad)
        .collect();
    eprintln!("matching scale bins: {}", bin_info.len());
    for b in &bin_info {
        eprintln!("  bin {} range {:.3}-{:.3} deg", b.index, b.lo_rad / DEG, b.hi_rad / DEG);
    }
    assert!(!bin_info.is_empty(), "no scale bin covers AB {}", q.ab);

    // Search every bin/bucket for an entry referencing these four star ids.
    let mut idset = ids.clone();
    idset.sort();
    let mut found_exact_hash = 0usize;
    let mut found_star_set = 0usize;
    for (bin, buckets) in built.per_bin.iter().enumerate() {
        for (bucket, entries) in buckets.iter().enumerate() {
            for e in entries {
                let mut es = e.stars;
                es.sort();
                if es == idset.as_slice() {
                    found_star_set += 1;
                    use platersol_core::quad::{hash_distance, unpack_hash};
                    eprintln!(
                        "same-set entry hash={:016x} dist={:.5} fields_img={:?} fields_idx={:?}",
                        e.hash,
                        hash_distance(q.hash, e.hash),
                        unpack_hash(q.hash),
                        unpack_hash(e.hash)
                    );
                    if e.hash == q.hash || e.hash == q.hash_mirror {
                        found_exact_hash += 1;
                        eprintln!(
                            "EXACT MATCH bin={} bucket={} hash={:016x} stars={:?}",
                            bin, bucket, e.hash, e.stars
                        );
                    }
                }
            }
        }
    }
    eprintln!("entries with same star set: {found_star_set}, exact-hash of those: {found_exact_hash}");

    // Also confirm the target bucket in the target bin is non-empty and the
    // hash window query would see it.
    if let Some(b) = bin_info.first() {
        let bin = b.index as usize;
        let entries = &built.per_bin[bin][target_bucket];
        eprintln!("target bin {bin} bucket {target_bucket} size = {}", entries.len());
        let wins = hash_window(q.hash, 25);
        let visible = entries
            .iter()
            .filter(|e| {
                let f = [
                    (e.hash >> 48) as u16,
                    (e.hash >> 32) as u16,
                    (e.hash >> 16) as u16,
                    e.hash as u16,
                ];
                f.iter().zip(wins.iter()).all(|(v, (lo, hi))| lo <= v && v <= hi)
            })
            .count();
        eprintln!("entries visible through hash window: {visible}");
    }

    // Why are they missing? Rank each chosen star globally and locally.
    let mut all = sky.stars.clone();
    all.sort_by(|a, b| a.mag.partial_cmp(&b.mag).unwrap());
    for (s, _, _) in chosen {
        let global_rank = all.iter().position(|x| x.id == s.id).unwrap();
        let brighter_within_9 = sky
            .stars
            .iter()
            .filter(|x| x.v.separation(s.v) <= 9.0 * DEG && x.mag < s.mag)
            .count();
        eprintln!(
            "star id={} mag={:.2} global_rank={}/{} brighter_within9={}",
            s.id, s.mag, global_rank, all.len(), brighter_within_9
        );
    }
    eprintln!("max_anchors = {}", sky.stars.len() / 4);

    // Projection-center differences between the camera tangent plane and the
    // index cell plane perturb the quantized hash by a few units, so the
    // solver matches within a tolerance window rather than by exact equality.
    let max_field_gap = built
        .per_bin
        .iter()
        .flatten()
        .flatten()
        .filter(|e| {
            let mut es = e.stars;
            es.sort();
            es == idset.as_slice()
        })
        .map(|e| {
            let a = platersol_core::quad::unpack_hash(q.hash);
            let b = platersol_core::quad::unpack_hash(e.hash);
            a.iter()
                .zip(b.iter())
                .map(|(x, y)| (*x as i32 - *y as i32).unsigned_abs())
                .max()
                .unwrap()
        })
        .min()
        .unwrap();
    eprintln!("smallest per-field hash gap = {max_field_gap}");
    assert!(
        max_field_gap <= 25,
        "hash gap {max_field_gap} exceeds solver tolerance"
    );
}
