//! Blind plate solver: image quads → geometric hash lookup → hypothesis
//! voting → WCS fit.
//!
//! The solver is storage-agnostic: an [`IndexStore`] implementation supplies
//! index buckets on demand. Native code uses a directory-backed store; the
//! WASM layer supplies chunks from browser cache / HTTP range requests.

use crate::astro::ARCSEC;
use crate::catalog::CatStar;
use crate::image::Source;
use crate::index::{IndexManifest, LoadedBin, ENTRY_SIZE};
use crate::quad::{make_quad, hash_window, Pt, Quad};
use crate::wcs::{fit_wcs, Correspondence, FitConfig, WcsSolution};
use std::collections::HashMap;

/// Abstract chunk storage. `load_bucket` returns the raw bytes of one bucket
/// run (the store may return a whole-bin slice and let the solver slice it).
pub trait IndexStore {
    fn manifest(&self) -> &IndexManifest;
    /// Raw bytes for bucket `bucket` of scale bin `bin`.
    fn load_bucket(&mut self, bin: usize, bucket: usize) -> std::io::Result<Vec<u8>>;
    /// Catalog star by id (needed to turn matches into correspondences).
    fn star_by_id(&self, id: u32) -> Option<CatStar>;
    /// All catalog ids held by this store (full catalog for native stores;
    /// the resident set for streaming WASM stores).
    fn known_ids(&self) -> Vec<u32>;
    /// Total chunk bytes fetched since open, for diagnostics.
    fn bytes_fetched(&self) -> u64;
}

/// Directory-backed store for native use.
pub struct DirStore {
    pub manifest: IndexManifest,
    pub dir: String,
    pub cache: HashMap<(usize, usize), Vec<u8>>,
    pub catalog: HashMap<u32, CatStar>,
    pub reads: u64,
    pub bytes_read: u64,
}

impl DirStore {
    pub fn open(dir: &str, catalog: &[CatStar]) -> std::io::Result<Self> {
        let manifest = crate::index::parse_manifest(&std::fs::read(format!("{dir}/index.json"))?)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        let catalog = catalog.iter().map(|s| (s.id, *s)).collect();
        Ok(DirStore {
            manifest,
            dir: dir.to_string(),
            cache: HashMap::new(),
            catalog,
            reads: 0,
            bytes_read: 0,
        })
    }
}

impl IndexStore for DirStore {
    fn manifest(&self) -> &IndexManifest {
        &self.manifest
    }

    fn load_bucket(&mut self, bin: usize, bucket: usize) -> std::io::Result<Vec<u8>> {
        if let Some(v) = self.cache.get(&(bin, bucket)) {
            return Ok(v.clone());
        }
        let offs = &self.manifest.bucket_offsets[bin];
        let (start, end) = (offs[bucket], offs[bucket + 1]);
        let path = format!("{}/bin_{:04}.bin", self.dir, bin);
        use std::io::{Read, Seek, SeekFrom};
        let mut f = std::fs::File::open(path)?;
        f.seek(SeekFrom::Start(start))?;
        let len = (end - start) as usize;
        let mut buf = vec![0u8; len];
        f.read_exact(&mut buf)?;
        self.reads += 1;
        self.bytes_read += len as u64;
        self.cache.insert((bin, bucket), buf.clone());
        Ok(buf)
    }

    fn star_by_id(&self, id: u32) -> Option<CatStar> {
        self.catalog.get(&id).copied()
    }

    fn known_ids(&self) -> Vec<u32> {
        self.catalog.keys().copied().collect()
    }

    fn bytes_fetched(&self) -> u64 {
        self.bytes_read
    }
}

/// Solver tuning parameters.
#[derive(Clone, Debug)]
pub struct SolveConfig {
    /// Brightest N image sources to build quads from.
    pub brightest: usize,
    /// Image A–B separations to probe, as fractions of the image diagonal.
    pub min_ab_frac: f64,
    pub max_ab_frac: f64,
    /// Hash window in 16-bit field units (both sides).
    pub hash_tol: u16,
    /// Max catalog candidates accepted per image quad.
    pub max_candidates_per_quad: usize,
    /// Stop once a mapping pair collects this many agreeing quads.
    pub quorum: usize,
    /// Max image quads attempted.
    pub max_image_quads: usize,
    /// Hard cap on bucket fetches per solve (the 30 s budget guard).
    pub max_bucket_reads: usize,
    /// WCS fit configuration.
    pub fit: FitConfig,
}

impl Default for SolveConfig {
    fn default() -> Self {
        SolveConfig {
            brightest: 40,
            // Compact quads: A-B between 3% and 55% of the diagonal. Small
            // baselines stay inside the index's local-cell coverage and are
            // less distorted by the tangent-plane nonlinearity.
            min_ab_frac: 0.03,
            max_ab_frac: 0.55,
            hash_tol: 25,
            max_candidates_per_quad: 8,
            quorum: 4,
            max_image_quads: 120,
            max_bucket_reads: 600,
            fit: FitConfig::default(),
        }
    }
}

/// Outcome of a successful solve plus diagnostics.
#[derive(Clone, Debug)]
pub struct SolveOutcome {
    pub wcs: WcsSolution,
    /// Matched pairs used for the final fit.
    pub matches: Vec<(usize, u32)>,
    pub quads_tried: usize,
    pub buckets_read: usize,
    pub index_bytes_read: u64,
}

/// Why a solve failed.
#[derive(Clone, Debug)]
pub enum SolveError {
    TooFewSources(usize),
    NoScaleBins,
    NoCandidates,
    NoConsistentHypothesis,
    FitFailed,
}

impl std::fmt::Display for SolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SolveError::TooFewSources(n) => {
                write!(f, "only {n} sources extracted (need >= 4)")
            }
            SolveError::NoScaleBins => write!(f, "no index scale bins overlap the image geometry"),
            SolveError::NoCandidates => write!(f, "no catalog candidates for any image quad"),
            SolveError::NoConsistentHypothesis => {
                write!(f, "candidates never reached a consistent hypothesis")
            }
            SolveError::FitFailed => write!(f, "WCS fit failed to converge"),
        }
    }
}

/// A single quad match hypothesis: image quad + catalog entry + the star
/// correspondence implied by the shared A,B,C,D ordering.
#[derive(Clone, Debug)]
struct Hypothesis {
    img_stars: [u32; 4],
    cat_stars: [u32; 4],
}

/// Run the full blind solve.
pub fn solve<S: IndexStore>(
    sources: &[Source],
    img_w: f64,
    img_h: f64,
    store: &mut S,
    cfg: &SolveConfig,
) -> Result<SolveOutcome, SolveError> {
    if sources.len() < 4 {
        return Err(SolveError::TooFewSources(sources.len()));
    }
    let n = cfg.brightest.min(sources.len());
    let pts: Vec<Pt> = sources[..n]
        .iter()
        .enumerate()
        .map(|(i, s)| Pt {
            x: s.x,
            y: s.y,
            id: i as u32,
        })
        .collect();

    let diag = (img_w * img_w + img_h * img_h).sqrt();
    let ab_lo = cfg.min_ab_frac * diag;
    let ab_hi = cfg.max_ab_frac * diag;

    // Image quads grouped by A-B separation (pixels).
    let img_quads = enumerate_image_quads(&pts, ab_lo, ab_hi, cfg.max_image_quads);
    if img_quads.is_empty() {
        return Err(SolveError::TooFewSources(n));
    }

    let manifest = store.manifest().clone();
    let n_bins = manifest.bins.len();
    let mut quads_tried = 0usize;
    let mut buckets_read = 0usize;
    let mut any_candidate = false;

    // Pair-level voting: (image source id, catalog star id) -> agreeing quads.
    let mut pair_votes: HashMap<(u32, u32), u32> = HashMap::new();
    let mut hypo_log: Vec<Hypothesis> = Vec::new();

    'outer: for q in &img_quads {
        // Both chiralities: direct hash keeps A,B,C,D order; the mirror hash
        // swaps the C and D roles to align with reflected index entries.
        let probes = [
            (q.hash, q.stars),
            (
                q.hash_mirror,
                [q.stars[0], q.stars[1], q.stars[3], q.stars[2]],
            ),
        ];

        // Blind scale search: walk every scale bin, but order them by how
        // close their implied plate scale is to a typical deep-sky image
        // (~2.5 arcsec/px). The ordering adapts per quad via its pixel A–B.
        let target = 2.5 * ARCSEC;
        let mut order: Vec<usize> = (0..n_bins).collect();
        order.sort_by_cached_key(|&b| {
            let bin = &manifest.bins[b];
            let mid = (bin.lo_rad + bin.hi_rad) * 0.5;
            let implied = mid / q.ab; // rad per pixel
            ((implied - target).abs() * 1e7) as i64
        });

        for &(probe_hash, probe_stars) in &probes {
            let windows = hash_window(probe_hash, cfg.hash_tol);
            // The tolerance window can straddle bucket boundaries (a bucket
            // spans 64 values of the first 16-bit field). Probe every bucket
            // the window's first field touches, not just the exact-hash one.
            let field0 = (probe_hash >> 48) as u16;
            let lo_bucket = field0.saturating_sub(cfg.hash_tol) >> 6;
            let hi_bucket = field0.saturating_add(cfg.hash_tol) >> 6;
            let bucket_ids: Vec<usize> = (lo_bucket..=hi_bucket)
                .map(|b| b as usize)
                .collect();
            for &b in &order {
                for &bucket in &bucket_ids {
                    let raw = match store.load_bucket(b, bucket) {
                        Ok(v) => v,
                        Err(_) => continue,
                    };
                    if raw.is_empty() {
                        continue;
                    }
                    buckets_read += 1;
                    let n_ent = raw.len() / ENTRY_SIZE;
                    let mut cands = 0usize;
                    for i in 0..n_ent {
                        let off = i * ENTRY_SIZE;
                        let hash = u64::from_le_bytes(raw[off..off + 8].try_into().unwrap());
                        if !hash_in_window(hash, windows) {
                            continue;
                        }
                        let mut cs = [0u32; 4];
                        for (k, slot) in cs.iter_mut().enumerate() {
                            let p = off + 8 + 4 * k;
                            *slot = u32::from_le_bytes(raw[p..p + 4].try_into().unwrap());
                        }
                        cands += 1;
                        any_candidate = true;
                        let hyp = Hypothesis {
                            img_stars: probe_stars,
                            cat_stars: cs,
                        };
                        for k in 0..4 {
                            *pair_votes
                                .entry((hyp.img_stars[k], hyp.cat_stars[k]))
                                .or_insert(0) += 1;
                        }
                        hypo_log.push(hyp);
                        if cands >= cfg.max_candidates_per_quad {
                            break;
                        }
                    }
                    quads_tried += 1;

                    // Early exit as soon as a quorum-strength mapping exists.
                    let strongest = pair_votes.values().copied().max().unwrap_or(0);
                    if strongest >= cfg.quorum as u32
                        && pair_votes.values().filter(|v| **v >= 2).count() >= 4
                    {
                        break 'outer;
                    }
                    if buckets_read >= cfg.max_bucket_reads {
                        break 'outer;
                    }
                }
            }
        }
    }
    let _ = hypo_log;

    if !any_candidate {
        return Err(SolveError::NoCandidates);
    }

    // Assemble the strongest consistent mapping: take pairs above quorum,
    // requiring a unique catalog star per image source and vice versa.
    let mut ranked: Vec<((u32, u32), u32)> = pair_votes
        .iter()
        .map(|(k, v)| (*k, *v))
        .filter(|(_, v)| *v >= 2)
        .collect();
    ranked.sort_by_key(|(_, v)| std::cmp::Reverse(*v));

    let mut img_to_cat: HashMap<u32, u32> = HashMap::new();
    let mut cat_used: HashMap<u32, u32> = HashMap::new();
    let mut votes: Vec<u32> = Vec::new();
    for ((im, ca), v) in ranked {
        if img_to_cat.contains_key(&im) || cat_used.contains_key(&ca) {
            continue;
        }
        img_to_cat.insert(im, ca);
        cat_used.insert(ca, im);
        votes.push(v);
    }

    let strong = votes.iter().filter(|v| **v >= cfg.quorum as u32).count();
    if img_to_cat.len() < 4 || strong < 2 {
        return Err(SolveError::NoConsistentHypothesis);
    }

    // Build correspondences and fit.
    let mut corrs: Vec<Correspondence> = Vec::new();
    let mut matches: Vec<(usize, u32)> = Vec::new();
    for (&im, &ca) in &img_to_cat {
        let Some(star) = store.star_by_id(ca) else {
            continue;
        };
        let s = &sources[im as usize];
        corrs.push(Correspondence {
            img_x: s.x,
            img_y: s.y,
            cat_ra: star.ra,
            cat_dec: star.dec,
        });
        matches.push((im as usize, ca));
    }
    if corrs.len() < 4 {
        return Err(SolveError::NoConsistentHypothesis);
    }
    let wcs = fit_wcs(&corrs, img_w, img_h, &cfg.fit).ok_or(SolveError::FitFailed)?;

    // Verify against *all* extracted sources via nearest catalog stars, not
    // just quad members: grow the match set under the current WCS and refit.
    let (wcs, matches, _n_inliers) =
        verify_and_extend(store, sources, img_w, img_h, &wcs, &cfg.fit);

    Ok(SolveOutcome {
        wcs,
        matches,
        quads_tried,
        buckets_read,
        index_bytes_read: store.bytes_fetched(),
    })
}

/// Enumerate image quads using the same local-neighborhood strategy the index
/// builder uses, so image and catalog quads share the same scale regime.
///
/// For each bright anchor, the nearby bright points form a spatial group
/// (radius = the largest acceptable A-B); its four-combinations yield
/// compact quads. This mirrors the overlapping-cell index construction so an
/// image quad built from real stars meets the corresponding catalog quad.
fn enumerate_image_quads(
    pts: &[Pt],
    ab_lo: f64,
    ab_hi: f64,
    cap: usize,
) -> Vec<Quad> {
    let mut out: Vec<Quad> = Vec::new();
    let m = pts.len();
    if m < 4 {
        return out;
    }
    let r2 = ab_hi * ab_hi;
    let group_cap = 12usize; // bound C(n,4) per neighborhood
    let mut seen_sets: std::collections::HashSet<[u32; 4]> = std::collections::HashSet::new();

    'outer: for i in 0..m {
        if out.len() >= cap {
            break;
        }
        let mut group: Vec<usize> = Vec::with_capacity(group_cap);
        for j in 0..m {
            let dx = pts[i].x - pts[j].x;
            let dy = pts[i].y - pts[j].y;
            if dx * dx + dy * dy <= r2 {
                group.push(j);
                if group.len() >= group_cap {
                    break;
                }
            }
        }
        if group.len() < 4 {
            continue;
        }
        for a in 0..group.len() {
            for b in (a + 1)..group.len() {
                for c in (b + 1)..group.len() {
                    for d in (c + 1)..group.len() {
                        let Some(q) = make_quad([
                            pts[group[a]],
                            pts[group[b]],
                            pts[group[c]],
                            pts[group[d]],
                        ]) else {
                            continue;
                        };
                        if q.ab < ab_lo || q.ab > ab_hi {
                            continue;
                        }
                        let mut key = q.stars;
                        key.sort();
                        if seen_sets.insert(key) {
                            out.push(q);
                            if out.len() >= cap {
                                break 'outer;
                            }
                        }
                    }
                }
            }
        }
    }
    out
}

#[inline]
fn hash_in_window(hash: u64, w: [(u16, u16); 4]) -> bool {
    let f = [
        (hash >> 48) as u16,
        (hash >> 32) as u16,
        (hash >> 16) as u16,
        hash as u16,
    ];
    f.iter().zip(w.iter()).all(|(v, (lo, hi))| *lo <= *v && *v <= *hi)
}

/// Verify a WCS against every extracted source: project each source, find the
/// catalog star nearest to that sky position, accept within a tightening
/// tolerance, then refit on the enlarged inlier set.
fn verify_and_extend<S: IndexStore>(
    store: &S,
    sources: &[Source],
    img_w: f64,
    img_h: f64,
    wcs: &WcsSolution,
    fit_cfg: &FitConfig,
) -> (WcsSolution, Vec<(usize, u32)>, usize) {
    // Build a local KD-tree from catalog stars reachable through the store.
    // Stores expose point lookups; the DirStore holds the full catalog map.
    // For generic stores we re-project via star_by_id over ids seen in quads,
    // which the WASM store also keeps resident.
    let all_ids: Vec<u32> = store.known_ids();
    let stars: Vec<CatStar> = all_ids
        .iter()
        .filter_map(|id| store.star_by_id(*id))
        .collect();
    let tree = crate::catalog::KdTree::build(stars);

    let tol = (wcs.rms_arcsec * 3.0).clamp(2.0, 30.0) * ARCSEC;
    let mut corrs = Vec::new();
    let mut matches = Vec::new();
    for (i, s) in sources.iter().enumerate() {
        let (ra, dec) = wcs.pixel_to_radec(s.x, s.y);
        let v = crate::astro::Vec3::from_radec(ra, dec);
        if let Some(near) = tree.nearest(v) {
            if near.v.separation(v) <= tol {
                corrs.push(Correspondence {
                    img_x: s.x,
                    img_y: s.y,
                    cat_ra: near.ra,
                    cat_dec: near.dec,
                });
                matches.push((i, near.id));
            }
        }
    }
    if corrs.len() >= 6 {
        if let Some(fit2) = fit_wcs(&corrs, img_w, img_h, fit_cfg) {
            return (fit2, matches, corrs.len());
        }
    }
    (wcs.clone(), matches, corrs.len())
}

/// Helper used by the WASM store and tests: parse a whole-bin byte blob
/// (header included) into entries.
pub fn parse_full_bin(bytes: &[u8]) -> std::io::Result<LoadedBin> {
    LoadedBin::parse(bytes)
}
