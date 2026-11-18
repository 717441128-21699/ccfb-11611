//! Geometric-hash quad index: binary chunk format, in-memory lookup and the
//! catalog-side builder.
//!
//! ## Scale bins
//! Quads are grouped by the angular A–B separation into logarithmic bins
//! (factor ~1.5 per bin). A blind solver does not know the field of view, so
//! it probes bins across a configurable range instead of touching the whole
//! index.
//!
//! ## Chunking (on-disk layout, friendly to HTTP Range requests)
//! Within one scale bin, entries are sorted by hash and split into buckets by
//! the top `BUCKET_BITS` bits of the hash. Every bucket is one contiguous run
//! in a single `.bin` file; the manifest records each bucket's byte offset and
//! length, so the UI can fetch exactly one bucket of one scale bin with one
//! range request.
//!
//! ## Entry (24 bytes, little-endian)
//! ```text
//! u64 hash
//! u32 star_id A, B, C, D
//! ```

use crate::astro::{TanProjection, Vec3, ARCSEC, DEG};
use crate::catalog::{CatStar, KdTree};
use crate::quad::{make_quad, Pt};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::{Read, Write};

pub const MAGIC: u32 = 0x50_53_51_31; // "PSQ1"
pub const BUCKET_BITS: u32 = 10;
pub const N_BUCKETS: usize = 1 << BUCKET_BITS; // 1024
pub const ENTRY_SIZE: usize = 24;
pub const HEADER_SIZE: usize = 36;

#[inline]
pub fn hash_bucket(hash: u64) -> usize {
    (hash >> (64 - BUCKET_BITS)) as usize
}

/// One scale bin: A–B separation in radians lies in `[lo, hi)`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ScaleBin {
    pub index: u32,
    pub lo_rad: f64,
    pub hi_rad: f64,
}

/// Manifest describing every chunk file of a built index.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct IndexManifest {
    pub magic: u32,
    pub version: u32,
    /// Catalog identifier / provenance.
    pub catalog: String,
    pub n_stars: u32,
    pub n_quads: u64,
    pub bins: Vec<ScaleBin>,
    /// `bucket_offsets[bin][bucket]` = byte offset inside that bin's `.bin`.
    /// Length is `bins.len() * N_BUCKETS`; a zero-length bucket has
    /// `off[i] == off[i+1]`.
    pub bucket_offsets: Vec<Vec<u64>>,
    /// Hash tolerance baked into the build (fraction of one 16-bit field),
    /// purely advisory — the solver picks its own tolerance.
    pub build_tol_frac: f64,
}

/// An index entry as stored on disk.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IndexEntry {
    pub hash: u64,
    pub stars: [u32; 4],
}

impl IndexEntry {
    pub fn write_to(&self, w: &mut impl Write) -> std::io::Result<()> {
        w.write_all(&self.hash.to_le_bytes())?;
        for s in self.stars {
            w.write_all(&s.to_le_bytes())?;
        }
        Ok(())
    }

    pub fn read_from(r: &mut impl Read) -> std::io::Result<Self> {
        let mut hb = [0u8; 8];
        r.read_exact(&mut hb)?;
        let hash = u64::from_le_bytes(hb);
        let mut stars = [0u32; 4];
        for s in stars.iter_mut() {
            let mut b = [0u8; 4];
            r.read_exact(&mut b)?;
            *s = u32::from_le_bytes(b);
        }
        Ok(IndexEntry { hash, stars })
    }
}

/// In-memory view of one scale bin's entries (typically only the buckets
/// touched by a solve are ever loaded).
#[derive(Clone, Debug, Default)]
pub struct LoadedBin {
    pub entries: Vec<IndexEntry>,
}

impl LoadedBin {
    /// Parse a raw chunk byte slice (header + entries as laid out on disk).
    pub fn parse(bytes: &[u8]) -> std::io::Result<Self> {
        if bytes.len() < HEADER_SIZE {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "chunk too short",
            ));
        }
        let magic = u32::from_le_bytes(bytes[0..4].try_into().unwrap());
        if magic != MAGIC {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "bad chunk magic",
            ));
        }
        let n = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
        let mut cur = &bytes[HEADER_SIZE..];
        let mut entries = Vec::with_capacity(n);
        for _ in 0..n {
            entries.push(IndexEntry::read_from(&mut cur)?);
        }
        Ok(LoadedBin { entries })
    }

    /// Range query: all entries whose four hash fields fall inside the given
    /// per-field windows. Entries are hash-sorted so the scan is linear over
    /// a typically tiny bucket.
    pub fn query(
        &self,
        windows: [(u16, u16); 4],
        max_results: usize,
    ) -> Vec<IndexEntry> {
        let mut out = Vec::new();
        for e in &self.entries {
            let h = e.hash;
            let f = [
                (h >> 48) as u16,
                (h >> 32) as u16,
                (h >> 16) as u16,
                h as u16,
            ];
            let ok = f
                .iter()
                .zip(windows.iter())
                .all(|(v, (lo, hi))| *lo <= *v && *v <= *hi);
            if ok {
                out.push(*e);
                if out.len() >= max_results {
                    break;
                }
            }
        }
        out
    }
}

/// Builder configuration.
#[derive(Clone, Debug)]
pub struct BuildConfig {
    /// Log-scale bins span this range in radians. Defaults cover ~10
    /// arcminutes to ~20 degrees of A–B separation.
    pub min_ab_rad: f64,
    pub max_ab_rad: f64,
    pub bins_per_decade: u32,
    /// Brightest N stars to consider inside each local sky cell.
    pub stars_per_cell: usize,
    /// Cell size (radians) for local tangent-plane enumeration.
    pub cell_size_rad: f64,
    /// Skip a quad when its smallest canonical inner-coordinate magnitude is
    /// below this — rejects skinny, numerically fragile quads.
    pub min_inner: f64,
}

impl Default for BuildConfig {
    fn default() -> Self {
        BuildConfig {
            // Quad A-B separations the overlapping cells can cover.
            min_ab_rad: 10.0 * 60.0 * ARCSEC,
            max_ab_rad: 1.5 * DEG,
            bins_per_decade: 5,
            // Brightest neighbors enumerated per small cell. C(10,4)=210
            // quads per cell before scale/dedup filtering.
            stars_per_cell: 10,
            cell_size_rad: 1.0 * DEG,
            min_inner: 0.01,
        }
    }
}

/// Built but not yet serialized: entries grouped by scale bin and bucket.
pub struct BuiltIndex {
    pub manifest: IndexManifest,
    /// `per_bin[bin][bucket]` entries.
    pub per_bin: Vec<Vec<Vec<IndexEntry>>>,
}

/// Build the full quad index from a catalog.
///
/// The sphere is tiled into RA/Dec cells; within each cell the brightest
/// `stars_per_cell` stars are projected onto a tangent plane centered on the
/// cell, and every 4-combination is hashed. Quads whose A–B separation falls
/// in a configured scale bin are stored. Duplicate quads (same four stars,
/// produced by overlapping cells) are deduplicated.
pub fn build_index(catalog: &[CatStar], cfg: &BuildConfig) -> BuiltIndex {
    // Scale bins.
    let log_lo = cfg.min_ab_rad.log10();
    let log_hi = cfg.max_ab_rad.log10();
    let n_bins = ((log_hi - log_lo) * cfg.bins_per_decade as f64).ceil() as usize;
    let mut bins = Vec::with_capacity(n_bins);
    for i in 0..n_bins {
        let lo = 10f64.powf(log_lo + i as f64 / cfg.bins_per_decade as f64);
        let hi = 10f64.powf(log_lo + (i + 1) as f64 / cfg.bins_per_decade as f64);
        bins.push(ScaleBin {
            index: i as u32,
            lo_rad: lo,
            hi_rad: hi,
        });
    }

    let tree = KdTree::build(catalog.to_vec());
    // Overlapping small-cell grid. Cell spacing is roughly the largest quad
    // A-B we want to cover; the KD query radius extends ~0.7 cell beyond each
    // center so a four-set spanning a cell border is still enumerated by a
    // neighboring cell. Small cells keep per-cell star counts bounded
    // without any magnitude truncation, preserving coverage of faint stars.
    let cell = cfg.cell_size_rad;
    // Query radius just past the cell half-width so a four-set straddling a
    // cell border is still captured, while minimizing the overlap that
    // regenerates the same quad in adjacent cells.
    let query_r = cell * 1.05;
    let n_dec = (std::f64::consts::PI / cell).ceil() as usize;
    let n_ra_equator = (2.0 * std::f64::consts::PI / cell).ceil() as usize;

    // (bin, 4 star ids) dedup set.
    let mut seen: std::collections::HashSet<(u32, [u32; 4])> = std::collections::HashSet::new();
    let mut per_bin: Vec<Vec<IndexEntry>> = vec![Vec::new(); n_bins];

    for idec in 0..n_dec {
        let dec_c = -std::f64::consts::FRAC_PI_2 + (idec as f64 + 0.5) * cell;
        // Narrow RA spacing toward the poles to keep physical cell area sane.
        let cos_dec = dec_c.cos().max(0.05);
        let n_ra = (n_ra_equator as f64 * cos_dec).ceil().max(1.0) as usize;
        for ira in 0..n_ra {
            let ra_c = (ira as f64 + 0.5) * (2.0 * std::f64::consts::PI / n_ra as f64);
            let center = Vec3::from_radec(ra_c, dec_c);
            let mut nearby = tree.within(center, query_r);
            // Bound per-cell combinatorics by keeping the locally brightest
            // stars; the small cell size preserves co-visibility coverage.
            nearby.sort_by(|a, b| a.mag.partial_cmp(&b.mag).unwrap());
            nearby.truncate(cfg.stars_per_cell);
            if nearby.len() < 4 {
                continue;
            }
            let proj = TanProjection::new(ra_c, dec_c);
            let pts: Vec<Pt> = nearby
                .iter()
                .filter_map(|s| {
                    proj.project(s.ra, s.dec).map(|(xi, eta)| Pt {
                        x: xi,
                        y: eta,
                        id: s.id,
                    })
                })
                .collect();
            if pts.len() < 4 {
                continue;
            }

            for i in 0..pts.len() {
                for j in (i + 1)..pts.len() {
                    for k in (j + 1)..pts.len() {
                        for l in (k + 1)..pts.len() {
                            let Some(q) = make_quad([pts[i], pts[j], pts[k], pts[l]]) else {
                                continue;
                            };
                            if q.c.0.abs() < cfg.min_inner && q.c.1.abs() < cfg.min_inner {
                                continue;
                            }
                            let Some(bin) = bins
                                .iter()
                                .position(|b| q.ab >= b.lo_rad && q.ab < b.hi_rad)
                            else {
                                continue;
                            };
                            // Both chiralities: direct A,B,C,D and mirrored
                            // A,B,D,C so reflected images map with correct
                            // per-role star identities.
                            let mut key = q.stars;
                            key.sort();
                            if seen.insert((bin as u32, key)) {
                                per_bin[bin].push(IndexEntry {
                                    hash: q.hash,
                                    stars: q.stars,
                                });
                                per_bin[bin].push(IndexEntry {
                                    hash: q.hash_mirror,
                                    stars: [q.stars[0], q.stars[1], q.stars[3], q.stars[2]],
                                });
                            }
                        }
                    }
                }
            }
        }
    }

    let mut total = 0u64;
    let mut bucketed = Vec::with_capacity(n_bins);
    let mut offsets = Vec::with_capacity(n_bins);
    for bin_entries in per_bin.iter_mut() {
        bin_entries.sort_by_key(|a| a.hash);
        bin_entries.dedup_by(|a, b| a.hash == b.hash && a.stars == b.stars);
        total += bin_entries.len() as u64;

        let mut buckets: Vec<Vec<IndexEntry>> = vec![Vec::new(); N_BUCKETS];
        for e in bin_entries.drain(..) {
            buckets[hash_bucket(e.hash)].push(e);
        }
        let mut off = Vec::with_capacity(N_BUCKETS + 1);
        let mut cursor = HEADER_SIZE as u64;
        for b in &buckets {
            off.push(cursor);
            cursor += (b.len() * ENTRY_SIZE) as u64;
        }
        off.push(cursor);
        offsets.push(off);
        bucketed.push(buckets);
    }

    let manifest = IndexManifest {
        magic: MAGIC,
        version: 1,
        catalog: "default".to_string(),
        n_stars: catalog.len() as u32,
        n_quads: total,
        bins,
        bucket_offsets: offsets,
        build_tol_frac: 0.004,
    };
    BuiltIndex {
        manifest,
        per_bin: bucketed,
    }
}

impl BuiltIndex {
    /// Serialize one scale bin to its binary chunk:
    /// 36-byte header (`magic u32, bin u32, reserved u32, n_entries u32,
    /// total_bytes u64, reserved 12B`) followed by sorted entries.
    pub fn write_bin(&self, bin: usize) -> Vec<u8> {
        let buckets = &self.per_bin[bin];
        let n: usize = buckets.iter().map(|b| b.len()).sum();
        let mut out = Vec::with_capacity(HEADER_SIZE + n * ENTRY_SIZE);
        out.extend_from_slice(&MAGIC.to_le_bytes());
        out.extend_from_slice(&(bin as u32).to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes());
        out.extend_from_slice(&(n as u32).to_le_bytes());
        out.extend_from_slice(&((HEADER_SIZE + n * ENTRY_SIZE) as u64).to_le_bytes());
        out.extend_from_slice(&[0u8; 12]);
        debug_assert_eq!(out.len(), HEADER_SIZE);
        for bucket in buckets {
            for e in bucket {
                e.write_to(&mut out).unwrap();
            }
        }
        out
    }

    /// Write manifest (`index.json`) and one `bin_XXXX.bin` per scale bin
    /// into `dir`.
    pub fn write_to_dir(&self, dir: &str) -> std::io::Result<()> {
        std::fs::create_dir_all(dir)?;
        let manifest_path = format!("{dir}/index.json");
        std::fs::write(
            &manifest_path,
            serde_json::to_vec_pretty(&self.manifest).unwrap(),
        )?;
        for bin in 0..self.manifest.bins.len() {
            let path = format!("{dir}/bin_{:04}.bin", bin);
            std::fs::write(path, self.write_bin(bin))?;
        }
        Ok(())
    }
}

/// Read a manifest from parsed JSON bytes.
pub fn parse_manifest(bytes: &[u8]) -> serde_json::Result<IndexManifest> {
    serde_json::from_slice(bytes)
}

/// Choose scale bins whose `[lo, hi)` overlaps `[ab_lo, ab_hi]` (radians).
pub fn bins_for_range(
    manifest: &IndexManifest,
    ab_lo: f64,
    ab_hi: f64,
) -> Vec<&ScaleBin> {
    manifest
        .bins
        .iter()
        .filter(|b| b.hi_rad > ab_lo && b.lo_rad < ab_hi)
        .collect()
}

/// Aggregate counts per bucket for diagnostics / load planning.
pub fn bucket_counts(built: &BuiltIndex) -> HashMap<usize, usize> {
    let mut m = HashMap::new();
    for bins in &built.per_bin {
        for (i, b) in bins.iter().enumerate() {
            if !b.is_empty() {
                *m.entry(i).or_insert(0) += b.len();
            }
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::DEG;

    fn uniform_catalog(n: u32) -> Vec<CatStar> {
        (0..n)
            .map(|i| {
                let ra = (i as f64 * 2.39996) % (2.0 * std::f64::consts::PI);
                let dec = ((i as f64 * 1.71347) % std::f64::consts::PI)
                    - std::f64::consts::FRAC_PI_2;
                let v = Vec3::from_radec(ra, dec);
                CatStar {
                    id: i,
                    ra,
                    dec,
                    mag: 5.0 + (i % 20) as f32 * 0.2,
                    v,
                }
            })
            .collect()
    }

    #[test]
    fn entry_roundtrip() {
        let e = IndexEntry {
            hash: 0x1234_5678_9abc_def0,
            stars: [1, 2, 3, 4],
        };
        let mut buf = Vec::new();
        e.write_to(&mut buf).unwrap();
        let got = IndexEntry::read_from(&mut &buf[..]).unwrap();
        assert_eq!(e, got);
    }

    #[test]
    fn build_and_query_roundtrip() {
        let cats = uniform_catalog(3000);
        let cfg = BuildConfig::default();
        let built = build_index(&cats, &cfg);
        assert!(built.manifest.n_quads > 100);
        // Serialize + parse every non-empty bin and verify a window query
        // returns at least the exact-hash entries.
        for bin in 0..built.manifest.bins.len() {
            let raw = built.write_bin(bin);
            let loaded = LoadedBin::parse(&raw).unwrap();
            // Two chirality entries are stored per unique quad.
            let expect: usize = built.per_bin[bin].iter().map(|b| b.len()).sum();
            assert_eq!(loaded.entries.len(), expect);
            if let Some(first) = loaded.entries.first() {
                let f = [
                    (first.hash >> 48) as u16,
                    (first.hash >> 32) as u16,
                    (first.hash >> 16) as u16,
                    first.hash as u16,
                ];
                let wins = f.map(|q| (q, q));
                let got = loaded.query(wins, 1000);
                assert!(got.iter().any(|e| e.hash == first.hash));
            }
        }
    }

    #[test]
    fn scale_bins_span_configured_range() {
        let cfg = BuildConfig::default();
        let cats = uniform_catalog(200);
        let built = build_index(&cats, &cfg);
        let first = &built.manifest.bins[0];
        let last = built.manifest.bins.last().unwrap();
        assert!((first.lo_rad - cfg.min_ab_rad).abs() < 1e-12);
        assert!(last.hi_rad >= cfg.max_ab_rad * 0.99);
        let _ = DEG;
    }
}
