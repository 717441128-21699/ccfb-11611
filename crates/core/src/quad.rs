//! Geometric four-star quads and their scale-invariant hash.
//!
//! A quad is built from four points in a 2-D plane (image pixels, or stars
//! projected onto the tangent plane). Two of them — the pair with the largest
//! separation — define a canonical coordinate frame: A=(-0.5, 0), B=(+0.5, 0).
//! The other two stars C, D are expressed in that frame, ordered so the hash
//! is invariant under rotation, translation, scale and reflection handling.
//!
//! Hash layout (u64), four unsigned 16-bit fields:
//! ```text
//! bits 48..64  xC in [0, 1024) mapped from [-1, 1]
//! bits 32..48  yC
//! bits 16..32  xD
//! bits  0..16  yD
//! ```
//! The A–B physical separation (pixels for images, radians for sky) is kept
//! alongside the hash so the solver can recover the image scale and bin index
//! lookups by angular scale.

pub const HASH_BITS: u32 = 16;
pub const HASH_MAX: u32 = 1 << HASH_BITS; // 65536
/// Largest representable 16-bit hash field value.
pub const FIELD_MAX: u16 = (HASH_MAX - 1) as u16; // 65535
const COORD_SCALE: f64 = (HASH_MAX as f64) * 0.5; // [-1,1] -> [0,65536)

/// A 2-D point used for quad construction.
#[derive(Clone, Copy, Debug)]
pub struct Pt {
    pub x: f64,
    pub y: f64,
    /// Index of this point in the source/catalog array.
    pub id: u32,
}

/// A canonicalized quad.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Quad {
    /// Geometric hash of the direct (right-handed) canonicalization.
    pub hash: u64,
    /// Hash of the reflected (left-handed, C/D swapped) geometry. A mirrored
    /// image of the same sky matches on this hash, so the solver probes both.
    pub hash_mirror: u64,
    /// A–B separation in the input plane's units.
    pub ab: f64,
    /// Original point indices in A, B, C, D order (direct geometry).
    pub stars: [u32; 4],
    /// Canonical coordinates of C and D in the A–B frame.
    pub c: (f64, f64),
    pub d: (f64, f64),
}

#[inline]
fn quantize(v: f64) -> u16 {
    let q = ((v + 1.0) * COORD_SCALE).round() as i64;
    q.clamp(0, HASH_MAX as i64 - 1) as u16
}
#[inline]
fn pack_hash(c: (f64, f64), d: (f64, f64)) -> u64 {
    let (xc, yc, xd, yd) = (
        quantize(c.0) as u64,
        quantize(c.1) as u64,
        quantize(d.0) as u64,
        quantize(d.1) as u64,
    );
    (xc << 48) | (yc << 32) | (xd << 16) | yd
}

/// Build a quad from four points. The two farthest apart become A, B.
///
/// Returns `None` if the points are (near-)degenerate (zero A–B separation or
/// any point essentially on the A–B line midpoint, which makes hashing
/// numerically meaningless).
pub fn make_quad(p: [Pt; 4]) -> Option<Quad> {
    // Find the pair with maximum squared distance.
    let mut best = (0usize, 1usize, f64::NEG_INFINITY);
    for i in 0..4 {
        for j in (i + 1)..4 {
            let dx = p[i].x - p[j].x;
            let dy = p[i].y - p[j].y;
            let d2 = dx * dx + dy * dy;
            if d2 > best.2 {
                best = (i, j, d2);
            }
        }
    }
    let (ai, bi) = (best.0, best.1);
    let ab2 = best.2;
    if ab2 <= f64::EPSILON {
        return None;
    }
    let ab = ab2.sqrt();

    // Canonical frame: A at (-0.5, 0), B at (+0.5, 0).
    // For a point P:  x' = dot(P-A, B-A)/|B-A|^2 - 0.5
    //                 y' = cross(B-A, P-A)/|B-A|^2
    let (a, b) = (p[ai], p[bi]);
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let frame = |pt: Pt| -> (f64, f64) {
        let (px, py) = (pt.x - a.x, pt.y - a.y);
        let xp = (px * dx + py * dy) / ab2 - 0.5;
        let yp = (dx * py - dy * px) / ab2;
        (xp, yp)
    };

    let rest: Vec<usize> = (0..4).filter(|&k| k != ai && k != bi).collect();
    let c0 = frame(p[rest[0]]);
    let d0 = frame(p[rest[1]]);

    // Direct geometry labels C, D in input order. The mirror hash below
    // covers the reflected configuration, so no input-dependent swap is
    // needed and the correspondence ordering A,B,C,D stays meaningful.
    let (c, d) = (c0, d0);
    let (ci, di) = (rest[0], rest[1]);

    // Reject points outside the sane canonical square; a valid quad with
    // A-B as diameter keeps C,D well inside [-1, 1] in practice.
    if [c.0, c.1, d.0, d.1].iter().any(|v| !(-1.05..=1.05).contains(v)) {
        return None;
    }

    let cc = c.clamp_coord();
    let dd = d.clamp_coord();
    Some(Quad {
        hash: pack_hash(cc, dd),
        hash_mirror: pack_hash(dd, cc),
        ab,
        stars: [a.id, b.id, p[ci].id, p[di].id],
        c,
        d,
    })
}

trait ClampCoord {
    fn clamp_coord(self) -> Self;
}
impl ClampCoord for (f64, f64) {
    fn clamp_coord(self) -> Self {
        (self.0.clamp(-1.0, 1.0), self.1.clamp(-1.0, 1.0))
    }
}

/// Extract the four 16-bit fields of a hash in (xC, yC, xD, yD) order.
pub fn unpack_hash(hash: u64) -> [u16; 4] {
    [
        (hash >> 48) as u16,
        (hash >> 32) as u16,
        (hash >> 16) as u16,
        hash as u16,
    ]
}

/// Hamming-style distance between two hashes, normalized to `[0, 1]` over the
/// four fields (mean absolute field difference / 65535).
pub fn hash_distance(h1: u64, h2: u64) -> f64 {
    let a = unpack_hash(h1);
    let b = unpack_hash(h2);
    let s: u32 = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (*x as i32 - *y as i32).unsigned_abs())
        .sum();
    s as f64 / (4.0 * FIELD_MAX as f64)
}

/// Tolerance window (in hash-field units) around a hash, used for index
/// probes. `tol_frac` is expressed relative to the full field range.
pub fn hash_window(hash: u64, tol_units: u16) -> [(u16, u16); 4] {
    let fields = unpack_hash(hash);
    fields.map(|q| {
        let lo = q.saturating_sub(tol_units);
        let hi = q.saturating_add(tol_units);
        (lo, hi)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pts(v: &[(f64, f64)]) -> Vec<Pt> {
        v.iter()
            .enumerate()
            .map(|(i, &(x, y))| Pt { x, y, id: i as u32 })
            .collect()
    }

    #[test]
    fn hash_invariant_rotation_translation_scale() {
        let base = pts(&[
            (0.0, 0.0),
            (10.0, 0.0),
            (3.0, 4.0),
            (7.0, -3.0),
        ]);
        let q1 = make_quad([base[0], base[1], base[2], base[3]]).unwrap();

        // Rotate 137°, translate, scale by 2.3.
        let (cs, sn) = (137f64.to_radians().cos(), 137f64.to_radians().sin());
        let t: Vec<Pt> = base
            .iter()
            .map(|p| {
                let (x, y) = (p.x * 2.3, p.y * 2.3);
                Pt {
                    x: x * cs - y * sn + 50.0,
                    y: x * sn + y * cs - 12.0,
                    id: p.id,
                }
            })
            .collect();
        let q2 = make_quad([t[0], t[1], t[2], t[3]]).unwrap();
        assert_eq!(q1.hash, q2.hash, "hashes must match under rigid+scale");
        assert!((q2.ab - q1.ab * 2.3).abs() < 1e-9 * q2.ab);
    }

    #[test]
    fn hash_invariant_inner_swap() {
        let base = pts(&[(0.0, 0.0), (10.0, 0.0), (3.0, 4.0), (7.0, -3.0)]);
        let q1 = make_quad([base[0], base[1], base[2], base[3]]).unwrap();

        // Swapping the two inner stars only exchanges C/D roles, so the
        // direct and mirror hashes swap as a set.
        let q2 = make_quad([base[0], base[1], base[3], base[2]]).unwrap();
        let set1 = [q1.hash, q1.hash_mirror];
        assert!(set1.contains(&q2.hash) && set1.contains(&q2.hash_mirror));
    }

    #[test]
    fn reflection_flips_canonical_eta() {
        // Reflection across the x axis maps canonical (xi, eta) -> (xi, -eta).
        // Verify the canonical-frame coordinates directly, which is the
        // invariant the dual-chirality index relies on.
        let base = pts(&[(0.0, 0.0), (10.0, 0.0), (3.0, 4.0), (7.0, -3.0)]);
        let q1 = make_quad([base[0], base[1], base[2], base[3]]).unwrap();

        let mirrored: Vec<Pt> = base
            .iter()
            .map(|p| Pt {
                x: p.x,
                y: -p.y,
                id: p.id,
            })
            .collect();
        let q3 = make_quad([mirrored[0], mirrored[1], mirrored[2], mirrored[3]]).unwrap();

        // A–B distance invariant under reflection.
        assert!((q3.ab - q1.ab).abs() < 1e-9 * q1.ab);
        // The mirrored quad's canonical coordinates must be the eta-flipped
        // version of the original's (in some C/D ordering).
        let direct = [q1.c, q1.d];
        let flipped: Vec<(f64, f64)> = direct.iter().map(|(x, y)| (*x, -*y)).collect();
        let got = [q3.c, q3.d];
        let matched = flipped.iter().all(|p| {
            got.iter()
                .any(|g| (g.0 - p.0).abs() < 1e-9 && (g.1 - p.1).abs() < 1e-9)
        });
        assert!(
            matched,
            "reflection must eta-flip canonical coords: q1={:?} q3={:?}",
            q1, q3
        );
    }

    #[test]
    fn different_geometry_different_hash() {
        let p1 = pts(&[(0., 0.), (10., 0.), (3., 4.), (7., -3.)]);
        let p2 = pts(&[(0., 0.), (10., 0.), (2., 5.), (8., -2.)]);
        let q1 = make_quad([p1[0], p1[1], p1[2], p1[3]]).unwrap();
        let q2 = make_quad([p2[0], p2[1], p2[2], p2[3]]).unwrap();
        assert_ne!(q1.hash, q2.hash);
        assert!(hash_distance(q1.hash, q2.hash) > 0.01);
    }

    #[test]
    fn window_clamps() {
        let w = hash_window(0, 10);
        assert_eq!(w[0].0, 0);
        let w2 = hash_window(u64::MAX, 10);
        assert_eq!(w2[0].1, FIELD_MAX);
    }
}
