//! On-disk / in-memory star catalog and a static 3-D KD-tree for fast
//! nearest-neighbour lookup on the unit sphere.

use crate::astro::Vec3;

/// A catalog star. `id` is stable across the catalog lifetime and is what quad
/// index entries reference.
#[derive(Clone, Copy, Debug)]
pub struct CatStar {
    pub id: u32,
    pub ra: f64,
    pub dec: f64,
    /// Apparent magnitude (smaller = brighter).
    pub mag: f32,
    /// Unit vector, cached for spatial queries.
    pub v: Vec3,
}

/// Serialization-friendly row (serde JSON is used for the small catalog; the
/// binary index stores everything else).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct CatStarRow {
    pub id: u32,
    pub ra_deg: f64,
    pub dec_deg: f64,
    pub mag: f32,
    /// Optional display name (HIP/Bayer/common), absent for anonymous stars.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
}

/// Static KD-tree over catalog stars, keyed by their unit 3-vectors.
/// Explicit node array: each node stores its star plus the indices of its
/// left/right children (`u32::MAX` = absent), split axis cycles x → y → z.
pub struct KdTree {
    nodes: Vec<KdNode>,
    stars: Vec<CatStar>,
}

#[derive(Clone, Copy)]
struct KdNode {
    star: CatStar,
    left: u32,
    right: u32,
    axis: u8,
}

impl KdTree {
    /// Build the tree by recursive median partitioning.
    pub fn build(stars: Vec<CatStar>) -> Self {
        struct Item {
            star: CatStar,
        }
        let mut items: Vec<Item> = stars.iter().map(|s| Item { star: *s }).collect();
        let mut nodes: Vec<KdNode> = Vec::with_capacity(items.len());

        fn axis_key(s: &CatStar, axis: u8) -> f64 {
            match axis {
                0 => s.v.x,
                1 => s.v.y,
                _ => s.v.z,
            }
        }

        fn build(items: &mut [Item], depth: u32, nodes: &mut Vec<KdNode>) -> u32 {
            let axis = (depth % 3) as u8;
            items.sort_by(|a, b| {
                axis_key(&a.star, axis)
                    .partial_cmp(&axis_key(&b.star, axis))
                    .unwrap()
            });
            let mid = items.len() / 2;
            // Reserve this node's slot up front so children get higher
            // indices and the returned id stays valid after pushes.
            let id = nodes.len() as u32;
            nodes.push(KdNode {
                star: items[mid].star,
                left: u32::MAX,
                right: u32::MAX,
                axis,
            });
            let (left_part, right_all) = items.split_at_mut(mid);
            let right_part = &mut right_all[1..];
            if !left_part.is_empty() {
                nodes[id as usize].left = build(left_part, depth + 1, nodes);
            }
            if !right_part.is_empty() {
                nodes[id as usize].right = build(right_part, depth + 1, nodes);
            }
            id
        }

        if items.is_empty() {
            return KdTree {
                nodes: Vec::new(),
                stars,
            };
        }
        build(&mut items, 0, &mut nodes);
        KdTree { nodes, stars }
    }

    pub fn len(&self) -> usize {
        self.stars.len()
    }

    pub fn is_empty(&self) -> bool {
        self.stars.is_empty()
    }

    pub fn stars(&self) -> &[CatStar] {
        &self.stars
    }

    /// All stars within angular `radius` (radians) of direction `center`.
    /// Chord distance `|a-b|² = 2(1-cos r)` replaces trig in the hot loop.
    pub fn within(&self, center: Vec3, radius: f64) -> Vec<CatStar> {
        let chord2 = 2.0 * (1.0 - radius.cos());
        let mut out = Vec::new();
        if self.nodes.is_empty() {
            return out;
        }
        fn walk(
            tree: &[KdNode],
            i: u32,
            c: Vec3,
            chord2: f64,
            out: &mut Vec<CatStar>,
        ) {
            let node = &tree[i as usize];
            let (dx, dy, dz) = (
                node.star.v.x - c.x,
                node.star.v.y - c.y,
                node.star.v.z - c.z,
            );
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 <= chord2 {
                out.push(node.star);
            }
            let diff = match node.axis {
                0 => dx,
                1 => dy,
                _ => dz,
            };
            let (near, far) = if diff > 0.0 {
                (node.left, node.right)
            } else {
                (node.right, node.left)
            };
            if near != u32::MAX {
                walk(tree, near, c, chord2, out);
            }
            if far != u32::MAX && diff * diff <= chord2 {
                walk(tree, far, c, chord2, out);
            }
        }
        walk(&self.nodes, 0, center, chord2, &mut out);
        out
    }

    /// Nearest star to `center`, if the tree is non-empty.
    pub fn nearest(&self, center: Vec3) -> Option<CatStar> {
        if self.nodes.is_empty() {
            return None;
        }
        fn d2_of(s: &CatStar, c: Vec3) -> f64 {
            let dx = s.v.x - c.x;
            let dy = s.v.y - c.y;
            let dz = s.v.z - c.z;
            dx * dx + dy * dy + dz * dz
        }
        fn walk(tree: &[KdNode], i: u32, c: Vec3) -> (f64, CatStar) {
            let node = &tree[i as usize];
            let mut best = (d2_of(&node.star, c), node.star);
            let diff = match node.axis {
                0 => node.star.v.x - c.x,
                1 => node.star.v.y - c.y,
                _ => node.star.v.z - c.z,
            };
            let (near, far) = if diff > 0.0 {
                (node.left, node.right)
            } else {
                (node.right, node.left)
            };
            if near != u32::MAX {
                let cand = walk(tree, near, c);
                if cand.0 < best.0 {
                    best = cand;
                }
            }
            if far != u32::MAX && diff * diff < best.0 {
                let cand = walk(tree, far, c);
                if cand.0 < best.0 {
                    best = cand;
                }
            }
            best
        }
        Some(walk(&self.nodes, 0, center).1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::{DEG, Vec3 as V};

    fn star(id: u32, ra_deg: f64, dec_deg: f64) -> CatStar {
        let (ra, dec) = (ra_deg * DEG, dec_deg * DEG);
        CatStar {
            id,
            ra,
            dec,
            mag: 6.0,
            v: V::from_radec(ra, dec),
        }
    }

    #[test]
    fn radius_query_sphere() {
        let stars: Vec<CatStar> = (0..200)
            .map(|i| {
                let ra = (i * 37 % 360) as f64;
                let dec = ((i * 53 % 180) as i64 - 90) as f64;
                star(i, ra, dec)
            })
            .collect();
        let tree = KdTree::build(stars);
        let center = V::from_radec(45.0 * DEG, 30.0 * DEG);
        let r = 15.0 * DEG;
        let got = tree.within(center, r);
        // Brute-force reference.
        let expect: Vec<u32> = tree
            .stars()
            .iter()
            .filter(|s| s.v.separation(center) <= r + 1e-12)
            .map(|s| s.id)
            .collect();
        let mut got_ids: Vec<u32> = got.iter().map(|s| s.id).collect();
        got_ids.sort();
        let mut expect_ids = expect;
        expect_ids.sort();
        assert_eq!(got_ids, expect_ids);
    }

    #[test]
    fn nearest_matches_brute_force() {
        let stars: Vec<CatStar> = (0..300)
            .map(|i| {
                star(
                    i,
                    (i * 23 % 360) as f64,
                    ((i * 41 % 160) as i64 - 80) as f64,
                )
            })
            .collect();
        let tree = KdTree::build(stars);
        let q = V::from_radec(123.0 * DEG, -40.0 * DEG);
        let got = tree.nearest(q).unwrap();
        let expect = tree
            .stars()
            .iter()
            .min_by(|a, b| {
                a.v.separation(q)
                    .partial_cmp(&b.v.separation(q))
                    .unwrap()
            })
            .unwrap();
        assert_eq!(got.id, expect.id);
    }
}
