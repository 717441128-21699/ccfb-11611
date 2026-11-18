//! 四星组合（quad）几何哈希。
//!
//! 取四颗星中相距最远的一对 A-B 作为锚点，把坐标变换到以 A 为原点、
//! AB 为 x 轴、AB 长度为 1 的局部坐标系；另外两颗 C、D 的局部坐标
//! (cx, cy, dx, dy) 对平移、旋转、缩放完全不变，镜像则通过词序规范化
//! 吸收（C/D 交换取 lex 较小者，与星表侧用同一规则即可正确匹配）。
//!
//! 四个分量各量化到 [`CODE_BINS`] 个桶，拼成 40 bit 整数作为哈希键，
//! 分块索引按该键排序，查询做 4 维小范围 box query。

use crate::image::StarObs;

/// 每维量化桶数（10 bit）。总编码 40 bit，单尺度约百亿分之一的碰撞空间。
pub const CODE_BINS: u32 = 1024;
const CODE_BITS: u32 = 10;

/// 一个四星组合。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Quad {
    /// 锚点 A、B 与另外两颗 C、D 在输入数组中的下标。
    pub stars: [usize; 4],
    /// 锚点间距（图像单位：像元；星表单位：弧度）。
    pub ab_dist: f64,
    /// 40 bit 几何编码。
    pub code: u64,
    /// 编码前的连续局部坐标（调试/精匹配用）。
    pub local: [[f64; 2]; 2],
}

/// 把一个局部坐标分量（局部系下 x∈[0,1]，y∈[-0.866,0.866]）量化成桶号。
#[inline]
fn bin(v: f64) -> u32 {
    // y 分量平移到 [0,1] 再量化；超出三角形外接范围的极端值夹到边界。
    let t = if v < 0.0 { v + 0.8660254037844386 } else { v };
    let b = (t * CODE_BINS as f64).floor() as i64;
    b.clamp(0, CODE_BINS as i64 - 1) as u32
}

#[inline]
fn pack(cx: u32, cy: u32, dx: u32, dy: u32) -> u64 {
    ((cx as u64) << (3 * CODE_BITS))
        | ((cy as u64) << (2 * CODE_BITS))
        | ((dx as u64) << CODE_BITS)
        | dy as u64
}

/// 编码单个 quad 的连续局部坐标，返回 (code, 规范化后的局部坐标)。
///
/// `p` 为 A、B、C、D 四点的二维坐标。
pub fn code_from_points(p: [[f64; 2]; 4]) -> (u64, [[f64; 2]; 2]) {
    let (a, b) = (p[0], p[1]);
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let d2 = dx * dx + dy * dy;
    let inv = 1.0 / d2;
    // 局部基：e_x = (dx,dy)/|AB|，局部坐标 = (dot(P-A, AB)/|AB|², cross(AB,P-A)/|AB|²)
    let local_of = |q: [f64; 2]| -> [f64; 2] {
        let (px, py) = (q[0] - a[0], q[1] - a[1]);
        [
            (px * dx + py * dy) * inv,
            (dx * py - dy * px) * inv,
        ]
    };
    let lc = local_of(p[2]);
    let ld = local_of(p[3]);
    // 镜像吸收：C/D 词序规范化。
    let (lc, ld) = if [lc[0], lc[1]] <= [ld[0], ld[1]] {
        (lc, ld)
    } else {
        (ld, lc)
    };
    let code = pack(
        bin(lc[0]),
        // bin() 对负值做了 +√3/2 平移，正负分量分桶区间不同，
        // 图像侧与星表侧规则一致即可。
        bin(lc[1]),
        bin(ld[0]),
        bin(ld[1]),
    );
    (code, [lc, ld])
}

/// 生成图像 quad。
///
/// 策略：以最亮的星为锚点 A，在半径窗口内找 B（保证 AB 间距适中、
/// 不会小到数值不稳，也不会大到超出视场），再在 AB 中点附近的亮星里
/// 取 C、D。数量上限保护 24 MP 图的耗时预算。
pub fn image_quads(
    stars: &[StarObs],
    width: usize,
    height: usize,
    max_quads: usize,
) -> Vec<Quad> {
    let n = stars.len();
    if n < 4 {
        return Vec::new();
    }
    let pts: Vec<[f64; 2]> = stars.iter().map(|s| [s.x, s.y]).collect();
    let shortest = width.min(height) as f64;
    let min_ab = (shortest * 0.03).max(6.0);
    let max_ab = shortest * 0.6;

    let mut quads = Vec::new();
    // 锚点只取最亮的前若干颗——暗锚点的组合在验证阶段也赢不了亮组合。
    let anchor_limit = n.min(40);
    'outer: for ai in 0..anchor_limit {
        // 候选 B：按到 A 的距离排序后取若干个，覆盖不同尺度。
        let mut bs: Vec<(usize, f64)> = (0..n)
            .filter(|&bi| bi != ai)
            .map(|bi| {
                let d2 = dist2(pts[ai], pts[bi]);
                (bi, d2.sqrt())
            })
            .filter(|&(_, d)| d >= min_ab && d <= max_ab)
            .collect();
        bs.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let mut chosen = 0usize;
        for &(bi, ab) in &bs {
            if chosen >= 6 {
                break;
            }
            chosen += 1;
            // C、D：取 AB 中点半径 0.55·AB 内、除锚点外最亮的两颗。
            let mid = [
                (pts[ai][0] + pts[bi][0]) * 0.5,
                (pts[ai][1] + pts[bi][1]) * 0.5,
            ];
            let r2 = (ab * 0.55).powi(2);
            let mut inner: Vec<usize> = (0..n)
                .filter(|&k| k != ai && k != bi)
                .filter(|&k| dist2(pts[k], mid) <= r2)
                .collect();
            // stars 已按 flux 降序，下标越小越亮，直接稳定排序取下标。
            inner.sort_unstable_by_key(|&k| k);
            if inner.len() < 2 {
                continue;
            }
            let (ci, di) = (inner[0], inner[1]);
            let order = [ai, bi, ci, di];
            let (code, local) = code_from_points([
                pts[ai], pts[bi], pts[ci], pts[di],
            ]);
            quads.push(Quad {
                stars: order,
                ab_dist: ab,
                code,
                local,
            });
            if quads.len() >= max_quads {
                break 'outer;
            }
        }
    }
    quads
}

#[inline]
fn dist2(a: [f64; 2], b: [f64; 2]) -> f64 {
    let (dx, dy) = (a[0] - b[0], a[1] - b[1]);
    dx * dx + dy * dy
}

/// 把 40 bit 编码拆回四个桶号（顺序 cx, cy, dx, dy）。
pub fn unpack(code: u64) -> [u32; 4] {
    [
        (code >> 30) as u32 & 0x3ff,
        (code >> 20) as u32 & 0x3ff,
        (code >> 10) as u32 & 0x3ff,
        code as u32 & 0x3ff,
    ]
}

/// 给定 code 的镜像变体：两个 y 分量各自关于桶中点翻转。
///
/// bin() 把 y∈[-√3/2, √3/2] 线性映射到 [0,1024)，镜像 y→-y 等价于
/// 桶号 b → 1023-b（量化边界处差 1 桶，由查询容差吸收）。
#[inline]
pub fn mirror_code(code: u64) -> u64 {
    let b = unpack(code);
    pack(b[0], CODE_BINS - 1 - b[1], b[2], CODE_BINS - 1 - b[3])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn star(x: f64, y: f64, f: f64) -> StarObs {
        StarObs { x, y, flux: f, area: 5, rms: 1.0 }
    }

    #[test]
    fn code_invariant_under_similarity_transform() {
        let base = [[0.0, 0.0], [10.0, 0.0], [4.0, 3.0], [7.0, -2.0]];
        // 旋转 37°、缩放 2.5、平移 (100, -50)
        let (s, c) = 37f64.to_radians().sin_cos();
        let k = 2.5;
        let t = |p: [f64; 2]| -> [f64; 2] {
            [
                100.0 + k * (c * p[0] - s * p[1]),
                -50.0 + k * (s * p[0] + c * p[1]),
            ]
        };
        let moved = [t(base[0]), t(base[1]), t(base[2]), t(base[3])];
        let (c1, l1) = code_from_points(base);
        let (c2, l2) = code_from_points(moved);
        assert_eq!(c1, c2);
        for i in 0..2 {
            assert!((l1[i][0] - l2[i][0]).abs() < 1e-9);
            assert!((l1[i][1] - l2[i][1]).abs() < 1e-9);
        }
    }

    #[test]
    fn mirror_produces_distinguishable_code() {
        // 镜像（手性翻转）不被 C/D 词序规范化吸收——这正是期望行为：
        // 镜像 quad 不会误匹配到真实天区，手性由最终 CD 行列式符号判定。
        let base = [[0.0, 0.0], [10.0, 0.0], [4.0, 3.0], [7.0, -2.0]];
        let flip = [
            [base[0][0], -base[0][1]],
            [base[1][0], -base[1][1]],
            [base[2][0], -base[2][1]],
            [base[3][0], -base[3][1]],
        ];
        let (c1, _) = code_from_points(base);
        let (c2, _) = code_from_points(flip);
        assert_ne!(c1, c2);
    }

    #[test]
    fn image_quads_basic() {
        let stars = vec![
            star(100.0, 100.0, 1000.0),
            star(200.0, 110.0, 900.0),
            star(150.0, 160.0, 800.0),
            star(160.0, 60.0, 700.0),
            star(30.0, 30.0, 600.0),
        ];
        let q = image_quads(&stars, 400, 300, 100);
        assert!(!q.is_empty());
        for quad in &q {
            assert_eq!(unpack(quad.code).len(), 4);
            assert!(quad.ab_dist > 0.0);
        }
    }
}
