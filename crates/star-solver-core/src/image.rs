//! 星点提取：背景/噪声估计 → 阈值分割 → 连通域标记 → 流量加权亚像素质心。
//!
//! 输入是 8 bit 灰度（浏览器端用 Canvas 解码 JPEG/PNG 后下采样传入；
//! CLI 端也走同一接口）。内部全用 f64 累积，避免亮星溢出。

use crate::angle::median;

/// 单通道灰度图（行优先，`data[y*width + x]`，取值 0..=255）。
#[derive(Clone)]
pub struct GrayImage {
    pub width: usize,
    pub height: usize,
    pub data: Vec<u8>,
}

impl GrayImage {
    pub fn new(width: usize, height: usize, data: Vec<u8>) -> Self {
        assert_eq!(data.len(), width * height);
        Self { width, height, data }
    }

    #[inline]
    pub fn pixel(&self, x: usize, y: usize) -> u8 {
        self.data[y * self.width + x]
    }

    /// 按比例 `scale`（输出边长 = 输入边长 / scale）做面积平均降采样。
    /// 降采样同时是抗混叠，亮星质心精度损失远小于单像素尺度。
    pub fn downsample(&self, scale: usize) -> GrayImage {
        if scale <= 1 {
            return self.clone();
        }
        let (ow, oh) = (self.width / scale, self.height / scale);
        let mut out = vec![0u8; ow * oh];
        let area = (scale * scale) as f64;
        for oy in 0..oh {
            for ox in 0..ow {
                let mut sum = 0u32;
                for dy in 0..scale {
                    let row = (oy * scale + dy) * self.width;
                    for dx in 0..scale {
                        sum += self.data[row + ox * scale + dx] as u32;
                    }
                }
                out[oy * ow + ox] = (sum as f64 / area).round() as u8;
            }
        }
        GrayImage::new(ow, oh, out)
    }
}

/// 星点检测参数。默认值针对普通数码/冷冻相机的深空图调过。
#[derive(Clone, Debug)]
pub struct DetectParams {
    /// 阈值 = 背景 + k_sigma × 噪声估计。
    pub k_sigma: f64,
    /// 连通域最小像素数（滤热噪单像素）。
    pub min_area: usize,
    /// 连通域最大像素数（滤过曝大斑、云迹）。
    pub max_area: usize,
    /// 长短轴比上限（圆度筛选，1.0=完美圆；过曝星拉长时可放宽）。
    pub max_elongation: f64,
    /// 返回流量最高的星数（后续 quad 只需要最亮的一小批）。
    pub max_stars: usize,
}

impl Default for DetectParams {
    fn default() -> Self {
        Self {
            k_sigma: 5.0,
            min_area: 3,
            max_area: 400,
            max_elongation: 2.5,
            max_stars: 500,
        }
    }
}

/// 检测到的星点。位置为亚像素像元坐标（原点在像元 0 中心）。
#[derive(Clone, Debug)]
pub struct StarObs {
    pub x: f64,
    pub y: f64,
    /// 背景扣除后的总计数（流量，跨图可比，无绝对定标）。
    pub flux: f64,
    /// 连通域面积（像素数）。
    pub area: usize,
    /// 流量加权质心的标准差（像元），用于评估质心质量。
    pub rms: f64,
}

/// 背景模型：`level + n(0, sigma)`。
#[derive(Clone, Copy, Debug)]
pub struct Background {
    pub level: f64,
    pub sigma: f64,
}

/// 从整幅图的像素分布估计背景（中位数 + 归一化 MAD）。
///
/// 深空图里背景像素占绝大多数，中位数稳健；MAD 对亮星不敏感，
/// 0.6745 是把 MAD 换算成高斯标准差的归一化常数。
pub fn estimate_background(img: &GrayImage) -> Background {
    let mut sorted: Vec<f64> = img.data.iter().map(|&p| p as f64).collect();
    sorted.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let level = median(&sorted);
    let mut dev: Vec<f64> = img
        .data
        .iter()
        .map(|&p| (p as f64 - level).abs())
        .collect();
    dev.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let sigma = (median(&dev) * 1.4826).max(0.5);
    Background { level, sigma }
}

/// 并查集（路径压缩 + 按秩合并），用于连通域标记。
struct UnionFind {
    parent: Vec<u32>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n as u32).collect(),
            rank: vec![0; n],
        }
    }
    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            self.parent[x as usize] = self.parent[self.parent[x as usize] as usize];
            x = self.parent[x as usize];
        }
        x
    }
    fn union(&mut self, a: u32, b: u32) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        match self.rank[ra as usize].cmp(&self.rank[rb as usize]) {
            std::cmp::Ordering::Less => self.parent[ra as usize] = rb,
            std::cmp::Ordering::Greater => self.parent[rb as usize] = ra,
            std::cmp::Ordering::Equal => {
                self.parent[rb as usize] = ra;
                self.rank[ra as usize] += 1;
            }
        }
    }
}

/// 连通域累积量（零阶/一阶/二阶矩，f64 防止 24 MP 亮星 u32 溢出）。
#[derive(Default)]
struct Component {
    area: u32,
    sum_w: f64,
    sum_wx: f64,
    sum_wy: f64,
    sum_wxx: f64,
    sum_wyy: f64,
}

/// 全流程星点检测。
pub fn detect_stars(img: &GrayImage, params: &DetectParams) -> (Vec<StarObs>, Background) {
    let bg = estimate_background(img);
    let threshold = bg.level + params.k_sigma * bg.sigma;
    let (w, h) = (img.width, img.height);
    let npix = w * h;

    // 第一遍：阈值像素打临时标号，4-邻接并查集合并。
    let mut label = vec![u32::MAX; npix];
    let mut uf = UnionFind::new(npix);
    let mut next: u32 = 0;
    for y in 0..h {
        let row = y * w;
        for x in 0..w {
            let idx = row + x;
            if (img.data[idx] as f64) < threshold {
                continue;
            }
            let mut root: Option<u32> = None;
            if x > 0 && label[idx - 1] != u32::MAX {
                root = Some(uf.find(label[idx - 1]));
            }
            if y > 0 && label[idx - w] != u32::MAX {
                let up = uf.find(label[idx - w]);
                root = Some(match root {
                    Some(r) => {
                        uf.union(r, up);
                        uf.find(r)
                    }
                    None => up,
                });
            }
            match root {
                Some(r) => label[idx] = r,
                None => {
                    label[idx] = next;
                    next += 1;
                }
            }
        }
    }

    // 第二遍：按根标号聚合矩。
    let mut comps: Vec<Component> = (0..next as usize).map(|_| Component::default()).collect();
    for y in 0..h {
        let row = y * w;
        for x in 0..w {
            let idx = row + x;
            if label[idx] == u32::MAX {
                continue;
            }
            let root = uf.find(label[idx]) as usize;
            let weight = (img.data[idx] as f64 - bg.level).max(0.0);
            let (fx, fy) = (x as f64, y as f64);
            let c = &mut comps[root];
            c.area += 1;
            c.sum_w += weight;
            c.sum_wx += weight * fx;
            c.sum_wyy += weight * fy * fy;
            c.sum_wxx += weight * fx * fx;
            c.sum_wy += weight * fy;
        }
    }

    // 第三遍：形态筛选 + 流量加权质心。
    let mut stars = Vec::new();
    for c in comps {
        if c.area == 0 {
            continue; // 被合并掉的非根槽位
        }
        let area = c.area as usize;
        if area < params.min_area || area > params.max_area || c.sum_w <= 0.0 {
            continue;
        }
        let (mx, my) = (c.sum_wx / c.sum_w, c.sum_wy / c.sum_w);
        let var_x = (c.sum_wxx / c.sum_w - mx * mx).max(0.0);
        let var_y = (c.sum_wyy / c.sum_w - my * my).max(0.0);
        let (sx, sy) = (var_x.sqrt(), var_y.sqrt());
        let long = sx.max(sy);
        let short = sx.min(sy).max(1e-6);
        if long / short > params.max_elongation {
            continue;
        }
        let rms = ((var_x + var_y) * 0.5).sqrt();
        stars.push(StarObs {
            x: mx,
            y: my,
            flux: c.sum_w,
            area,
            rms,
        });
    }

    stars.sort_unstable_by(|a, b| b.flux.partial_cmp(&a.flux).unwrap());
    stars.truncate(params.max_stars);
    (stars, bg)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn background_of_flat_field() {
        let img = GrayImage::new(16, 16, vec![100u8; 256]);
        let bg = estimate_background(&img);
        assert!((bg.level - 100.0).abs() < 1e-9);
    }

    #[test]
    fn detects_synthetic_star_at_subpixel_center() {
        // 32x32 暗场，中心 (12.3, 7.8) 放一个高斯星斑。
        let (w, h) = (32usize, 32usize);
        let mut data = vec![20u8; w * h];
        let (cx, cy, sigma) = (12.3f64, 7.8f64, 1.2f64);
        for y in 0..h {
            for x in 0..w {
                let r2 = (x as f64 - cx).powi(2) + (y as f64 - cy).powi(2);
                let v = 20.0 + 200.0 * (-r2 / (2.0 * sigma * sigma)).exp();
                data[y * w + x] = v.round() as u8;
            }
        }
        let img = GrayImage::new(w, h, data);
        let (stars, _) = detect_stars(
            &img,
            &DetectParams {
                min_area: 2,
                ..Default::default()
            },
        );
        assert_eq!(stars.len(), 1, "应只检出一颗星");
        let s = &stars[0];
        assert!((s.x - cx).abs() < 0.05, "x 质心偏差 {:.4} 像元", s.x - cx);
        assert!((s.y - cy).abs() < 0.05, "y 质心偏差 {:.4} 像元", s.y - cy);
    }

    #[test]
    fn rejects_hot_pixel_and_elongated_trail() {
        let (w, h) = (40usize, 20usize);
        let mut data = vec![10u8; w * h];
        // 单像素热点（min_area=3 时应被滤掉）
        data[2 * w + 3] = 255;
        // 横向拖线
        for x in 5..15 {
            data[10 * w + x] = 250;
        }
        let img = GrayImage::new(w, h, data);
        let (stars, _) = detect_stars(&img, &DetectParams::default());
        assert!(stars.is_empty(), "热点和拖线都应被形态过滤，实际 {:?}",
                stars.iter().map(|s| (s.x, s.y)).collect::<Vec<_>>());
    }
}
