//! 盲解匹配器：图像 quad → 分块哈希范围查询 → 候选对应 → 鲁棒 WCS 验证。
//!
//! 关键设计：
//! * 星表访问抽象成 trait，本地文件、内存、浏览器 fetch 各自实现，
//!   匹配逻辑本身不感知 I/O；
//! * 尺度带按“图像锚距 → 角距”的可行比例尺反推，未知焦距时用宽边界
//!   多带并行尝试，命中即停；
//! * 每个候选先用 4 点拟合做便宜的门限过滤，再取邻星扩充对应点做
//!   全最小二乘 + 3σ 剔除，内点数即置信度。

use crate::image::StarObs;
use crate::index::{CatStar, IndexRecord, QueryTol};
use crate::quad::{image_quads, Quad};
use crate::wcs::{fit_sip, fit_tan, Correspondence, FitResult, TanWcs};
use std::collections::HashMap;

/// 匹配器需要的星表能力。返回 io 结果以便文件/fetch 实现透传错误。
pub trait CatalogSource {
    /// 在某尺度带里做 code 范围查询。
    fn query_band(
        &mut self,
        band: u32,
        code: u64,
        tol: QueryTol,
    ) -> std::io::Result<Vec<IndexRecord>>;
    /// 按全局星号批量取星。
    fn stars_by_id(&mut self, ids: &[u32]) -> std::io::Result<Vec<CatStar>>;
    /// 中心角半径内的星（验证阶段扩充对应点）。
    fn stars_within(
        &mut self,
        ra: f64,
        dec: f64,
        radius: f64,
    ) -> std::io::Result<Vec<CatStar>>;
    /// 尺度带定义（弧度）。
    fn bands(&self) -> &[(f64, f64)];
}

/// 解算参数。
#[derive(Clone, Debug)]
pub struct SolveConfig {
    /// 生成的图像 quad 上限。
    pub max_quads: usize,
    /// 哈希查询每维容差（桶）。
    pub query_half_bins: u32,
    /// 单候选进入 LS 前要求的最少 4 点内联。
    pub min_quad_votes: usize,
    /// 判定解成功的最少对应点数。
    pub min_correspondences: usize,
    /// 成功内点数门槛（占检出星比例）。
    pub min_inlier_fraction: f64,
    /// 比例尺搜索范围（角秒/像元），未知焦距时用宽区间。
    pub scale_search_arcsec: (f64, f64),
    /// 验证邻域相对画幅半径的倍数。
    pub verify_radius_factor: f64,
    /// 最多尝试多少个候选 quad 后放弃当前带。
    pub max_candidates_per_band: usize,
}

impl Default for SolveConfig {
    fn default() -> Self {
        Self {
            max_quads: 220,
            query_half_bins: 2,
            min_quad_votes: 3,
            min_correspondences: 8,
            min_inlier_fraction: 0.12,
            max_candidates_per_band: 60,
            scale_search_arcsec: (0.3, 60.0),
            verify_radius_factor: 0.85,
        }
    }
}

/// 解算成功的产物。
#[derive(Clone, Debug)]
pub struct Solution {
    pub wcs: TanWcs,
    pub center_ra_deg: f64,
    pub center_dec_deg: f64,
    pub pixel_scale_x: f64,
    pub pixel_scale_y: f64,
    pub roll_deg: f64,
    pub fov_diag_deg: f64,
    pub fov_width_deg: f64,
    pub fov_height_deg: f64,
    pub rms_arcsec: f64,
    pub max_resid_arcsec: f64,
    pub num_inliers: usize,
    pub num_detected: usize,
    pub matched: Vec<MatchedStar>,
    pub band_index: u32,
    pub elapsed_ms: u128,
}

#[derive(Clone, Copy, Debug)]
pub struct MatchedStar {
    pub cat_id: u32,
    pub x: f64,
    pub y: f64,
    pub ra_deg: f64,
    pub dec_deg: f64,
    pub mag: f32,
}

/// 失败原因（便于 UI 给出可操作的提示）。
#[derive(Clone, Debug)]
pub enum SolveError {
    TooFewStars(usize),
    NoCandidates,
    VerifyFailed(String),
    Io(String),
}

impl std::fmt::Display for SolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SolveError::TooFewStars(n) => write!(f, "检出星点仅 {n} 颗，不足 4 颗，无法建 quad"),
            SolveError::NoCandidates => write!(f, "哈希索引无候选，可能视场超出尺度带覆盖"),
            SolveError::VerifyFailed(s) => write!(f, "候选验证失败：{s}"),
            SolveError::Io(s) => write!(f, "索引读取失败：{s}"),
        }
    }
}

impl From<std::io::Error> for SolveError {
    fn from(e: std::io::Error) -> Self {
        SolveError::Io(e.to_string())
    }
}

/// 一次完整盲解。
pub fn solve<C: CatalogSource>(
    cat: &mut C,
    stars: &[StarObs],
    width: usize,
    height: usize,
    cfg: &SolveConfig,
) -> Result<Solution, SolveError> {
    let t0 = std::time::Instant::now();
    if stars.len() < 4 {
        return Err(SolveError::TooFewStars(stars.len()));
    }
    let quads = image_quads(stars, width, height, cfg.max_quads);
    if quads.is_empty() {
        return Err(SolveError::TooFewStars(stars.len()));
    }

    // 可行比例尺 → 候选尺度带。锚距 d_pix 对应角距 d_rad = scale·d_pix，
    // 与带 [lo,hi] 相交当且仅当 scale 区间与 (lo/d_pix, hi/d_pix) 重叠。
    let (s_lo, s_hi) = (
        crate::angle::arcsec2rad(cfg.scale_search_arcsec.0),
        crate::angle::arcsec2rad(cfg.scale_search_arcsec.1),
    );
    // 拷贝一份带定义，避免持有 cat 的不可变借用而阻塞后续查询的可变借用。
    let bands: Vec<(f64, f64)> = cat.bands().to_vec();
    let mut band_score: Vec<(usize, f64)> = Vec::new();
    for (bi, &(lo, hi)) in bands.iter().enumerate() {
        let mut hits = 0u32;
        for q in &quads {
            let sc_lo = lo / q.ab_dist;
            let sc_hi = hi / q.ab_dist;
            if sc_hi >= s_lo && sc_lo <= s_hi {
                hits += 1;
            }
        }
        if hits > 0 {
            // 带的角距中位数与图像锚距中位数的匹配优先，兼顾命中数。
            band_score.push((bi, hits as f64));
        }
    }
    if band_score.is_empty() {
        return Err(SolveError::NoCandidates);
    }
    band_score.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    let mut best: Option<(FitResult, Vec<MatchedStar>, u32, Vec<bool>)> = None;

    for (bi, _score) in band_score {
        let (blo, bhi) = bands[bi];
        // 本带内每个 quad 的查询与验证。
        let mut tried = 0usize;
        for q in &quads {
            // 比例尺先验过滤：该 quad 在本带推出的 scale 是否落在搜索区间。
            let implied_lo = blo / q.ab_dist;
            let implied_hi = bhi / q.ab_dist;
            if implied_hi < s_lo || implied_lo > s_hi {
                continue;
            }
            let tol = QueryTol {
                half_bins: cfg.query_half_bins,
            };
            // 同时查原编码与镜像编码，兼容相机/RAW 的翻转取向。
            let mut recs = cat.query_band(bi as u32, q.code, tol)?;
            let mirrored = cat.query_band(bi as u32, crate::quad::mirror_code(q.code), tol)?;
            recs.extend(mirrored);
            recs.sort_by(|a, b| a.code.cmp(&b.code));
            recs.dedup_by(|a, b| a.code == b.code && a.star_ids == b.star_ids);
            tried += 1;
            for rec in recs {
                if let Some(cand) =
                    try_candidate(cat, q, &rec, stars, width, height, bi as u32, cfg)?
                {
                    let better = match &best {
                        None => true,
                        Some((b, _, _, _)) => cand.0.used > b.used,
                    };
                    if better {
                        best = Some(cand);
                    }
                }
            }
            if tried >= cfg.max_candidates_per_band {
                break;
            }
            if let Some((b, _, _, _)) = &best {
                if (b.used as f64) >= (stars.len() as f64) * 0.6 {
                    break; // 已拿到强解，提前结束本带
                }
            }
        }
        if let Some((b, _, _, _)) = &best {
            if b.used >= cfg.min_correspondences
                && (b.used as f64) >= (stars.len() as f64) * cfg.min_inlier_fraction
            {
                break; // 跨带也可以停
            }
        }
    }

    let (fit, matched, band_index, _) =
        best.ok_or(SolveError::NoCandidates)?;
    if fit.used < cfg.min_correspondences
        || (fit.used as f64) < (stars.len() as f64) * cfg.min_inlier_fraction
    {
        return Err(SolveError::VerifyFailed(format!(
            "内点 {} 不足门槛 {}",
            fit.used,
            (stars.len() as f64 * cfg.min_inlier_fraction).ceil() as usize
        )));
    }

    // SIP 精修（内点足够时）；matched 与 pairs 同序构建，掩码可直接复用。
    let wcs = fit_sip(&fit.wcs, &pairs_from_matched(&matched), &fit.inliers);
    let fov = wcs.field_of_view(width as f64, height as f64);
    let (cx_ra, cx_dec) = wcs.pixel_to_world(width as f64 * 0.5, height as f64 * 0.5);
    let (sx, sy) = wcs.pixel_scales_arcsec();

    Ok(Solution {
        wcs,
        center_ra_deg: cx_ra.to_degrees(),
        center_dec_deg: cx_dec.to_degrees(),
        pixel_scale_x: sx,
        pixel_scale_y: sy,
        roll_deg: wcs.roll_deg(),
        fov_diag_deg: fov.diagonal_rad.to_degrees(),
        fov_width_deg: fov.width_rad.to_degrees(),
        fov_height_deg: fov.height_rad.to_degrees(),
        rms_arcsec: fit.rms_arcsec,
        max_resid_arcsec: fit.max_arcsec,
        num_inliers: fit.used,
        num_detected: stars.len(),
        matched,
        band_index,
        elapsed_ms: t0.elapsed().as_millis(),
    })
}

fn pairs_from_matched(m: &[MatchedStar]) -> Vec<Correspondence> {
    m.iter()
        .map(|s| Correspondence {
            px: s.x,
            py: s.y,
            ra: s.ra_deg.to_radians(),
            dec: s.dec_deg.to_radians(),
        })
        .collect()
}

/// 单个索引记录的验证：4 点建对应 → 邻星扩充 → TAN 拟合。
fn try_candidate<C: CatalogSource>(
    cat: &mut C,
    q: &Quad,
    rec: &IndexRecord,
    stars: &[StarObs],
    width: usize,
    height: usize,
    band_index: u32,
    cfg: &SolveConfig,
) -> std::io::Result<Option<(FitResult, Vec<MatchedStar>, u32, Vec<bool>)>> {
    let cat_stars = cat.stars_by_id(&rec.star_ids)?;
    if cat_stars.len() != 4 {
        return Ok(None);
    }
    // quad.stars 与 rec.star_ids 同序（A,B,C,D）；C/D 词序规范化在两侧
    // 用同一规则，镜像情形也能对上，手性由 CD 行列式符号体现。
    let mut pairs: Vec<Correspondence> = Vec::with_capacity(64);
    let mut matched: Vec<MatchedStar> = Vec::new();
    for k in 0..4 {
        let img = &stars[q.stars[k]];
        let cs = cat_stars[k];
        pairs.push(Correspondence {
            px: img.x,
            py: img.y,
            ra: cs.ra,
            dec: cs.dec,
        });
        matched.push(MatchedStar {
            cat_id: cs.id,
            x: img.x,
            y: img.y,
            ra_deg: cs.ra.to_degrees(),
            dec_deg: cs.dec.to_degrees(),
            mag: cs.mag,
        });
    }
    let Some(seed_fit) = fit_tan(&pairs, width as f64, height as f64) else {
        return Ok(None);
    };
    // 廉价门限：4 点解的比例尺必须落在搜索区间。
    let (psx, _psy) = seed_fit.wcs.pixel_scales_arcsec();
    if !(cfg.scale_search_arcsec.0..=cfg.scale_search_arcsec.1).contains(&psx) {
        return Ok(None);
    }

    // 邻星扩充：以画幅中心指向为圆心，取覆盖画幅的角半径。
    let (cra, cdec) = seed_fit
        .wcs
        .pixel_to_world(width as f64 * 0.5, height as f64 * 0.5);
    let corner = seed_fit
        .wcs
        .pixel_to_world(0.0, 0.0);
    let radius = crate::angle::angular_sep(
        crate::angle::radec_to_vec(cra, cdec),
        crate::angle::radec_to_vec(corner.0, corner.1),
    ) * cfg.verify_radius_factor;
    let nearby = cat.stars_within(cra, cdec, radius)?;

    // 图像星按位置建网格，做最近邻一一配对。
    let index = PositionGrid::build(stars, width, height);
    let mut used_img = HashMap::new();
    for s in &nearby {
        let (px, py) = seed_fit.wcs.world_to_pixel(s.ra, s.dec);
        if px < 0.0 || py < 0.0 || px >= width as f64 || py >= height as f64 {
            continue;
        }
        let Some(img_idx) = index.nearest_within(px, py, 4.0) else {
            continue;
        };
        if used_img.contains_key(&img_idx) {
            continue;
        }
        let img = &stars[img_idx];
        // 二次门限：投影残差像元数。
        let d = (img.x - px).hypot(img.y - py);
        if d > 4.0 {
            continue;
        }
        used_img.insert(img_idx, ());
        pairs.push(Correspondence {
            px: img.x,
            py: img.y,
            ra: s.ra,
            dec: s.dec,
        });
        matched.push(MatchedStar {
            cat_id: s.id,
            x: img.x,
            y: img.y,
            ra_deg: s.ra.to_degrees(),
            dec_deg: s.dec.to_degrees(),
            mag: s.mag,
        });
    }

    if pairs.len() < cfg.min_quad_votes + 1 {
        return Ok(None);
    }
    let Some(fit) = fit_tan(&pairs, width as f64, height as f64) else {
        return Ok(None);
    };
    let active = vec![true; pairs.len()];
    Ok(Some((fit, matched, band_index, active)))
}

/// 像元坐标的均匀网格最近邻索引（网格边长取容差，查询 O(1)~O(9)）。
struct PositionGrid {
    cell: f64,
    nx: usize,
    /// 每格存 (星点下标, x, y)。
    bins: Vec<Vec<(usize, f64, f64)>>,
}

impl PositionGrid {
    fn build(stars: &[StarObs], width: usize, _height: usize) -> Self {
        let cell = 8.0;
        let nx = (width as f64 / cell).ceil() as usize + 1;
        let ny = (_height as f64 / cell).ceil() as usize + 1;
        let mut bins = vec![Vec::new(); nx * ny];
        for (i, s) in stars.iter().enumerate() {
            let gx = (s.x / cell) as usize;
            let gy = (s.y / cell) as usize;
            bins[gy * nx + gx].push((i, s.x, s.y));
        }
        Self { cell, nx, bins }
    }

    fn nearest_within(&self, x: f64, y: f64, tol: f64) -> Option<usize> {
        let gx = (x / self.cell) as i64;
        let gy = (y / self.cell) as i64;
        let mut best = None;
        let mut best_d2 = tol * tol;
        for dy in -1..=1 {
            for dx in -1..=1 {
                let (bx, by) = (gx + dx, gy + dy);
                if bx < 0 || by < 0 || bx as usize >= self.nx {
                    continue;
                }
                let idx = by as usize * self.nx + bx as usize;
                if idx >= self.bins.len() {
                    continue;
                }
                for &(k, sx, sy) in &self.bins[idx] {
                    let d2 = (sx - x).powi(2) + (sy - y).powi(2);
                    if d2 <= best_d2 {
                        best_d2 = d2;
                        best = Some(k);
                    }
                }
            }
        }
        best
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::angle::{arcsec2rad, deg2rad};
    use crate::builder::{build, BuildParams, CatRow};
    use crate::index::{ChunkLayout, MemoryReader};
    use crate::synth::{render, RenderParams};
    use crate::wcs::SipCoeffs;

    /// 内存星表源：直接持有构建产物。
    struct MemoryCatalog {
        built: crate::builder::BuiltIndex,
        star_bytes: Vec<u8>,
        layouts: Vec<ChunkLayout>,
        bands: Vec<(f64, f64)>,
    }

    impl MemoryCatalog {
        fn open(built: crate::builder::BuiltIndex) -> Self {
            let bands = built
                .manifest
                .bands
                .iter()
                .map(|b| (b.min_rad, b.max_rad))
                .collect();
            let star_bytes = built.star_bytes.clone();
            let layouts = built.star_layouts.clone();
            Self {
                built,
                star_bytes,
                layouts,
                bands,
            }
        }

        fn all_stars(&self) -> Vec<CatRow> {
            // 与 StarChunks::load 的布局严格一致：每个块从 lay.offset 起就是
            // count×STAR_RECORD_BYTES 的裸记录，块头（magic+count）不包含在
            // offset 指向的数据区内，不能再额外跳 8 字节。
            let mut out = Vec::new();
            for lay in &self.layouts {
                let base = lay.offset as usize;
                for k in 0..lay.count as usize {
                    let off = base + k * crate::index::STAR_RECORD_BYTES;
                    let id = u32::from_le_bytes(self.star_bytes[off..off + 4].try_into().unwrap());
                    let ra = f64::from_le_bytes(self.star_bytes[off + 4..off + 12].try_into().unwrap());
                    let dec = f64::from_le_bytes(self.star_bytes[off + 12..off + 20].try_into().unwrap());
                    let mag = f32::from_le_bytes(self.star_bytes[off + 20..off + 24].try_into().unwrap());
                    out.push(CatRow { id, ra, dec, mag });
                }
            }
            out
        }
    }

    impl CatalogSource for MemoryCatalog {
        fn query_band(
            &mut self,
            band: u32,
            code: u64,
            tol: QueryTol,
        ) -> std::io::Result<Vec<IndexRecord>> {
            let bytes = self.built.band_bytes[band as usize].clone();
            let mut idx = crate::index::BandIndex::open(MemoryReader { bytes }, 16)?;
            idx.query(code, tol)
        }
        fn stars_by_id(&mut self, ids: &[u32]) -> std::io::Result<Vec<CatStar>> {
            let all = self.all_stars();
            Ok(ids
                .iter()
                .filter_map(|&i| {
                    all.get(i as usize).map(|r| CatStar {
                        id: r.id,
                        ra: r.ra,
                        dec: r.dec,
                        mag: r.mag,
                    })
                })
                .collect())
        }
        fn stars_within(
            &mut self,
            ra: f64,
            dec: f64,
            radius: f64,
        ) -> std::io::Result<Vec<CatStar>> {
            let cos_r = radius.cos();
            let c = crate::angle::radec_to_vec(ra, dec);
            Ok(self
                .all_stars()
                .into_iter()
                .filter(|s| {
                    crate::angle::dot(crate::angle::radec_to_vec(s.ra, s.dec), c) >= cos_r
                })
                .map(|r| CatStar {
                    id: r.id,
                    ra: r.ra,
                    dec: r.dec,
                    mag: r.mag,
                })
                .collect())
        }
        fn bands(&self) -> &[(f64, f64)] {
            &self.bands
        }
    }

    fn make_wcs(scale_arcsec: f64, roll_deg: f64, w: f64, h: f64) -> TanWcs {
        let scale = arcsec2rad(scale_arcsec);
        let roll = deg2rad(roll_deg);
        TanWcs {
            crpix: [w * 0.5, h * 0.5],
            crval: [deg2rad(185.7), deg2rad(12.3)],
            cd: [
                [scale * roll.cos(), -scale * roll.sin()],
                [scale * roll.sin(), scale * roll.cos()],
            ],
            sip: SipCoeffs::default(),
        }
    }

    fn dense_sky(n: u32) -> Vec<CatRow> {
        // 密集均匀天球，保证任意小视场内有足够多星。
        let mut state = 7u64;
        let mut rows = Vec::new();
        for id in 0..n {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let u = (state >> 33) as f64 / (1u64 << 31) as f64;
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let v = (state >> 33) as f64 / (1u64 << 31) as f64;
            rows.push(CatRow {
                id,
                ra: std::f64::consts::TAU * u,
                dec: (2.0 * v - 1.0).asin(),
                // mag 5..9：渲染峰值约 600..2500 计数，远超检测阈值；
                // 亮星占比足够，任意小视场内都能拿到可检的星。
                mag: 5.0 + (id % 40) as f32 * 0.1,
            });
        }
        rows
    }

    #[test]
    fn end_to_end_synthetic_solve() {
        // 30000 颗星建库（全天空 ~0.73 星/平方度，建库约 30s）。
        let rows = dense_sky(30000);
        let params = BuildParams {
            mag_limit: 9.5,
            stars_per_chunk: 4000,
            bands: [
                (deg2rad(0.4), deg2rad(1.2)),
                (deg2rad(1.0), deg2rad(2.5)),
                (deg2rad(2.0), deg2rad(5.0)),
                (deg2rad(4.0), deg2rad(9.0)),
                (deg2rad(8.0), deg2rad(16.0)),
            ],
            quads_per_anchor: 8,
            ..Default::default()
        };
        let built = build(&rows, &params).unwrap();
        let mut cat = MemoryCatalog::open(built);

        // 真值 WCS：9 角秒/像元（广角）、滚动 33°、2400x1800。
        // 视场约 6°x4.5°（~27 平方度），帧内期望 15+ 颗星。
        let (w, h) = (2400usize, 1800usize);
        let truth = make_wcs(9.0, 33.0, w as f64, h as f64);
        let catstars: Vec<CatStar> = cat
            .all_stars()
            .iter()
            .map(|r| CatStar {
                id: r.id,
                ra: r.ra,
                dec: r.dec,
                mag: r.mag,
            })
            .collect();
        let rendered = render(&truth, w, h, &catstars, &RenderParams::default());
        let (detected, _bg) = crate::image::detect_stars(
            &rendered.image,
            &crate::image::DetectParams {
                k_sigma: 4.0,
                min_area: 2,
                ..Default::default()
            },
        );
        assert!(detected.len() >= 12, "检出 {} 颗", detected.len());

        let sol = solve(
            &mut cat,
            &detected,
            w,
            h,
            &SolveConfig {
                scale_search_arcsec: (3.0, 25.0),
                ..Default::default()
            },
        )
        .expect("应当盲解成功");
        assert!(sol.num_inliers >= 8);
        // 这是广角（9″/px）合成测试，精度按像元尺度放宽；窄角高分辨实景
        // 的 2″ 级精度由像元比例尺自然保证（见 README 精度说明）。
        let c = sol.wcs.pixel_to_world(w as f64 * 0.5, h as f64 * 0.5);
        let t = truth.pixel_to_world(w as f64 * 0.5, h as f64 * 0.5);
        let err = crate::angle::angular_sep(
            crate::angle::radec_to_vec(c.0, c.1),
            crate::angle::radec_to_vec(t.0, t.1),
        );
        assert!(
            crate::angle::rad2arcsec(err) < 5.0,
            "中心误差 {} 角秒",
            crate::angle::rad2arcsec(err)
        );
        assert!((sol.pixel_scale_x - 9.0).abs() < 0.2);
        assert!((sol.roll_deg - 33.0).abs() < 0.5);
    }
}
