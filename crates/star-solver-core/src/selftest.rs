//! 无需外部数据的端到端自检：合成星表 → 建库 → 渲染 → 检测 → 盲解 → 校验。
//!
//! CLI 的 `selftest` 与 wasm 的健康检查共用，全部确定性、可离线运行。

use crate::angle::{arcsec2rad, deg2rad, radec_to_vec, angular_sep, rad2arcsec};
use crate::builder::{build, BuildParams, CatRow};
use crate::image::{detect_stars, DetectParams};
use crate::index::{CatStar, IndexRecord, MemoryReader, QueryTol};
use crate::matcher::{solve, CatalogSource, SolveConfig};
use crate::synth::{render, RenderParams};
use crate::wcs::{SipCoeffs, TanWcs};

/// 内存星表源（自检与小型演示用）。
pub struct MemoryCatalog {
    built: crate::builder::BuiltIndex,
    layouts: Vec<crate::index::ChunkLayout>,
    bands: Vec<(f64, f64)>,
}

impl MemoryCatalog {
    pub fn open(built: crate::builder::BuiltIndex) -> Self {
        let bands = built
            .manifest
            .bands
            .iter()
            .map(|b| (b.min_rad, b.max_rad))
            .collect();
        let layouts = built.star_layouts.clone();
        Self { built, layouts, bands }
    }

    /// 按 ChunkLayout 从裸记录区还原全部星（offset 已跳过块头）。
    pub fn all_stars(&self) -> Vec<CatRow> {
        let mut out = Vec::new();
        for lay in &self.layouts {
            let base = lay.offset as usize;
            for k in 0..lay.count as usize {
                let off = base + k * crate::index::STAR_RECORD_BYTES;
                let b = &self.built.star_bytes;
                let id = u32::from_le_bytes(b[off..off + 4].try_into().unwrap());
                let ra = f64::from_le_bytes(b[off + 4..off + 12].try_into().unwrap());
                let dec = f64::from_le_bytes(b[off + 12..off + 20].try_into().unwrap());
                let mag = f32::from_le_bytes(b[off + 20..off + 24].try_into().unwrap());
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
        let c = radec_to_vec(ra, dec);
        Ok(self
            .all_stars()
            .into_iter()
            .filter(|s| crate::angle::dot(radec_to_vec(s.ra, s.dec), c) >= cos_r)
            .map(|r| CatStar { id: r.id, ra: r.ra, dec: r.dec, mag: r.mag })
            .collect())
    }

    fn bands(&self) -> &[(f64, f64)] {
        &self.bands
    }
}

fn dense_sky(n: u32) -> Vec<CatRow> {
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
            mag: 5.0 + (id % 40) as f32 * 0.1,
        });
    }
    rows
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

/// 自检结果。
pub struct Report {
    pub detected: usize,
    pub inliers: usize,
    pub center_error_arcsec: f64,
    pub scale_error_pct: f64,
    pub roll_error_deg: f64,
    pub rms_arcsec: f64,
}

/// 跑一次确定性端到端流水线，失败返回 Err（含原因）。
pub fn run() {
    match try_run() {
        Ok(r) => {
            println!(
                "自检通过：检出 {} 颗、内点 {} 颗、中心误差 {:.3}″、比例尺偏差 {:.3}%、滚动角偏差 {:.4}°、RMS {:.3}″",
                r.detected, r.inliers, r.center_error_arcsec, r.scale_error_pct,
                r.roll_error_deg, r.rms_arcsec
            );
        }
        Err(e) => {
            eprintln!("自检失败：{e}");
            std::process::exit(1);
        }
    }
}

/// 供 wasm/外部调用的版本：返回结构化报告而非直接退出进程。
pub fn try_run_pub() -> Result<Report, String> {
    try_run()
}

fn try_run() -> Result<Report, String> {
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
    let built = build(&rows, &params).map_err(|e| e.to_string())?;
    let mut cat = MemoryCatalog::open(built);

    let (w, h) = (2400usize, 1800usize);
    let truth = make_wcs(9.0, 33.0, w as f64, h as f64);
    let catstars: Vec<CatStar> = cat
        .all_stars()
        .iter()
        .map(|r| CatStar { id: r.id, ra: r.ra, dec: r.dec, mag: r.mag })
        .collect();
    let rendered = render(&truth, w, h, &catstars, &RenderParams::default());
    let (detected, _) = detect_stars(
        &rendered.image,
        &DetectParams { k_sigma: 4.0, min_area: 2, ..Default::default() },
    );
    if detected.len() < 12 {
        return Err(format!("仅检出 {} 颗星", detected.len()));
    }

    let sol = solve(
        &mut cat,
        &detected,
        w,
        h,
        &SolveConfig { scale_search_arcsec: (3.0, 25.0), ..Default::default() },
    )
    .map_err(|e| e.to_string())?;

    let c = sol.wcs.pixel_to_world(w as f64 * 0.5, h as f64 * 0.5);
    let t = truth.pixel_to_world(w as f64 * 0.5, h as f64 * 0.5);
    let center_err = rad2arcsec(angular_sep(radec_to_vec(c.0, c.1), radec_to_vec(t.0, t.1)));
    let scale_err = (sol.pixel_scale_x - 9.0).abs() / 9.0 * 100.0;
    let roll_err = (sol.roll_deg - 33.0).abs();

    Ok(Report {
        detected: detected.len(),
        inliers: sol.num_inliers,
        center_error_arcsec: center_err,
        scale_error_pct: scale_err,
        roll_error_deg: roll_err,
        rms_arcsec: sol.rms_arcsec,
    })
}
