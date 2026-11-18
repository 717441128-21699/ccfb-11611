//! star-solver 的 WebAssembly 绑定。
//!
//! 设计原则：所有重活（检测、哈希、匹配、WCS）在 Rust 内完成；JS 侧只负责
//! 图像解码（Canvas/ImageBitmap）、把灰度字节与尺寸传进来，以及承载索引块的
//! 按需读取（fetch + Range，见 web 包的 BlockReader 适配）。
//!
//! 本文件提供与具体 I/O 解耦的同步 API：JS 把需要的索引/星表字节一次性或
//! 分块喂入，Rust 在内存中解算。零网络依赖，满足离线要求。

use star_solver_core::builder::{build, parse_catalog_csv, BuildParams, CatRow, BuiltIndex};
use star_solver_core::image::{detect_stars, estimate_background, DetectParams, GrayImage};
use star_solver_core::index::{CatStar, IndexRecord, MemoryReader, QueryTol};
use star_solver_core::matcher::{solve, SolveConfig, Solution};
use star_solver_core::quad::mirror_code;
use star_solver_core::selftest::MemoryCatalog;
use star_solver_core::wcs::{SipCoeffs, TanWcs};
use wasm_bindgen::prelude::*;

/// 检测到的星点（JS 友好的扁平结构）。
#[wasm_bindgen]
pub struct DetectedStar {
    pub x: f64,
    pub y: f64,
    pub flux: f64,
    pub area: u32,
}

#[wasm_bindgen]
pub struct BackgroundStats {
    pub level: f64,
    pub sigma: f64,
}

/// 星点检测入口：灰度字节（行优先，0..255）+ 宽高。
#[wasm_bindgen]
pub fn detect(
    gray: &[u8],
    width: usize,
    height: usize,
    k_sigma: f64,
    max_stars: usize,
) -> Result<Vec<DetectedStar>, JsValue> {
    if gray.len() != width * height {
        return Err(JsValue::from_str("灰度字节数与宽高不符"));
    }
    let img = GrayImage::new(width, height, gray.to_vec());
    let params = DetectParams {
        k_sigma: if k_sigma > 0.0 { k_sigma } else { 5.0 },
        max_stars: if max_stars > 0 { max_stars } else { 500 },
        ..Default::default()
    };
    let (stars, _) = detect_stars(&img, &params);
    Ok(stars
        .into_iter()
        .map(|s| DetectedStar {
            x: s.x,
            y: s.y,
            flux: s.flux,
            area: s.area as u32,
        })
        .collect())
}

/// 背景估计（供 UI 显示与手动调阈值）。
#[wasm_bindgen]
pub fn background(gray: &[u8], width: usize, height: usize) -> Result<BackgroundStats, JsValue> {
    if gray.len() != width * height {
        return Err(JsValue::from_str("灰度字节数与宽高不符"));
    }
    let img = GrayImage::new(width, height, gray.to_vec());
    let bg = estimate_background(&img);
    Ok(BackgroundStats { level: bg.level, sigma: bg.sigma })
}

/// 由 CSV 星表文本在内存中构建索引（小星表/演示用，返回各文件的字节表）。
#[wasm_bindgen]
pub fn build_index(csv_text: &str, mag_limit: f64) -> Result<WasmIndex, JsValue> {
    let rows = parse_catalog_csv(csv_text).map_err(|e| JsValue::from_str(&e))?;
    let params = BuildParams {
        mag_limit: if mag_limit > 0.0 { mag_limit } else { 12.0 },
        ..Default::default()
    };
    let built = build(&rows, &params).map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(WasmIndex { built })
}

/// 构建产物的持有者（不透明句柄，JS 不解码内部布局）。
#[wasm_bindgen]
pub struct WasmIndex {
    built: BuiltIndex,
}

#[wasm_bindgen]
impl WasmIndex {
    #[wasm_bindgen(getter)]
    pub fn star_count(&self) -> f64 {
        self.built.manifest.star_count as f64
    }

    #[wasm_bindgen(getter)]
    pub fn band_count(&self) -> usize {
        self.built.manifest.bands.len()
    }

    /// manifest 的 JSON 文本（供 UI 展示/下载）。
    pub fn manifest_json(&self) -> String {
        star_solver_core::index::manifest_json::to_json(
            &self.built.manifest,
            &self.built.star_layouts,
        )
    }

    /// 合并后的星表块字节（含块头，可直接落盘为 stars.bin）。
    pub fn star_bytes(&self) -> Vec<u8> {
        self.built.star_bytes.clone()
    }

    /// 第 band 个尺度带文件字节。
    pub fn band_bytes(&self, band: usize) -> Option<Vec<u8>> {
        self.built.band_bytes.get(band).cloned()
    }
}

/// 完整解算的 JS 结果。
#[wasm_bindgen]
pub struct SolveResult {
    ok: bool,
    error: String,
    json: String,
}

#[wasm_bindgen]
impl SolveResult {
    #[wasm_bindgen(getter)]
    pub fn ok(&self) -> bool {
        self.ok
    }
    #[wasm_bindgen(getter)]
    pub fn error(&self) -> String {
        self.error.clone()
    }
    /// 成功时的完整 JSON（中心、比例尺、滚动、视场、内点、FITS 头、匹配星）。
    #[wasm_bindgen(getter)]
    pub fn json(&self) -> String {
        self.json.clone()
    }
}

/// 用内存索引对一张灰度图做盲解。
///
/// `index` 为 [`build_index`] 返回的句柄；`scale_lo/hi` 是先验比例尺范围
/// （角秒/像元），未知焦距时给宽区间（如 0.3~60）。
#[wasm_bindgen]
pub fn solve_image(
    index: &WasmIndex,
    gray: &[u8],
    width: usize,
    height: usize,
    scale_lo: f64,
    scale_hi: f64,
) -> SolveResult {
    let run = || -> Result<Solution, String> {
        if gray.len() != width * height {
            return Err("灰度字节数与宽高不符".into());
        }
        // BuiltIndex 被 WasmIndex 借用，这里克隆出一份独立的内存星表。
        let cloned = clone_built(&index.built);
        let mut cat = MemoryCatalog::open(cloned);
        let img = GrayImage::new(width, height, gray.to_vec());
        let (stars, _) = detect_stars(&img, &DetectParams::default());
        let cfg = SolveConfig {
            scale_search_arcsec: (
                if scale_lo > 0.0 { scale_lo } else { 0.3 },
                if scale_hi > 0.0 { scale_hi } else { 60.0 },
            ),
            ..Default::default()
        };
        solve(&mut cat, &stars, width, height, &cfg).map_err(|e| e.to_string())
    };
    match run() {
        Ok(sol) => SolveResult {
            ok: true,
            error: String::new(),
            json: solution_json(&sol),
        },
        Err(e) => SolveResult {
            ok: false,
            error: e,
            json: String::new(),
        },
    }
}

/// 内置合成星空自检（无需任何输入数据）。
#[wasm_bindgen]
pub fn run_selftest() -> SolveResult {
    match star_solver_core::selftest::try_run_pub() {
        Ok(report) => SolveResult {
            ok: true,
            error: String::new(),
            json: serde_lite::stringify(&report),
        },
        Err(e) => SolveResult {
            ok: false,
            error: e,
            json: String::new(),
        },
    }
}

/// 暴露镜像编码（JS 侧若要自行实现查询可用）。
#[wasm_bindgen]
pub fn mirror_code_of(code: f64) -> f64 {
    mirror_code(code as u64) as f64
}

/// 用解出的 WCS 把 DSO 表投影到画幅，返回落在图内天体的 JSON 数组字符串。
///
/// 入参 `wcs_json` 取自 solve_image 结果里的 `wcs` 字段；`dso_csv` 为
/// `name,ra_deg,dec_deg,kind,mag,major_arcmin,minor_arcmin,pa_deg` 格式。
/// 输出元素形如 `{"name","kind","x","y","ra_deg","dec_deg","mag"}`。
#[wasm_bindgen]
pub fn project_dso(
    wcs_json: &str,
    dso_csv: &str,
    width: f64,
    height: f64,
) -> Result<String, JsValue> {
    use star_solver_core::dso::{dso_in_frame, parse_dso_csv};
    use star_solver_core::json::json_escape;

    let w = parse_wcs_json(wcs_json)?;
    let objects = parse_dso_csv(dso_csv).map_err(|e| JsValue::from_str(&e))?;
    let plotted = dso_in_frame(&w, &objects, width, height, 8.0);
    let items: Vec<String> = plotted
        .iter()
        .map(|p| {
            format!(
                "{{\"name\":\"{}\",\"kind\":\"{}\",\"x\":{:.3},\"y\":{:.3},\
                 \"ra_deg\":{:.7},\"dec_deg\":{:.7},\"mag\":{:.3}}}",
                json_escape(&p.dso.name),
                p.dso.kind.code(),
                p.x,
                p.y,
                p.dso.ra.to_degrees(),
                p.dso.dec.to_degrees(),
                p.dso.mag
            )
        })
        .collect();
    Ok(format!("[{}]", items.join(",")))
}

/// 从 solve_image 结果的 wcs 子对象重建 TanWcs（SIP 留空，DSO 标注用线性 TAN 足够）。
fn parse_wcs_json(text: &str) -> Result<TanWcs, JsValue> {
    let v = star_solver_core::json::parse(text).map_err(|e| JsValue::from_str(&e))?;
    let getf = |k: &str| -> Result<f64, JsValue> {
        v[k]
            .as_f64()
            .ok_or_else(|| JsValue::from_str(&format!("wcs 缺字段 {k}")))
    };
    Ok(TanWcs {
        crpix: [getf("crpix1")?, getf("crpix2")?],
        crval: [
            getf("crval1_deg")?.to_radians(),
            getf("crval2_deg")?.to_radians(),
        ],
        cd: [
            [getf("cd11")?, getf("cd12")?],
            [getf("cd21")?, getf("cd22")?],
        ],
        sip: SipCoeffs::default(),
    })
}

// ---- 内部辅助 ----

/// BuiltIndex 没有派生 Clone（含大 Vec），这里按字段显式克隆。
fn clone_built(b: &BuiltIndex) -> BuiltIndex {
    BuiltIndex {
        manifest: b.manifest.clone(),
        star_layouts: b.star_layouts.clone(),
        star_bytes: b.star_bytes.clone(),
        band_bytes: b.band_bytes.clone(),
    }
}

fn solution_json(sol: &Solution) -> String {
    use star_solver_core::json::json_escape;
    let matched: Vec<String> = sol
        .matched
        .iter()
        .map(|m| {
            format!(
                "{{\"cat_id\":{},\"x\":{:.4},\"y\":{:.4},\"ra_deg\":{:.7},\"dec_deg\":{:.7},\"mag\":{:.3}}}",
                m.cat_id, m.x, m.y, m.ra_deg, m.dec_deg, m.mag
            )
        })
        .collect();
    format!(
        "{{\"center_ra_deg\":{:.7},\"center_dec_deg\":{:.7},\
         \"pixel_scale_x\":{:.4},\"pixel_scale_y\":{:.4},\"roll_deg\":{:.4},\
         \"fov_width_deg\":{:.5},\"fov_height_deg\":{:.5},\"fov_diag_deg\":{:.5},\
         \"rms_arcsec\":{:.4},\"max_resid_arcsec\":{:.4},\
         \"num_inliers\":{},\"num_detected\":{},\"elapsed_ms\":{},\
         \"wcs\":{{\"crpix1\":{:.6},\"crpix2\":{:.6},\"crval1_deg\":{:.10},\
         \"crval2_deg\":{:.10},\"cd11\":{:.12e},\"cd12\":{:.12e},\
         \"cd21\":{:.12e},\"cd22\":{:.12e}}},\
         \"fits_header\":\"{}\",\"matched\":[{}]}}",
        sol.center_ra_deg,
        sol.center_dec_deg,
        sol.pixel_scale_x,
        sol.pixel_scale_y,
        sol.roll_deg,
        sol.fov_width_deg,
        sol.fov_height_deg,
        sol.fov_diag_deg,
        sol.rms_arcsec,
        sol.max_resid_arcsec,
        sol.num_inliers,
        sol.num_detected,
        sol.elapsed_ms,
        sol.wcs.crpix[0],
        sol.wcs.crpix[1],
        sol.wcs.crval[0].to_degrees(),
        sol.wcs.crval[1].to_degrees(),
        sol.wcs.cd[0][0],
        sol.wcs.cd[0][1],
        sol.wcs.cd[1][0],
        sol.wcs.cd[1][1],
        json_escape(&sol.wcs.to_fits_header()).replace('\n', "\\n"),
        matched.join(",")
    )
}

/// 极简 JSON 序列化（避免给 wasm crate 引 serde）。
mod serde_lite {
    pub fn stringify(r: &star_solver_core::selftest::Report) -> String {
        format!(
            "{{\"detected\":{},\"inliers\":{},\"center_error_arcsec\":{:.4},\
             \"scale_error_pct\":{:.4},\"roll_error_deg\":{:.5},\"rms_arcsec\":{:.4}}}",
            r.detected,
            r.inliers,
            r.center_error_arcsec,
            r.scale_error_pct,
            r.roll_error_deg,
            r.rms_arcsec
        )
    }
}

// 让未直接使用的导入保持在文档中（避免警告噪音）。
#[allow(dead_code)]
fn _unused(_: CatRow, _: CatStar, _: IndexRecord, _: QueryTol, _: MemoryReader) {}
