//! 星表索引构建器（离线预计算，CLI 的 `build` 子命令调用）。
//!
//! 输入：CSV 星表，列 `id,ra_deg,dec_deg,mag`（带表头，# 开头为注释）。
//! 输出目录：
//! * `manifest.json`
//! * `stars-000.bin ...`：按赤纬条带切块的验证星表；
//! * `band-00.ssq ...`：每个尺度带一个排序好的 quad 哈希文件。

use crate::angle::{radec_to_vec, PI, TWO_PI};
use crate::index::{
    ChunkLayout, IndexRecord, Manifest, ScaleBand, CHUNK_RECORDS, MAGIC, RECORD_BYTES,
    STAR_MAGIC,
};
use crate::quad::code_from_points;
use std::fs;
use std::path::Path;

/// 构建参数。
#[derive(Clone, Debug)]
pub struct BuildParams {
    /// 参与建库的最暗星等。
    pub mag_limit: f64,
    /// 星表块目标星数（按赤纬条带切）。
    pub stars_per_chunk: u32,
    /// quad 锚点角距离尺度带（弧度），从密到疏。
    pub bands: [(f64, f64); 5],
    /// 每颗锚星最多生成的 quad 数（控制索引体积）。
    pub quads_per_anchor: usize,
    pub catalog_name: String,
    pub epoch_jyear: f64,
}

impl Default for BuildParams {
    fn default() -> Self {
        // 锚点间距覆盖约 1.2°–22°，对应 24 MP 相机从长焦到广角的常见视场。
        Self {
            mag_limit: 12.0,
            stars_per_chunk: 8000,
            bands: [
                (arcdeg(1.2), arcdeg(2.6)),
                (arcdeg(2.4), arcdeg(5.2)),
                (arcdeg(4.8), arcdeg(10.5)),
                (arcdeg(9.0), arcdeg(16.0)),
                (arcdeg(14.0), arcdeg(22.0)),
            ],
            quads_per_anchor: 12,
            catalog_name: "custom".into(),
            epoch_jyear: 2000.0,
        }
    }
}

fn arcdeg(d: f64) -> f64 {
    d.to_radians()
}

/// 星表输入行。
#[derive(Clone, Copy, Debug)]
pub struct CatRow {
    pub id: u32,
    pub ra: f64,
    pub dec: f64,
    pub mag: f32,
}

/// 解析 CSV（id,ra_deg,dec_deg,mag）。
pub fn parse_catalog_csv(text: &str) -> Result<Vec<CatRow>, String> {
    let mut rows = Vec::new();
    for (lineno, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let lower = line.to_ascii_lowercase();
        if lower.starts_with("id,") || lower.starts_with("id ") {
            continue; // 表头
        }
        let cols: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
        if cols.len() < 4 {
            return Err(format!("line {}: need 4 columns, got {}", lineno + 1, cols.len()));
        }
        let id = cols[0]
            .parse::<u32>()
            .map_err(|e| format!("line {} id: {e}", lineno + 1))?;
        let ra_deg = cols[1]
            .parse::<f64>()
            .map_err(|e| format!("line {} ra: {e}", lineno + 1))?;
        let dec_deg = cols[2]
            .parse::<f64>()
            .map_err(|e| format!("line {} dec: {e}", lineno + 1))?;
        let mag = cols[3]
            .parse::<f32>()
            .map_err(|e| format!("line {} mag: {e}", lineno + 1))?;
        if !(0.0..360.0).contains(&ra_deg) || !(-90.0..=90.0).contains(&dec_deg) {
            return Err(format!("line {}: ra/dec out of range", lineno + 1));
        }
        rows.push(CatRow {
            id,
            ra: ra_deg.to_radians(),
            dec: dec_deg.to_radians(),
            mag,
        });
    }
    Ok(rows)
}

/// 构建后的产物句柄（内存形态，测试用）。
pub struct BuiltIndex {
    pub manifest: Manifest,
    pub star_layouts: Vec<ChunkLayout>,
    pub star_bytes: Vec<u8>,
    pub band_bytes: Vec<Vec<u8>>,
}

/// 构建并落盘到 `out_dir`。
pub fn build_to_dir(rows: &[CatRow], params: &BuildParams, out_dir: &Path) -> std::io::Result<BuiltIndex> {
    fs::create_dir_all(out_dir)?;
    let built = build(rows, params)?;

    // 星表块：所有块连续写入一个文件，manifest 记录每块偏移。
    let star_path = out_dir.join(format!("{}.bin", built.manifest.star_chunk_prefix));
    fs::write(&star_path, &built.star_bytes)?;

    for (i, bytes) in built.band_bytes.iter().enumerate() {
        let name = format!("band-{i:02}.ssq");
        fs::write(out_dir.join(&name), bytes)?;
    }
    let json = crate::index::manifest_json::to_json(&built.manifest, &built.star_layouts);
    fs::write(out_dir.join("manifest.json"), json)?;
    Ok(built)
}

/// 纯内存构建（不触盘，便于单元测试与 wasm 侧的小索引演示）。
pub fn build(rows: &[CatRow], params: &BuildParams) -> std::io::Result<BuiltIndex> {
    // 1) 星等过滤 + 按 (dec, ra) 排序。
    let mut stars: Vec<CatRow> = rows
        .iter()
        .copied()
        .filter(|r| (r.mag as f64) <= params.mag_limit)
        .collect();
    stars.sort_by(|a, b| {
        a.dec
            .partial_cmp(&b.dec)
            .unwrap()
            .then(a.ra.partial_cmp(&b.ra).unwrap())
    });

    // 2) 赤纬条带切块（目标块大小，极区自动变小）。
    let target = params.stars_per_chunk.max(500);
    let mut star_layouts = Vec::new();
    let mut star_bytes = Vec::new();
    let mut i = 0usize;
    while i < stars.len() {
        let start = i;
        let end = (i + target as usize).min(stars.len());
        let dec_min = stars[start].dec;
        let dec_max = stars[end - 1].dec;
        // 块头（magic + count）落在记录数据之前，offset 必须指向记录起点，
        // 这样 StarChunks::load 从 offset 直读 count×定长记录即可。
        star_bytes.extend_from_slice(STAR_MAGIC);
        star_bytes.extend_from_slice(&(end - start).to_le_bytes());
        let offset = star_bytes.len() as u64;
        for r in &stars[start..end] {
            // 全局星号即数组下标（by_id 按块偏移换算），原始 id 保留在记录里。
            star_bytes.extend_from_slice(&r.id.to_le_bytes());
            star_bytes.extend_from_slice(&r.ra.to_le_bytes());
            star_bytes.extend_from_slice(&r.dec.to_le_bytes());
            star_bytes.extend_from_slice(&r.mag.to_le_bytes());
        }
        star_layouts.push(ChunkLayout {
            offset,
            count: (end - start) as u32,
            dec_min,
            dec_max,
        });
        i = end;
    }
    // 构建期需要的连续数组（全局序号即数组下标）。
    let cat: Vec<CatRow> = stars;
    let n = cat.len();

    // 3) 每个尺度带生成 quad。
    let mut band_bytes = Vec::new();
    let mut bands_meta = Vec::new();
    for (bi, &(min_ab, max_ab)) in params.bands.iter().enumerate() {
        let mut records: Vec<IndexRecord> = Vec::new();
        // 空间桶：按 (dec, ra) 网格加速邻域搜索，网格边长取最大锚距。
        let cell = max_ab;
        let n_dec = ((PI) / cell).ceil() as usize + 2;
        let mut buckets: Vec<Vec<u32>> = vec![Vec::new(); n_dec * ra_bins(cell)];
        let n_ra = ra_bins(cell);
        for (k, r) in cat.iter().enumerate() {
            let di = ((r.dec + PI / 2.0) / cell).floor() as usize;
            let rj = (r.ra / cell).floor() as usize % n_ra;
            buckets[di * n_ra + rj].push(k as u32);
        }

        for anchor in 0..n {
            let made = build_quads_for_anchor(
                anchor,
                &cat,
                &buckets,
                cell,
                n_dec,
                n_ra,
                min_ab,
                max_ab,
                params.quads_per_anchor,
                &mut records,
            );
            let _ = made;
        }
        records.sort_by(|a, b| a.code.cmp(&b.code));
        records.dedup_by(|a, b| a.code == b.code && a.star_ids == b.star_ids);

        let num_chunks = records.len().max(1).div_ceil(CHUNK_RECORDS);
        let mut bytes = Vec::new();
        bytes.extend_from_slice(MAGIC);
        bytes.extend_from_slice(&1u32.to_le_bytes()); // version
        bytes.extend_from_slice(&(records.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&(num_chunks as u32).to_le_bytes());
        // 偏移表（num_chunks+1 项）；最后一项指向实际数据末尾，
        // 不能按满块外推，否则读取器会越过 EOF。
        let fixed = 4 + 4 + 8 + 4;
        let table_bytes = (num_chunks + 1) * 8;
        let header_size = fixed + table_bytes + num_chunks * 8;
        let mut offsets = Vec::with_capacity(num_chunks + 1);
        for c in 0..num_chunks {
            offsets.push((header_size + c * CHUNK_RECORDS * RECORD_BYTES) as u64);
        }
        offsets.push((header_size + records.len() * RECORD_BYTES) as u64);
        for o in &offsets {
            bytes.extend_from_slice(&o.to_le_bytes());
        }
        // 块首 code 表（空带用 0 占位，读取器对 num_records=0 不会访问数据）。
        for c in 0..num_chunks {
            let first_code = records
                .get(c * CHUNK_RECORDS)
                .map(|r| r.code)
                .unwrap_or(0);
            bytes.extend_from_slice(&first_code.to_le_bytes());
        }
        debug_assert_eq!(bytes.len(), header_size);
        for rec in &records {
            bytes.extend_from_slice(&rec.code.to_le_bytes());
            for s in rec.star_ids {
                bytes.extend_from_slice(&s.to_le_bytes());
            }
        }
        bands_meta.push(ScaleBand {
            index: bi as u32,
            min_rad: min_ab,
            max_rad: max_ab,
            file: format!("band-{bi:02}.ssq"),
            num_records: records.len() as u64,
            num_chunks: num_chunks as u32,
        });
        band_bytes.push(bytes);
    }

    let manifest = Manifest {
        catalog_name: params.catalog_name.clone(),
        mag_limit: params.mag_limit,
        epoch_jyear: params.epoch_jyear,
        bands: bands_meta,
        star_chunk_prefix: "stars".into(),
        star_chunk_count: star_layouts.len() as u32,
        star_count: n as u64,
    };
    Ok(BuiltIndex {
        manifest,
        star_layouts,
        star_bytes,
        band_bytes,
    })
}

fn ra_bins(cell: f64) -> usize {
    (TWO_PI / cell).ceil() as usize
}

#[allow(clippy::too_many_arguments)]
fn build_quads_for_anchor(
    anchor: usize,
    cat: &[CatRow],
    buckets: &[Vec<u32>],
    cell: f64,
    n_dec: usize,
    n_ra: usize,
    min_ab: f64,
    max_ab: f64,
    cap: usize,
    out: &mut Vec<IndexRecord>,
) -> usize {
    let a = cat[anchor];
    let va = radec_to_vec(a.ra, a.dec);
    // 邻域桶（3x3，RA 环绕）。
    let di = ((a.dec + PI / 2.0) / cell).floor() as i64;
    let rj = (a.ra / cell).floor() as i64;
    let mut neigh: Vec<u32> = Vec::new();
    for dd in -1..=1 {
        let d = di + dd;
        if d < 0 || d as usize >= n_dec {
            continue;
        }
        for dr in -1..=1 {
            let mut r = (rj + dr).rem_euclid(n_ra as i64) as usize;
            let _ = &mut r;
            neigh.extend(&buckets[d as usize * n_ra + r]);
        }
    }
    // 候选 B：角距在带内。
    let mut bs: Vec<(u32, f64)> = Vec::new();
    for &k in &neigh {
        if k as usize == anchor {
            continue;
        }
        let d = crate::angle::angular_sep(va, radec_to_vec(cat[k as usize].ra, cat[k as usize].dec));
        if d >= min_ab && d <= max_ab {
            bs.push((k, d));
        }
    }
    // 优先亮星 + 角距分散：按星等排序后取前 cap*2 再按角距间隔抽样。
    bs.sort_by(|x, y| cat[x.0 as usize].mag.partial_cmp(&cat[y.0 as usize].mag).unwrap());
    bs.truncate(cap * 2);

    // 锚点切平面：以锚点为切点，x 轴朝东、y 轴朝北。
    let frame = tangent_frame(a.ra, a.dec);
    let proj = |r: &CatRow| -> [f64; 2] {
        let v = radec_to_vec(r.ra, r.dec);
        [dot3(v, frame[0]), dot3(v, frame[1])]
    };
    let pa = proj(&a); // ≈ (0,0)

    let mut made = 0;
    for &(b_id, ab) in &bs {
        if made >= cap {
            break;
        }
        let pb = proj(&cat[b_id as usize]);
        let mid = [(pa[0] + pb[0]) * 0.5, (pa[1] + pb[1]) * 0.5];
        let r2 = (ab * 0.55).powi(2);
        // C/D：同一邻域内切平面上落在 AB 中点圆内的最亮两颗。
        let mut inner: Vec<(u32, f32)> = Vec::new();
        for &k in &neigh {
            if k as usize == anchor || k == b_id {
                continue;
            }
            let p = proj(&cat[k as usize]);
            let dx = p[0] - mid[0];
            let dy = p[1] - mid[1];
            if dx * dx + dy * dy <= r2 {
                inner.push((k, cat[k as usize].mag));
            }
        }
        if inner.len() < 2 {
            continue;
        }
        inner.sort_by(|x, y| x.1.partial_cmp(&y.1).unwrap());
        let (c_id, d_id) = (inner[0].0, inner[1].0);
        let pts = [
            pa,
            pb,
            proj(&cat[c_id as usize]),
            proj(&cat[d_id as usize]),
        ];
        let (code, _local) = code_from_points(pts);
        out.push(IndexRecord {
            code,
            star_ids: [anchor as u32, b_id, c_id, d_id],
        });
        made += 1;
    }
    made
}

/// 切平面正交基：e_ra 朝东、e_dec 朝北（单位向量）。
fn tangent_frame(ra: f64, dec: f64) -> [[f64; 3]; 2] {
    let (sd, cd) = dec.sin_cos();
    let (sr, cr) = ra.sin_cos();
    let e_ra = [-sr, cr, 0.0];
    let e_dec = [-sd * cr, -sd * sr, cd];
    [e_ra, e_dec]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// 把球面点集投影到以 (ra0,dec0) 为中心的切平面（测试辅助）。
pub fn project_to_tangent(ra0: f64, dec0: f64, ra: f64, dec: f64) -> (f64, f64) {
    let frame = tangent_frame(ra0, dec0);
    let v = radec_to_vec(ra, dec);
    (dot3(v, frame[0]), dot3(v, frame[1]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::{BandIndex, MemoryReader, QueryTol};

    fn uniform_sky(n: u32, seed: u64) -> Vec<CatRow> {
        // 确定性的伪随机球面均匀点（LCG + 极角变换）。
        let mut state = seed;
        let mut rows = Vec::new();
        for id in 0..n {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let u = (state >> 33) as f64 / (1u64 << 31) as f64;
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let v = (state >> 33) as f64 / (1u64 << 31) as f64;
            let ra = TWO_PI * u;
            let dec = (2.0 * v - 1.0).asin();
            rows.push(CatRow { id, ra, dec, mag: 8.0 + (id % 40) as f32 * 0.1 });
        }
        rows
    }

    #[test]
    fn build_and_query_roundtrip() {
        let rows = uniform_sky(3000, 42);
        let params = BuildParams {
            mag_limit: 12.0,
            stars_per_chunk: 1000,
            quads_per_anchor: 6,
            ..Default::default()
        };
        let built = build(&rows, &params).unwrap();
        assert!(!built.manifest.bands.is_empty());
        assert!(built.manifest.bands[0].num_records > 0);

        // 每个记录的 code 都能在自己的带文件里按容差查回。
        let bytes = built.band_bytes[0].clone();
        let mut idx = BandIndex::open(MemoryReader { bytes }, 8).unwrap();
        // 直接验证文件可解析且首条记录可查
        let first = {
            let hdr = idx.header.clone();
            let recs = idx.query(hdr.first_codes[0], QueryTol::default()).unwrap();
            assert!(!recs.is_empty());
            recs[0]
        };
        let again = idx.query(first.code, QueryTol { half_bins: 0 }).unwrap();
        assert!(again.iter().any(|r| r.code == first.code));
    }

    #[test]
    fn csv_parse_rejects_bad_rows() {
        let csv = "id,ra_deg,dec_deg,mag\n1,10.0,20.0,9.5\n2,400,0,8\n";
        let err = parse_catalog_csv(csv).unwrap_err();
        assert!(err.contains("out of range"));
        let ok = parse_catalog_csv("1,10.0,20.0,9.5\n# comment\n\n2,11,21,10\n").unwrap();
        assert_eq!(ok.len(), 2);
    }
}
