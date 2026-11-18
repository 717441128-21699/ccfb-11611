//! 星表哈希索引的磁盘格式与按需读取。
//!
//! 一个索引目录（或打包后的单个 `.ssidx` 文件）由两部分组成：
//!
//! * `manifest.json`：尺度带、星等极限、分块布局等元信息；
//! * 每个尺度带一个 `band-XX.ssq` 文件：文件头 + 按 40 bit 编码排序的
//!   定长记录（u64 code + 4×u32 星表星号），记录按每 256 条切成 chunk，
//!   chunk 偏移表在文件头里，查询时只加载相交的 chunk。
//!
//! 验证用星表（matcher 拿到候选后要在全天空里找邻星做投票）按天区切块：
//! `stars-N.bin`，每块含若干颗星的 (id, ra, dec, mag)，按 (dec band, ra)
//! 组织，访问器按中心角距离做范围查询，只触及视场覆盖的少数几块。

use crate::lru::LruCache;
use crate::quad::CODE_BINS;
use std::io::{self, Read, Seek, SeekFrom};

/// 魔数与格式版本。
pub const MAGIC: &[u8; 4] = b"SSQ1";
pub const STAR_MAGIC: &[u8; 4] = b"SSC1";
/// 每个数据 chunk 包含的记录数。
pub const CHUNK_RECORDS: usize = 256;
/// 单条索引记录字节数：u64 code + 4×u32 星号。
pub const RECORD_BYTES: usize = 8 + 4 * 4;

/// 一个尺度带：quad 锚点角间距落在 [min_rad, max_rad] 内。
#[derive(Clone, Debug)]
pub struct ScaleBand {
    pub index: u32,
    pub min_rad: f64,
    pub max_rad: f64,
    pub file: String,
    /// 文件内记录总数。
    pub num_records: u64,
    /// chunk 数（= 偏移表长度 - 1）。
    pub num_chunks: u32,
}

/// 索引清单（manifest.json 的等价物，二进制头里也冗余一份关键值）。
#[derive(Clone, Debug)]
pub struct Manifest {
    pub catalog_name: String,
    pub mag_limit: f64,
    pub epoch_jyear: f64,
    pub bands: Vec<ScaleBand>,
    /// 验证星表块文件名前缀与块数。
    pub star_chunk_prefix: String,
    pub star_chunk_count: u32,
    pub star_count: u64,
}

/// 一条命中记录：几何编码 + 组成该 quad 的四颗星在验证星表里的全局星号。
#[derive(Clone, Copy, Debug)]
pub struct IndexRecord {
    pub code: u64,
    pub star_ids: [u32; 4],
}

/// 随机读抽象：本地文件、内存缓冲、浏览器 fetch+Range 适配层都实现它。
pub trait BlockReader {
    /// 读取从 `offset` 开始的 `len` 字节（实现侧可自行缓存）。
    fn read_at(&mut self, offset: u64, len: usize) -> io::Result<Vec<u8>>;
}

/// 本地文件实现。
pub struct FileBlockReader {
    file: std::fs::File,
}

impl FileBlockReader {
    pub fn open(path: &std::path::Path) -> io::Result<Self> {
        Ok(Self {
            file: std::fs::File::open(path)?,
        })
    }
}

impl BlockReader for FileBlockReader {
    fn read_at(&mut self, offset: u64, len: usize) -> io::Result<Vec<u8>> {
        self.file.seek(SeekFrom::Start(offset))?;
        let mut buf = vec![0u8; len];
        self.file.read_exact(&mut buf)?;
        Ok(buf)
    }
}

/// 全内存实现（测试/小索引）。
pub struct MemoryReader {
    pub bytes: Vec<u8>,
}

impl BlockReader for MemoryReader {
    fn read_at(&mut self, offset: u64, len: usize) -> io::Result<Vec<u8>> {
        let o = offset as usize;
        if o + len > self.bytes.len() {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "MemoryReader OOB"));
        }
        Ok(self.bytes[o..o + len].to_vec())
    }
}

/// 尺度带文件头（解析后的形态）。
#[derive(Clone, Debug)]
pub struct BandHeader {
    pub num_records: u64,
    pub num_chunks: u32,
    pub chunk_offsets: Vec<u64>,
    /// 每个 chunk 首条记录的 code，用于 O(log n) 定位。
    pub first_codes: Vec<u64>,
    pub header_bytes: u64,
}

/// 读取尺度带文件头。布局：
/// `MAGIC | u32 version | u64 num_records | u32 num_chunks`
/// `| [u64; num_chunks+1] chunk_offsets | [u64; num_chunks] first_codes`
pub fn read_band_header<R: BlockReader>(r: &mut R) -> io::Result<BandHeader> {
    let fixed = 4 + 4 + 8 + 4;
    let head = r.read_at(0, fixed)?;
    if &head[0..4] != MAGIC {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "bad band magic"));
    }
    let num_records = u64::from_le_bytes(head[8..16].try_into().unwrap());
    let num_chunks = u32::from_le_bytes(head[16..20].try_into().unwrap());
    let table_bytes = (num_chunks as usize + 1) * 8;
    let table = r.read_at(fixed as u64, table_bytes)?;
    let chunk_offsets = (0..=num_chunks as usize)
        .map(|i| u64::from_le_bytes(table[i * 8..i * 8 + 8].try_into().unwrap()))
        .collect();
    let code_table = r.read_at(fixed as u64 + table_bytes as u64, num_chunks as usize * 8)?;
    let first_codes = (0..num_chunks as usize)
        .map(|i| u64::from_le_bytes(code_table[i * 8..i * 8 + 8].try_into().unwrap()))
        .collect();
    Ok(BandHeader {
        num_records,
        num_chunks,
        chunk_offsets,
        first_codes,
        header_bytes: (fixed + table_bytes + num_chunks as usize * 8) as u64,
    })
}

fn parse_record(buf: &[u8]) -> IndexRecord {
    let code = u64::from_le_bytes(buf[0..8].try_into().unwrap());
    let mut star_ids = [0u32; 4];
    for i in 0..4 {
        star_ids[i] = u32::from_le_bytes(buf[8 + i * 4..12 + i * 4].try_into().unwrap());
    }
    IndexRecord { code, star_ids }
}

/// 4 维 box 查询的半宽（桶）。质心噪声 + 量化各贡献一点，±2 桶经验覆盖。
#[derive(Clone, Copy, Debug)]
pub struct QueryTol {
    pub half_bins: u32,
}

impl Default for QueryTol {
    fn default() -> Self {
        Self { half_bins: 2 }
    }
}

/// 尺度带读取器：持有文件头与 chunk LRU，提供编码范围查询。
pub struct BandIndex<R: BlockReader> {
    reader: R,
    /// 尺度带文件头（同模块测试需要读取首码表）。
    pub header: BandHeader,
    cache: LruCache<u32, Vec<IndexRecord>>,
    /// 统计：实际发生的底层读取次数（性能诊断用）。
    pub chunk_loads: u64,
}

impl<R: BlockReader> BandIndex<R> {
    pub fn open(mut reader: R, cache_chunks: usize) -> io::Result<Self> {
        let header = read_band_header(&mut reader)?;
        Ok(Self {
            reader,
            header,
            cache: LruCache::new(cache_chunks.max(1)),
            chunk_loads: 0,
        })
    }

    pub fn num_records(&self) -> u64 {
        self.header.num_records
    }

    fn load_chunk(&mut self, chunk: u32) -> io::Result<&Vec<IndexRecord>> {
        if !self.cache.contains(&chunk) {
            let start = self.header.chunk_offsets[chunk as usize];
            let end = self.header.chunk_offsets[chunk as usize + 1];
            let raw = self.reader.read_at(start, (end - start) as usize)?;
            let n = raw.len() / RECORD_BYTES;
            let mut recs = Vec::with_capacity(n);
            for i in 0..n {
                recs.push(parse_record(&raw[i * RECORD_BYTES..(i + 1) * RECORD_BYTES]));
            }
            self.chunk_loads += 1;
            self.cache.put(chunk, recs);
        }
        Ok(self.cache.get(&chunk).unwrap())
    }

    /// 查询编码落在 `center` 各维 ±tol 内的全部记录。
    ///
    /// 记录按 code 全局排序，box 的最小/最大编码给出扫描区间；
    /// 由于桶号分段打包，区间内会混入边界外记录，逐条用四维判据过滤。
    pub fn query(&mut self, center: u64, tol: QueryTol) -> io::Result<Vec<IndexRecord>> {
        let c = crate::quad::unpack(center);
        let h = tol.half_bins as i64;
        let lo = [
            c[0] as i64 - h,
            c[1] as i64 - h,
            c[2] as i64 - h,
            c[3] as i64 - h,
        ];
        let hi = [
            c[0] as i64 + h,
            c[1] as i64 + h,
            c[2] as i64 + h,
            c[3] as i64 + h,
        ];
        let clamp_b = |v: i64| v.clamp(0, CODE_BINS as i64 - 1) as u32;
        let code_lo = pack4(lo.map(clamp_b));
        let code_hi = pack4(hi.map(clamp_b));

        // 二分定位首条记录所在 chunk。
        let first_chunk = self.locate_chunk(code_lo);
        let mut out = Vec::new();
        let mut chunk = first_chunk;
        loop {
            let recs = self.load_chunk(chunk)?;
            for rec in recs {
                if rec.code > code_hi {
                    return Ok(out);
                }
                if rec.code < code_lo {
                    continue;
                }
                let b = crate::quad::unpack(rec.code);
                if (0..4).all(|i| (b[i] as i64 - c[i] as i64).abs() <= h) {
                    out.push(*rec);
                }
            }
            chunk += 1;
            if chunk >= self.header.num_chunks {
                break;
            }
        }
        Ok(out)
    }

    /// 找第一条记录 code >= target 所在的 chunk 号（按 first_codes 二分）。
    fn locate_chunk(&self, target: u64) -> u32 {
        let codes = &self.header.first_codes;
        match codes.binary_search_by(|c| c.cmp(&target)) {
            Ok(i) => i as u32,
            Err(0) => 0,
            Err(i) => (i - 1) as u32,
        }
    }
}

#[inline]
fn pack4(b: [u32; 4]) -> u64 {
    ((b[0] as u64) << 30) | ((b[1] as u64) << 20) | ((b[2] as u64) << 10) | b[3] as u64
}

/// 验证星表中的一颗星（角度均为弧度）。
#[derive(Clone, Copy, Debug)]
pub struct CatStar {
    pub id: u32,
    pub ra: f64,
    pub dec: f64,
    pub mag: f32,
}

/// 星表块头：MAGIC | u32 count，随后是定长 24 字节记录
/// （u32 id, f64 ra, f64 dec, f32 mag）。
pub const STAR_RECORD_BYTES: usize = 4 + 8 + 8 + 4;

/// 星表块访问器：块按 dec 条带顺序排列，块内按 ra 排序。
pub struct StarChunks<R: BlockReader> {
    reader: R,
    /// 每块的字节偏移与条数。
    pub layouts: Vec<ChunkLayout>,
    cache: LruCache<u32, Vec<CatStar>>,
    pub chunk_loads: u64,
}

#[derive(Clone, Copy, Debug)]
pub struct ChunkLayout {
    pub offset: u64,
    pub count: u32,
    pub dec_min: f64,
    pub dec_max: f64,
}

impl<R: BlockReader> StarChunks<R> {
    /// `layouts` 由 builder 产出的清单给出。
    pub fn new(reader: R, layouts: Vec<ChunkLayout>, cache_chunks: usize) -> Self {
        Self {
            reader,
            layouts,
            cache: LruCache::new(cache_chunks.max(1)),
            chunk_loads: 0,
        }
    }

    fn load(&mut self, idx: u32) -> io::Result<&Vec<CatStar>> {
        if !self.cache.contains(&idx) {
            let lay = self.layouts[idx as usize];
            let raw = self
                .reader
                .read_at(lay.offset, lay.count as usize * STAR_RECORD_BYTES)?;
            let mut v = Vec::with_capacity(lay.count as usize);
            for i in 0..lay.count as usize {
                let o = i * STAR_RECORD_BYTES;
                let id = u32::from_le_bytes(raw[o..o + 4].try_into().unwrap());
                let ra = f64::from_le_bytes(raw[o + 4..o + 12].try_into().unwrap());
                let dec = f64::from_le_bytes(raw[o + 12..o + 20].try_into().unwrap());
                let mag = f32::from_le_bytes(raw[o + 20..o + 24].try_into().unwrap());
                v.push(CatStar { id, ra, dec, mag });
            }
            self.chunk_loads += 1;
            self.cache.put(idx, v);
        }
        Ok(self.cache.get(&idx).unwrap())
    }

    /// 取中心 (ra,dec)、角半径 `radius` 内的全部星。
    /// 跨 RA=0 的查询由对每个块做角距离判据自然处理。
    pub fn within(&mut self, ra: f64, dec: f64, radius: f64) -> io::Result<Vec<CatStar>> {
        let (dmin, dmax) = (dec - radius, dec + radius);
        let cos_r = radius.cos();
        // 先收集要访问的块号（只用不可变数据），再逐块可变加载。
        let wanted: Vec<u32> = self
            .layouts
            .iter()
            .enumerate()
            .filter(|(_, lay)| lay.dec_max >= dmin && lay.dec_min <= dmax)
            .map(|(i, _)| i as u32)
            .collect();
        let mut out = Vec::new();
        for i in wanted {
            let stars = self.load(i)?;
            for s in stars {
                let v = crate::angle::radec_to_vec(s.ra, s.dec);
                let c = crate::angle::radec_to_vec(ra, dec);
                if crate::angle::dot(v, c) >= cos_r {
                    out.push(*s);
                }
            }
        }
        Ok(out)
    }

    /// 按全局星号取星（builder 保证星号 = 全体星按块拼接后的序号）。
    pub fn by_id(&mut self, id: u32) -> io::Result<CatStar> {
        // 先定位块与块内序号，避免在 layouts 的不可变迭代中可变借用 self。
        let mut acc = 0u32;
        let mut hit: Option<(u32, u32)> = None;
        for (i, lay) in self.layouts.iter().enumerate() {
            if id < acc + lay.count {
                hit = Some((i as u32, id - acc));
                break;
            }
            acc += lay.count;
        }
        let Some((chunk, local)) = hit else {
            return Err(io::Error::new(io::ErrorKind::NotFound, "star id out of range"));
        };
        let stars = self.load(chunk)?;
        Ok(stars[local as usize])
    }
}

/// 清单的极简 JSON 解析/序列化（避免核心库引入 serde）。
/// 仅支持本项目 manifest 的扁平结构 + 嵌套数组。
pub mod manifest_json {
    use super::{ChunkLayout, Manifest, ScaleBand};
    use crate::json::json_escape;
    use std::collections::HashMap;
    use std::io::{self, ErrorKind};

    pub fn to_json(m: &Manifest, chunks: &[ChunkLayout]) -> String {
        let mut s = String::new();
        s.push_str("{\n");
        s.push_str(&format!("  \"catalog_name\": \"{}\",\n", json_escape(&m.catalog_name)));
        s.push_str(&format!("  \"mag_limit\": {},\n", m.mag_limit));
        s.push_str(&format!("  \"epoch_jyear\": {},\n", m.epoch_jyear));
        s.push_str(&format!("  \"star_chunk_count\": {},\n", m.star_chunk_count));
        s.push_str(&format!("  \"star_count\": {},\n", m.star_count));
        s.push_str(&format!("  \"star_chunk_prefix\": \"{}\",\n", json_escape(&m.star_chunk_prefix)));
        s.push_str("  \"bands\": [\n");
        for (i, b) in m.bands.iter().enumerate() {
            s.push_str(&format!(
                "    {{\"index\":{},\"min_rad\":{:.12e},\"max_rad\":{:.12e},\"file\":\"{}\",\"num_records\":{},\"num_chunks\":{}}}{}\n",
                b.index, b.min_rad, b.max_rad, json_escape(&b.file), b.num_records, b.num_chunks,
                if i + 1 == m.bands.len() { "" } else { "," }
            ));
        }
        s.push_str("  ],\n  \"star_chunks\": [\n");
        for (i, c) in chunks.iter().enumerate() {
            s.push_str(&format!(
                "    {{\"offset\":{},\"count\":{},\"dec_min\":{:.12e},\"dec_max\":{:.12e}}}{}\n",
                c.offset, c.count, c.dec_min, c.dec_max,
                if i + 1 == chunks.len() { "" } else { "," }
            ));
        }
        s.push_str("  ]\n}\n");
        s
    }

    pub fn parse(text: &str) -> io::Result<(Manifest, Vec<ChunkLayout>)> {
        let v = crate::json::parse(text)
            .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?;
        let obj = v.as_object().ok_or_else(|| bad("manifest root"))?;
        let get = |k: &str| obj.get(k).ok_or_else(|| bad(k));
        let bands_v = get("bands")?.as_array().ok_or_else(|| bad("bands"))?;
        let mut bands = Vec::new();
        for bv in bands_v {
            let o = bv.as_object().ok_or_else(|| bad("band"))?;
            bands.push(ScaleBand {
                index: o["index"].as_u64().ok_or_else(|| bad("index"))? as u32,
                min_rad: o["min_rad"].as_f64().ok_or_else(|| bad("min_rad"))?,
                max_rad: o["max_rad"].as_f64().ok_or_else(|| bad("max_rad"))?,
                file: o["file"].as_str().ok_or_else(|| bad("file"))?.to_string(),
                num_records: o["num_records"].as_u64().ok_or_else(|| bad("num_records"))?,
                num_chunks: o["num_chunks"].as_u64().ok_or_else(|| bad("num_chunks"))? as u32,
            });
        }
        let chunks_v = get("star_chunks")?.as_array().ok_or_else(|| bad("star_chunks"))?;
        let mut layouts = Vec::new();
        for cv in chunks_v {
            let o = cv.as_object().ok_or_else(|| bad("chunk"))?;
            layouts.push(ChunkLayout {
                offset: o["offset"].as_u64().ok_or_else(|| bad("offset"))?,
                count: o["count"].as_u64().ok_or_else(|| bad("count"))? as u32,
                dec_min: o["dec_min"].as_f64().ok_or_else(|| bad("dec_min"))?,
                dec_max: o["dec_max"].as_f64().ok_or_else(|| bad("dec_max"))?,
            });
        }
        let m = Manifest {
            catalog_name: get("catalog_name")?.as_str().unwrap_or("unknown").to_string(),
            mag_limit: get("mag_limit")?.as_f64().unwrap_or(12.0),
            epoch_jyear: get("epoch_jyear")?.as_f64().unwrap_or(2000.0),
            bands,
            star_chunk_prefix: get("star_chunk_prefix")?
                .as_str()
                .unwrap_or("stars")
                .to_string(),
            star_chunk_count: get("star_chunk_count")?.as_u64().unwrap_or(0) as u32,
            star_count: get("star_count")?.as_u64().unwrap_or(0),
        };
        Ok((m, layouts))
    }

    fn bad(what: &str) -> io::Error {
        io::Error::new(ErrorKind::InvalidData, format!("manifest: missing {what}"))
    }

    /// 运行时把文件名映射到已打开的 BlockReader 的小工具。
    pub type ReaderMap = HashMap<String, Vec<u8>>;
}
