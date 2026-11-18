//! 基于本地文件的星表源，实现核心库的 [`CatalogSource`]。
//!
//! 验证星表块与各尺度带分别打开文件，按需读取并由索引层 LRU 缓存。

use star_solver_core::index::{
    manifest_json, BandIndex, ChunkLayout, FileBlockReader, IndexRecord, Manifest,
    QueryTol, StarChunks,
};
use star_solver_core::matcher::CatalogSource;
use star_solver_core::index::CatStar;
use std::path::{Path, PathBuf};

pub struct FileCatalog {
    root: PathBuf,
    manifest: Manifest,
    star_layouts: Vec<ChunkLayout>,
    bands: Vec<(f64, f64)>,
    band_indexes: Vec<Option<BandIndex<FileBlockReader>>>,
    stars: Option<StarChunks<FileBlockReader>>,
}

impl FileCatalog {
    pub fn open(root: &Path) -> std::io::Result<Self> {
        let manifest_text = std::fs::read_to_string(root.join("manifest.json"))?;
        let (manifest, star_layouts) = manifest_json::parse(&manifest_text)?;
        let bands = manifest
            .bands
            .iter()
            .map(|b| (b.min_rad, b.max_rad))
            .collect();
        let n_bands = manifest.bands.len();
        Ok(Self {
            root: root.to_path_buf(),
            manifest,
            star_layouts,
            bands,
            band_indexes: (0..n_bands).map(|_| None).collect(),
            stars: None,
        })
    }

    fn band(&mut self, band: u32) -> std::io::Result<&mut BandIndex<FileBlockReader>> {
        let idx = band as usize;
        if self.band_indexes[idx].is_none() {
            let fname = &self.manifest.bands[idx].file;
            let reader = FileBlockReader::open(&self.root.join(fname))?;
            self.band_indexes[idx] = Some(BandIndex::open(reader, 64)?);
        }
        Ok(self.band_indexes[idx].as_mut().unwrap())
    }

    fn star_chunks(&mut self) -> std::io::Result<&mut StarChunks<FileBlockReader>> {
        if self.stars.is_none() {
            let fname = format!("{}.bin", self.manifest.star_chunk_prefix);
            let reader = FileBlockReader::open(&self.root.join(&fname))?;
            self.stars = Some(StarChunks::new(reader, self.star_layouts.clone(), 32));
        }
        Ok(self.stars.as_mut().unwrap())
    }
}

impl CatalogSource for FileCatalog {
    fn query_band(
        &mut self,
        band: u32,
        code: u64,
        tol: QueryTol,
    ) -> std::io::Result<Vec<IndexRecord>> {
        self.band(band)?.query(code, tol)
    }

    fn stars_by_id(&mut self, ids: &[u32]) -> std::io::Result<Vec<CatStar>> {
        let mut out = Vec::with_capacity(ids.len());
        for &id in ids {
            out.push(self.star_chunks()?.by_id(id)?);
        }
        Ok(out)
    }

    fn stars_within(
        &mut self,
        ra: f64,
        dec: f64,
        radius: f64,
    ) -> std::io::Result<Vec<CatStar>> {
        self.star_chunks()?.within(ra, dec, radius)
    }

    fn bands(&self) -> &[(f64, f64)] {
        &self.bands
    }
}
