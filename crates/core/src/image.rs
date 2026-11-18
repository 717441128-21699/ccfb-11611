//! Grayscale image container, PGM I/O (native tests) and source extraction.
//!
//! The browser path decodes JPEG/PNG to grayscale in TypeScript and passes a
//! raw `f32`/`u16` buffer in, so this module intentionally has no image codec
//! dependencies. PGM (binary P5) is provided for native round-trip tests.

use crate::stats::median;

#[derive(Clone, Debug)]
pub struct GrayImage {
    pub width: usize,
    pub height: usize,
    /// Row-major `f32` pixel values in arbitrary flux units (higher = brighter).
    pub buf: Vec<f32>,
}

impl GrayImage {
    pub fn new(width: usize, height: usize) -> Self {
        GrayImage {
            width,
            height,
            buf: vec![0.0; width * height],
        }
    }

    #[inline]
    pub fn get(&self, x: usize, y: usize) -> f32 {
        self.buf[y * self.width + x]
    }

    #[inline]
    pub fn set(&mut self, x: usize, y: usize, v: f32) {
        self.buf[y * self.width + x] = v;
    }

    /// Write 16-bit binary PGM (P5). Values are offset/rescaled to `[0, 65535]`.
    pub fn write_pgm(&self, path: &str) -> std::io::Result<()> {
        use std::io::Write;
        let mut f = std::fs::File::create(path)?;
        let (mut lo, mut hi) = (f32::INFINITY, f32::NEG_INFINITY);
        for &v in &self.buf {
            if v < lo {
                lo = v;
            }
            if v > hi {
                hi = v;
            }
        }
        let span = (hi - lo).max(1e-9);
        write!(f, "P5\n{} {}\n65535\n", self.width, self.height)?;
        let mut out = Vec::with_capacity(self.buf.len() * 2);
        for &v in &self.buf {
            let q = ((v - lo) / span * 65535.0).round().clamp(0.0, 65535.0) as u16;
            out.extend_from_slice(&q.to_be_bytes());
        }
        f.write_all(&out)?;
        Ok(())
    }

    /// Read 8/16-bit binary PGM (P5). Comments starting with `#` are skipped.
    pub fn read_pgm(path: &str) -> std::io::Result<Self> {
        let bytes = std::fs::read(path)?;
        let mut pos = 0usize;

        fn read_token(b: &[u8], pos: &mut usize) -> std::io::Result<Vec<u8>> {
            loop {
                while *pos < b.len() && b[*pos].is_ascii_whitespace() {
                    *pos += 1;
                }
                if *pos < b.len() && b[*pos] == b'#' {
                    while *pos < b.len() && b[*pos] != b'\n' {
                        *pos += 1;
                    }
                } else {
                    break;
                }
            }
            let start = *pos;
            while *pos < b.len() && !b[*pos].is_ascii_whitespace() && b[*pos] != b'#' {
                *pos += 1;
            }
            if start == *pos {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "malformed PGM header",
                ));
            }
            Ok(b[start..*pos].to_vec())
        }

        let magic = read_token(&bytes, &mut pos)?;
        if magic != b"P5" {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "not a binary PGM (P5)",
            ));
        }
        let width: usize = String::from_utf8_lossy(&read_token(&bytes, &mut pos)?)
            .parse()
            .map_err(|_| std::io::Error::new(std::io::ErrorKind::InvalidData, "bad width"))?;
        let height: usize = String::from_utf8_lossy(&read_token(&bytes, &mut pos)?)
            .parse()
            .map_err(|_| std::io::Error::new(std::io::ErrorKind::InvalidData, "bad height"))?;
        let maxval: u32 = String::from_utf8_lossy(&read_token(&bytes, &mut pos)?)
            .parse()
            .map_err(|_| std::io::Error::new(std::io::ErrorKind::InvalidData, "bad maxval"))?;
        // Exactly one whitespace byte separates header from raster.
        pos += 1;
        let n = width * height;
        let mut buf = vec![0f32; n];
        if maxval < 256 {
            if bytes.len() < pos + n {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    "truncated PGM",
                ));
            }
            for i in 0..n {
                buf[i] = bytes[pos + i] as f32;
            }
        } else {
            if bytes.len() < pos + 2 * n {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    "truncated PGM",
                ));
            }
            for i in 0..n {
                buf[i] = u16::from_be_bytes([bytes[pos + 2 * i], bytes[pos + 2 * i + 1]])
                    as f32;
            }
        }
        Ok(GrayImage { width, height, buf })
    }
}

/// A detected point source with sub-pixel photometry.
#[derive(Clone, Debug)]
pub struct Source {
    /// Flux-weighted centroid, 0-based pixel coordinates (pixel center = integer).
    pub x: f64,
    pub y: f64,
    /// Total background-subtracted flux in the detection aperture.
    pub flux: f64,
    /// Peak background-subtracted pixel value.
    pub peak: f64,
    /// Number of pixels in the connected region.
    pub npix: usize,
    /// Second central moments (pixels^2), useful for shape filtering.
    pub mxx: f64,
    pub myy: f64,
    pub mxy: f64,
}

/// Background statistics estimated from a coarse grid of masked tiles.
#[derive(Clone, Copy, Debug)]
pub struct Background {
    pub global_mode: f32,
    pub global_sigma: f32,
}

/// Estimate sky background and its noise via a tile-grid median-of-mediants.
///
/// Tiles containing bright sources are rejected (their median is pulled up),
/// the global background level is the median of tile medians and the noise
/// the median absolute deviation of tile medians scaled to a Gaussian sigma.
pub fn estimate_background(img: &GrayImage, tile: usize) -> Background {
    let tile = tile.max(16);
    let nx = img.width.div_ceil(tile);
    let ny = img.height.div_ceil(tile);
    let mut tile_med: Vec<f32> = Vec::with_capacity(nx * ny);
    for ty in 0..ny {
        for tx in 0..nx {
            let x0 = tx * tile;
            let y0 = ty * tile;
            let x1 = (x0 + tile).min(img.width);
            let y1 = (y0 + tile).min(img.height);
            let mut v: Vec<f32> = Vec::with_capacity((x1 - x0) * (y1 - y0));
            for y in y0..y1 {
                let row = &img.buf[y * img.width..y * img.width + img.width];
                v.extend_from_slice(&row[x0..x1]);
            }
            if v.len() >= 16 {
                tile_med.push(median(&v));
            }
        }
    }
    let mode = median(&tile_med);
    let dev: Vec<f32> = tile_med.iter().map(|v| (v - mode).abs()).collect();
    let mad = median(&dev);
    // 1.4826 converts MAD to Gaussian sigma; floor avoids zero-noise divides.
    let sigma = (1.4826 * mad).max(1e-6);
    Background {
        global_mode: mode,
        global_sigma: sigma,
    }
}

/// Source extraction configuration.
#[derive(Clone, Debug)]
pub struct ExtractConfig {
    /// Detection threshold in background-sigma units above the local level.
    pub threshold_sigma: f32,
    /// Reject regions smaller than this many pixels.
    pub min_pixels: usize,
    /// Reject regions larger than this many pixels (galaxies/blooming/edges).
    pub max_pixels: usize,
    /// Reject sources whose centroid is closer than this many pixels to a border.
    pub border: usize,
    /// Optional cap on returned sources (brightest kept), `None` = all.
    pub max_sources: Option<usize>,
}

impl Default for ExtractConfig {
    fn default() -> Self {
        ExtractConfig {
            threshold_sigma: 5.0,
            min_pixels: 2,
            max_pixels: 400,
            border: 4,
            max_sources: Some(2000),
        }
    }
}

/// Extract point sources from a grayscale image.
///
/// Two-pass approach:
/// 1. per-pixel mask `v > background + k·sigma` with a mild local background
///    correction (8x8 tile median) to handle gradients;
/// 2. 4-connected component labeling with run-length flood fill;
/// 3. flux-weighted centroid and second moments over background-subtracted
///    pixels inside each component.
pub fn extract_sources(img: &GrayImage, cfg: &ExtractConfig) -> Vec<Source> {
    let (w, h) = (img.width, img.height);
    let bg = estimate_background(img, 64);

    // Local background: tile medians on a coarse grid, bilinear-ish via lookup.
    let lt = 32usize;
    let gnx = w.div_ceil(lt);
    let gny = h.div_ceil(lt);
    let mut grid = vec![0f32; gnx * gny];
    for gy in 0..gny {
        for gx in 0..gnx {
            let x0 = gx * lt;
            let y0 = gy * lt;
            let x1 = (x0 + lt).min(w);
            let y1 = (y0 + lt).min(h);
            let mut v = Vec::with_capacity((x1 - x0) * (y1 - y0));
            for y in y0..y1 {
                for x in x0..x1 {
                    v.push(img.buf[y * w + x]);
                }
            }
            // Fall back to global mode for tiny/edge tiles.
            grid[gy * gnx + gx] = if v.len() >= 8 {
                median(&v)
            } else {
                bg.global_mode
            };
        }
    }
    #[inline]
    fn local_bg(grid: &[f32], gnx: usize, gny: usize, x: usize, y: usize, lt: usize) -> f32 {
        let gx = (x / lt).min(gnx - 1);
        let gy = (y / lt).min(gny - 1);
        grid[gy * gnx + gx]
    }

    let thr = bg.global_mode + cfg.threshold_sigma * bg.global_sigma;
    let _ = thr; // actual comparison uses local level below
    let mut mask = vec![false; w * h];
    for y in 0..h {
        let row_off = y * w;
        for x in 0..w {
            let v = img.buf[row_off + x];
            let level = local_bg(&grid, gnx, gny, x, y, lt);
            if v > level + cfg.threshold_sigma * bg.global_sigma {
                mask[row_off + x] = true;
            }
        }
    }

    // 4-connected components via iterative flood fill with a reusable stack.
    /// Accumulator for one connected component's flux-weighted moments.
    #[derive(Default)]
    struct CompAcc {
        sum: f64,
        sx: f64,
        sy: f64,
        sxx: f64,
        syy: f64,
        sxy: f64,
        peak: f64,
        npix: usize,
    }
    let mut labels = vec![0u32; w * h];
    let mut comps: Vec<CompAcc> = Vec::new();
    let mut stack: Vec<usize> = Vec::new();
    let mut next_label = 0u32;

    for seed in 0..w * h {
        if !mask[seed] || labels[seed] != 0 {
            continue;
        }
        next_label += 1;
        let label = next_label;
        labels[seed] = label;
        stack.clear();
        stack.push(seed);

        let mut sum = 0f64;
        let mut sx = 0f64;
        let mut sy = 0f64;
        let mut sxx = 0f64;
        let mut syy = 0f64;
        let mut sxy = 0f64;
        let mut peak = 0f64;
        let mut npix = 0usize;
        let mut minx = w;
        let mut maxx = 0usize;
        let mut miny = h;
        let mut maxy = 0usize;

        while let Some(idx) = stack.pop() {
            let x = idx % w;
            let y = idx / w;
            let level = local_bg(&grid, gnx, gny, x, y, lt) as f64;
            let f = (img.buf[idx] as f64 - level).max(0.0);
            let xf = x as f64;
            let yf = y as f64;
            sum += f;
            sx += f * xf;
            sy += f * yf;
            sxx += f * xf * xf;
            syy += f * yf * yf;
            sxy += f * xf * yf;
            if f > peak {
                peak = f;
            }
            npix += 1;
            if x < minx {
                minx = x;
            }
            if x > maxx {
                maxx = x;
            }
            if y < miny {
                miny = y;
            }
            if y > maxy {
                maxy = y;
            }
            // 4-neighbours.
            if x > 0 {
                let n = idx - 1;
                if mask[n] && labels[n] == 0 {
                    labels[n] = label;
                    stack.push(n);
                }
            }
            if x + 1 < w {
                let n = idx + 1;
                if mask[n] && labels[n] == 0 {
                    labels[n] = label;
                    stack.push(n);
                }
            }
            if y > 0 {
                let n = idx - w;
                if mask[n] && labels[n] == 0 {
                    labels[n] = label;
                    stack.push(n);
                }
            }
            if y + 1 < h {
                let n = idx + w;
                if mask[n] && labels[n] == 0 {
                    labels[n] = label;
                    stack.push(n);
                }
            }
        }

        if npix < cfg.min_pixels || npix > cfg.max_pixels {
            continue;
        }
        if minx < cfg.border
            || miny < cfg.border
            || maxx + cfg.border >= w
            || maxy + cfg.border >= h
        {
            continue;
        }
        if sum <= 0.0 {
            continue;
        }
        comps.push(CompAcc {
            sum,
            sx,
            sy,
            sxx,
            syy,
            sxy,
            peak,
            npix,
        });
    }

    let mut sources: Vec<Source> = comps
        .into_iter()
        .map(|c| {
            let cx = c.sx / c.sum;
            let cy = c.sy / c.sum;
            Source {
                x: cx,
                y: cy,
                flux: c.sum,
                peak: c.peak,
                npix: c.npix,
                mxx: (c.sxx / c.sum - cx * cx).max(0.0),
                myy: (c.syy / c.sum - cy * cy).max(0.0),
                mxy: c.sxy / c.sum - cx * cy,
            }
        })
        .collect();

    // Brightest first; quad generation only ever wants the top N anyway.
    sources.sort_by(|a, b| b.flux.partial_cmp(&a.flux).unwrap());
    if let Some(cap) = cfg.max_sources {
        sources.truncate(cap);
    }
    sources
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pgm_roundtrip() {
        let mut img = GrayImage::new(12, 8);
        img.set(3, 4, 1000.0);
        img.set(10, 1, 500.0);
        let dir = std::env::temp_dir().join("platersol_pgm_test.pgm");
        let p = dir.to_str().unwrap();
        img.write_pgm(p).unwrap();
        let img2 = GrayImage::read_pgm(p).unwrap();
        assert_eq!(img2.width, 12);
        assert_eq!(img2.height, 8);
        let _ = std::fs::remove_file(p);
    }

    #[test]
    fn finds_injected_sources() {
        let mut img = GrayImage::new(200, 200);
        // Uniform sky + read noise pattern.
        for v in img.buf.iter_mut() {
            *v = 100.0;
        }
        // Inject three PSF-like blobs.
        let spots = [(50.0f32, 60.0f32, 800.0), (150.0, 40.0, 600.0), (120.0, 150.0, 400.0)];
        for &(sx, sy, amp) in &spots {
            for dy in -3..=3 {
                for dx in -3..=3 {
                    let r2 = (dx * dx + dy * dy) as f32;
                    let x = (sx as i32 + dx) as usize;
                    let y = (sy as i32 + dy) as usize;
                    img.buf[y * 200 + x] += amp * (-r2 / 2.0).exp();
                }
            }
        }
        let srcs = extract_sources(
            &img,
            &ExtractConfig {
                threshold_sigma: 8.0,
                ..Default::default()
            },
        );
        assert_eq!(srcs.len(), 3, "should find exactly 3 sources, got {srcs:?}");
        for (s, (sx, sy, _)) in srcs.iter().zip(spots.iter()) {
            assert!((s.x - *sx as f64).abs() < 0.5, "x centroid {}", s.x);
            assert!((s.y - *sy as f64).abs() < 0.5, "y centroid {}", s.y);
        }
        // Flux ordering preserved.
        assert!(srcs[0].flux >= srcs[1].flux && srcs[1].flux >= srcs[2].flux);
    }
}
