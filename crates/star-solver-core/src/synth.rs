//! 合成星空图：给定 WCS 与星表，渲染带高斯 PSF、背景和噪声的灰度图。
//!
//! 用途有两个：端到端流水线的确定性测试数据；无真实照片时的交互式 demo。
//! 噪声使用固定种子的 LCG，保证 wasm/原生两侧结果可复现。

use crate::image::GrayImage;
use crate::index::CatStar;
use crate::wcs::TanWcs;

/// 渲染参数。
#[derive(Clone, Debug)]
pub struct RenderParams {
    /// 背景计数（0..255）。
    pub background: f64,
    /// 背景高斯噪声 sigma（计数单位）。
    pub noise_sigma: f64,
    /// 星等零点：mag=0 的星在 PSF 中心像元的峰值幅度（背景以上计数）。
    pub mag_zero_peak: f64,
    /// PSF 高斯 sigma（像元）。
    pub psf_sigma: f64,
    /// 最暗渲染星等。
    pub faint_mag: f32,
    /// 噪声种子。
    pub seed: u64,
}

impl Default for RenderParams {
    fn default() -> Self {
        Self {
            background: 28.0,
            noise_sigma: 3.2,
            // mag=9 的星峰值约 15 计数（≈5σ），亮星峰值更高，保证全星表可检。
            mag_zero_peak: 60000.0,
            psf_sigma: 1.3,
            faint_mag: 11.5,
            seed: 0x57a1_5eed,
        }
    }
}

struct Lcg(u64);
impl Lcg {
    fn next_f64(&mut self) -> f64 {
        // MMIX LCG
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 33) as f64 + 0.5) / (1u64 << 31) as f64
    }
    /// Box-Muller 标准正态。
    fn normal(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-12);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

/// 渲染结果：图像 + 真值对应（像元坐标 ↔ 星表）。
pub struct RenderedSky {
    pub image: GrayImage,
    pub truth: Vec<RenderedStar>,
}

#[derive(Clone, Copy, Debug)]
pub struct RenderedStar {
    pub id: u32,
    pub x: f64,
    pub y: f64,
    pub mag: f32,
}

/// 渲染一张合成图。落在画幅外或过暗的星自动跳过。
pub fn render(
    wcs: &TanWcs,
    width: usize,
    height: usize,
    stars: &[CatStar],
    params: &RenderParams,
) -> RenderedSky {
    let mut f32buf = vec![params.background; width * height];
    let mut truth = Vec::new();
    let cutoff_r = (params.psf_sigma * 3.5).ceil() as i64;

    for s in stars {
        if s.mag > params.faint_mag {
            continue;
        }
        let (x, y) = wcs.world_to_pixel(s.ra, s.dec);
        if x < -cutoff_r as f64
            || y < -cutoff_r as f64
            || x >= width as f64 + cutoff_r as f64
            || y >= height as f64 + cutoff_r as f64
        {
            continue;
        }
        // mag_zero_peak 即 PSF 中心点幅度；高斯形状直接乘 exp(−r²/2σ²)。
        let peak = params.mag_zero_peak * 10f64.powf(-0.4 * s.mag as f64);
        let x0 = x.floor() as i64;
        let y0 = y.floor() as i64;
        for dy in -cutoff_r..=cutoff_r {
            let py = y0 + dy;
            if py < 0 || py >= height as i64 {
                continue;
            }
            for dx in -cutoff_r..=cutoff_r {
                let px = x0 + dx;
                if px < 0 || px >= width as i64 {
                    continue;
                }
                let d2 = (px as f64 - x).powi(2) + (py as f64 - y).powi(2);
                let psf = (-d2 / (2.0 * params.psf_sigma * params.psf_sigma)).exp();
                f32buf[py as usize * width + px as usize] += peak * psf;
            }
        }
        truth.push(RenderedStar {
            id: s.id,
            x,
            y,
            mag: s.mag,
        });
    }

    // 噪声 + 量化
    let mut rng = Lcg(params.seed);
    let data: Vec<u8> = f32buf
        .iter()
        .map(|v| {
            let n = v + rng.normal() * params.noise_sigma;
            n.round().clamp(0.0, 255.0) as u8
        })
        .collect();

    RenderedSky {
        image: GrayImage::new(width, height, data),
        truth,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::angle::{arcsec2rad, deg2rad};
    use crate::wcs::SipCoeffs;

    fn test_wcs() -> TanWcs {
        let scale = arcsec2rad(2.0);
        let roll = deg2rad(12.0);
        TanWcs {
            crpix: [300.0, 200.0],
            crval: [deg2rad(90.0), deg2rad(20.0)],
            cd: [
                [scale * roll.cos(), -scale * roll.sin()],
                [scale * roll.sin(), scale * roll.cos()],
            ],
            sip: SipCoeffs::default(),
        }
    }

    #[test]
    fn render_is_deterministic() {
        let wcs = test_wcs();
        let stars = vec![CatStar {
            id: 1,
            ra: deg2rad(90.05),
            dec: deg2rad(20.05),
            mag: 7.0,
        }];
        let p = RenderParams::default();
        let a = render(&wcs, 600, 400, &stars, &p);
        let b = render(&wcs, 600, 400, &stars, &p);
        assert_eq!(a.image.data, b.image.data);
        assert_eq!(a.truth.len(), 1);
    }
}
