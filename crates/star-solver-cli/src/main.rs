//! star-solver 命令行：本地星表建库、照片盲解、合成自检。
//!
//! 子命令：
//! * `build <catalog.csv> <out-dir>`：从 CSV 星表构建分块哈希索引；
//! * `solve <index-dir> <image>`：对一张照片做离线板解，打印 WCS 与 FITS 头；
//! * `selftest`：用合成星空跑端到端流水线，无需任何外部数据。

mod catalog;

use star_solver_core::builder::{build_to_dir, parse_catalog_csv, BuildParams};
use star_solver_core::image::{detect_stars, estimate_background, DetectParams, GrayImage};
use star_solver_core::matcher::{solve, SolveConfig};
use std::path::Path;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        print_usage();
        std::process::exit(2);
    }
    let code = match args[1].as_str() {
        "build" => cmd_build(&args[2..]),
        "solve" => cmd_solve(&args[2..]),
        "selftest" => cmd_selftest(),
        "-h" | "--help" | "help" => {
            print_usage();
            0
        }
        other => {
            eprintln!("未知子命令：{other}");
            print_usage();
            2
        }
    };
    std::process::exit(code);
}

fn print_usage() {
    eprintln!(
        "star-solver {}\n\n\
         用法:\n  \
         star-solver build <catalog.csv> <index-dir>   构建分块哈希索引\n  \
         star-solver solve <index-dir> <image>         离线盲解一张照片\n  \
         star-solver selftest                           合成星空端到端自检\n",
        star_solver_core::VERSION
    );
}

fn load_gray(path: &Path) -> Result<GrayImage, String> {
    let dyn_img = image::open(path).map_err(|e| format!("无法解码图像 {path:?}: {e}"))?;
    // 大图按长边等比缩到 ≤4000 像元，控制检测耗时（亚像素精度仍足够）。
    let max_dim = 4000u32;
    let dyn_img = if dyn_img.width().max(dyn_img.height()) > max_dim {
        dyn_img.resize(
            max_dim,
            max_dim,
            image::imageops::FilterType::Triangle,
        )
    } else {
        dyn_img
    };
    let gray = image::imageops::grayscale(&dyn_img);
    let (width, height) = (gray.width() as usize, gray.height() as usize);
    Ok(GrayImage::new(width, height, gray.into_raw()))
}

fn cmd_build(args: &[String]) -> i32 {
    if args.len() != 2 {
        eprintln!("用法：star-solver build <catalog.csv> <index-dir>");
        return 2;
    }
    let csv_path = Path::new(&args[0]);
    let out_dir = Path::new(&args[1]);
    let text = match std::fs::read_to_string(csv_path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("读取星表失败：{e}");
            return 1;
        }
    };
    let rows = match parse_catalog_csv(&text) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("星表解析失败：{e}");
            return 1;
        }
    };
    eprintln!("读入 {} 颗星，开始构建索引……", rows.len());
    let t = Instant::now();
    let params = BuildParams::default();
    match build_to_dir(&rows, &params, out_dir) {
        Ok(built) => {
            eprintln!(
                "完成：{} 颗星、{} 个尺度带、{} 个星表块，用时 {:.1}s",
                built.manifest.star_count,
                built.manifest.bands.len(),
                built.manifest.star_chunk_count,
                t.elapsed().as_secs_f64()
            );
            0
        }
        Err(e) => {
            eprintln!("构建失败：{e}");
            1
        }
    }
}

fn cmd_solve(args: &[String]) -> i32 {
    if args.len() != 2 {
        eprintln!("用法：star-solver solve <index-dir> <image>");
        return 2;
    }
    let index_dir = Path::new(&args[0]);
    let img_path = Path::new(&args[1]);
    let mut cat = match catalog::FileCatalog::open(index_dir) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("打开索引失败：{e}");
            return 1;
        }
    };
    let img = match load_gray(img_path) {
        Ok(i) => i,
        Err(e) => {
            eprintln!("{e}");
            return 1;
        }
    };
    let bg = estimate_background(&img);
    eprintln!(
        "图像 {}x{}，背景 {:.1} ± {:.2}",
        img.width, img.height, bg.level, bg.sigma
    );
    let (stars, _) = detect_stars(&img, &DetectParams::default());
    eprintln!("检出 {} 颗星，开始盲解……", stars.len());

    let t = Instant::now();
    match solve(&mut cat, &stars, img.width, img.height, &SolveConfig::default()) {
        Ok(sol) => {
            println!("{}", "-".repeat(60));
            println!("center  RA = {:12.6} deg  DEC = {:11.6} deg", sol.center_ra_deg, sol.center_dec_deg);
            println!("scale   {:.3} x {:.3} arcsec/pixel", sol.pixel_scale_x, sol.pixel_scale_y);
            println!("roll    {:+.3} deg (CROTA2, north-to-east positive)", sol.roll_deg);
            println!("fov     {:.3} x {:.3} deg (diagonal {:.3} deg)", sol.fov_width_deg, sol.fov_height_deg, sol.fov_diag_deg);
            println!("inliers {} / {} stars, RMS {:.3} arcsec, max resid {:.3} arcsec", sol.num_inliers, sol.num_detected, sol.rms_arcsec, sol.max_resid_arcsec);
            println!("elapsed {:.2}s (scale band #{})", t.elapsed().as_secs_f64(), sol.band_index);
            println!("{}", "-".repeat(60));
            println!("{}", sol.wcs.to_fits_header());
            0
        }
        Err(e) => {
            eprintln!("解算失败：{e}");
            1
        }
    }
}

fn cmd_selftest() -> i32 {
    star_solver_core::selftest::run();
    0
}
