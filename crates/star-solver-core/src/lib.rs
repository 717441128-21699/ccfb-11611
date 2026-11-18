//! star-solver-core —— 纯 std 的离线天文板解（plate solving）核心库。
//!
//! 流水线：
//! 1. [`image`]：背景估计 → 阈值分割 → 连通域 → 流量加权亚像素质心；
//! 2. [`quad`]：亮星构造几何不变四星组合，编码为 40 bit 哈希；
//! 3. [`index`] / [`builder`]：星表预计算的分块哈希索引（构建端 + 按需读取端）；
//! 4. [`matcher`]：多尺度假设、分块范围查询、RANSAC 式验证；
//! 5. [`wcs`]：TAN 投影 + 二阶 SIP 畸变拟合，输出 CRVAL/CD/CRPIX 与 FITS 风格头。
//!
//! 所有角度在库内部一律使用弧度，仅在序列化/输入输出边界换算成角度。

pub mod angle;
pub mod builder;
pub mod dso;
pub mod image;
pub mod index;
pub mod json;
pub mod lru;
pub mod matcher;
pub mod quad;
pub mod selftest;
pub mod synth;
pub mod wcs;

/// 库语义化版本（随 workspace 包版本同步）。
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
