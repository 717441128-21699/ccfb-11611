//! 天球投影与 WCS 拟合：TAN（切平面）+ 二阶 SIP 畸变。
//!
//! 模型（FITS WCS 标准的简化实现）：
//! * 像元偏移 (u,v) = (x - CRPIX1, y - CRPIX2)；
//! * SIP 前向修正 q = (u + A(u,v), v + B(u,v))，A/B 为二阶三项式；
//! * 中间切平面坐标 (ξ,η) = CD · q（单位弧度）；
//! * TAN 投影在 (CRVAL1, CRVAL2) 处反投影到天球。
//!
//! 拟合走“切平面投影 + 线性最小二乘 + 鲁棒剔除 + 一轮 SIP 修正”，
//! 全部用正规方程（4×4 / 6×6），无外部线性代数库。

use crate::angle::{rad2arcsec, radec_to_vec, vec_to_radec, PI, TWO_PI};

/// 拟合/投影用的对应点：像元坐标 ↔ 星表赤经赤纬（弧度）。
#[derive(Clone, Copy, Debug)]
pub struct Correspondence {
    pub px: f64,
    pub py: f64,
    pub ra: f64,
    pub dec: f64,
}

/// 二阶 SIP 系数（A 作用于 u 方向，B 作用于 v 方向）。
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SipCoeffs {
    pub a20: f64,
    pub a21: f64,
    pub a02: f64,
    pub b20: f64,
    pub b21: f64,
    pub b02: f64,
}

/// 一个完整的 TAN+SIP 解。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TanWcs {
    pub crpix: [f64; 2],
    pub crval: [f64; 2], // 弧度 (ra, dec)
    pub cd: [[f64; 2]; 2],
    pub sip: SipCoeffs,
}

impl TanWcs {
    /// 像元 → 赤经赤纬（弧度）。SIP 反演用 3 次牛顿迭代，亚毫角秒收敛。
    pub fn pixel_to_world(&self, x: f64, y: f64) -> (f64, f64) {
        let (u0, v0) = (x - self.crpix[0], y - self.crpix[1]);
        // 求 q 使 F(q) = q + sip(q) − (u0,v0) = 0。
        // 雅可比 J = I + ∂sip，2x2 直接求逆。
        let (mut qu, mut qv) = (u0, v0);
        for _ in 0..4 {
            let (du, dv) = self.sip_forward(qu, qv);
            let (da_du, da_dv, db_du, db_dv) = self.sip_jacobian(qu, qv);
            let fu = qu + du - u0;
            let fv = qv + dv - v0;
            let j00 = 1.0 + da_du;
            let j01 = da_dv;
            let j10 = db_du;
            let j11 = 1.0 + db_dv;
            let det = j00 * j11 - j01 * j10;
            if det.abs() < 1e-18 {
                break;
            }
            // dq = −J⁻¹·F
            let dqu = -(j11 * fu - j01 * fv) / det;
            let dqv = -(-j10 * fu + j00 * fv) / det;
            qu += dqu;
            qv += dqv;
            if dqu.hypot(dqv) < 1e-11 {
                break;
            }
        }
        let xi = self.cd[0][0] * qu + self.cd[0][1] * qv;
        let eta = self.cd[1][0] * qu + self.cd[1][1] * qv;
        tan_deproject(self.crval[0], self.crval[1], xi, eta)
    }

    /// 赤经赤纬（弧度）→ 像元坐标。
    pub fn world_to_pixel(&self, ra: f64, dec: f64) -> (f64, f64) {
        let (xi, eta) = tan_project(self.crval[0], self.crval[1], ra, dec);
        let det = self.cd[0][0] * self.cd[1][1] - self.cd[0][1] * self.cd[1][0];
        if det.abs() < 1e-18 {
            return (self.crpix[0], self.crpix[1]);
        }
        // q = CD^-1 (xi,eta)
        let qu = (self.cd[1][1] * xi - self.cd[0][1] * eta) / det;
        let qv = (-self.cd[1][0] * xi + self.cd[0][0] * eta) / det;
        let (du, dv) = self.sip_forward(qu, qv);
        (self.crpix[0] + qu + du, self.crpix[1] + qv + dv)
    }

    fn sip_forward(&self, u: f64, v: f64) -> (f64, f64) {
        let s = &self.sip;
        (
            s.a20 * u * u + s.a21 * u * v + s.a02 * v * v,
            s.b20 * u * u + s.b21 * u * v + s.b02 * v * v,
        )
    }

    fn sip_jacobian(&self, u: f64, v: f64) -> (f64, f64, f64, f64) {
        let s = &self.sip;
        (
            2.0 * s.a20 * u + s.a21 * v,
            s.a21 * u + 2.0 * s.a02 * v,
            2.0 * s.b20 * u + s.b21 * v,
            s.b21 * u + 2.0 * s.b02 * v,
        )
    }

    /// 像元比例尺（角秒/像元），取 x、y 方向。
    pub fn pixel_scales_arcsec(&self) -> (f64, f64) {
        let sx = rad2arcsec(self.cd[0][0].hypot(self.cd[1][0]));
        let sy = rad2arcsec(self.cd[0][1].hypot(self.cd[1][1]));
        (sx, sy)
    }

    /// 相机滚动角（度），采用 FITS CROTA2 约定：图像 +y 轴相对天体北方向，
    /// 正向为"北转向东"（从北天极看向切点时为顺时针）。
    ///
    /// 切平面基为 ξ 朝东、η 朝北，CD 第二列 `(cd12, cd22)` 是 +y 像元
    /// 步进的 (东, 北) 分量，故 CROTA2 = atan2(−cd12, cd22)。
    pub fn roll_deg(&self) -> f64 {
        (-self.cd[0][1]).atan2(self.cd[1][1]).to_degrees()
    }

    /// 图像四角的天球坐标与对角视场角（弧度）。
    pub fn field_of_view(&self, width: f64, height: f64) -> Fov {
        let corners = [
            self.pixel_to_world(0.0, 0.0),
            self.pixel_to_world(width, 0.0),
            self.pixel_to_world(width, height),
            self.pixel_to_world(0.0, height),
        ];
        let sep = |a: (f64, f64), b: (f64, f64)| {
            crate::angle::angular_sep(radec_to_vec(a.0, a.1), radec_to_vec(b.0, b.1))
        };
        let diag = sep(corners[0], corners[2]).max(sep(corners[1], corners[3]));
        let hfov = sep(corners[0], corners[1]).max(sep(corners[2], corners[3]));
        let vfov = sep(corners[0], corners[3]).max(sep(corners[1], corners[2]));
        Fov {
            corners,
            width_rad: hfov,
            height_rad: vfov,
            diagonal_rad: diag,
        }
    }

    /// 导出 FITS WCS 头卡片（角度单位度）。
    pub fn to_fits_header(&self) -> String {
        let card = |key: &str, val: &str, comment: &str| {
            // FITS 卡片固定 80 列
            let mut left = format!("{key:<8}= {val}");
            if left.len() > 30 {
                left.truncate(30);
            }
            let c = if comment.is_empty() {
                String::new()
            } else {
                format!(" / {comment}")
            };
            let mut line = format!("{left:<30}{c}");
            if line.len() > 80 {
                line.truncate(80);
            }
            format!("{line:<80}")
        };
        let d = |r: f64| format!("{:>20.10}", r.to_degrees());
        let mut h = String::new();
        h.push_str(&card("CRPIX1", &format!("{:>20.10}", self.crpix[0]), "[pix]"));
        h.push('\n');
        h.push_str(&card("CRPIX2", &format!("{:>20.10}", self.crpix[1]), "[pix]"));
        h.push('\n');
        h.push_str(&card("CRVAL1", &d(self.crval[0]), "[deg] RA of reference"));
        h.push('\n');
        h.push_str(&card("CRVAL2", &d(self.crval[1]), "[deg] DEC of reference"));
        h.push('\n');
        h.push_str(&card("CD1_1", &format!("{:>20.12e}", self.cd[0][0].to_degrees()), "[deg/pix]"));
        h.push('\n');
        h.push_str(&card("CD1_2", &format!("{:>20.12e}", self.cd[0][1].to_degrees()), "[deg/pix]"));
        h.push('\n');
        h.push_str(&card("CD2_1", &format!("{:>20.12e}", self.cd[1][0].to_degrees()), "[deg/pix]"));
        h.push('\n');
        h.push_str(&card("CD2_2", &format!("{:>20.12e}", self.cd[1][1].to_degrees()), "[deg/pix]"));
        h.push('\n');
        h.push_str(&card("CTYPE1", "'RA---TAN-SIP'    ", "Gnomonic + SIP"));
        h.push('\n');
        h.push_str(&card("CTYPE2", "'DEC--TAN-SIP'    ", "Gnomonic + SIP"));
        h.push('\n');
        h.push_str(&card("CUNIT1", "'deg'              ", ""));
        h.push('\n');
        h.push_str(&card("CUNIT2", "'deg'              ", ""));
        h.push('\n');
        h.push_str(&card("RADESYS", "'ICRS'             ", "J2000"));
        h.push('\n');
        h
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Fov {
    pub corners: [(f64, f64); 4],
    pub width_rad: f64,
    pub height_rad: f64,
    pub diagonal_rad: f64,
}

/// TAN 正向投影：天球 → 以 (ra0,dec0) 为切点的中间坐标 (ξ,η)。
pub fn tan_project(ra0: f64, dec0: f64, ra: f64, dec: f64) -> (f64, f64) {
    let (sdec0, cdec0) = dec0.sin_cos();
    let dra = wrap_dra(ra - ra0);
    let (sdec, cdec) = dec.sin_cos();
    let (sdra, cdra) = dra.sin_cos();
    let denom = sdec * sdec0 + cdec * cdec0 * cdra;
    let xi = cdec * sdra / denom;
    let eta = (sdec * cdec0 - cdec * sdec0 * cdra) / denom;
    (xi, eta)
}

/// TAN 反向投影：中间坐标 (ξ,η) → 天球。
///
/// 标准 gnomonic 反投影（与 [`tan_project`] 严格互逆）：
/// u = (n0 + ξ·e_east + η·e_north) / √(1 + ξ² + η²)，
/// 即把切平面向量归一化回单位球面。切平面基：ξ 朝东、η 朝北。
pub fn tan_deproject(ra0: f64, dec0: f64, xi: f64, eta: f64) -> (f64, f64) {
    let (sdec0, cdec0) = dec0.sin_cos();
    let (sra0, cra0) = ra0.sin_cos();
    let n0 = [cdec0 * cra0, cdec0 * sra0, sdec0];
    let e_east = [-sra0, cra0, 0.0];
    let e_north = [-sdec0 * cra0, -sdec0 * sra0, cdec0];
    let norm = (1.0 + xi * xi + eta * eta).sqrt();
    let u = [
        (n0[0] + xi * e_east[0] + eta * e_north[0]) / norm,
        (n0[1] + xi * e_east[1] + eta * e_north[1]) / norm,
        (n0[2] + xi * e_east[2] + eta * e_north[2]) / norm,
    ];
    vec_to_radec(u)
}

/// RA 差归一化到 (-π, π]。
fn wrap_dra(d: f64) -> f64 {
    let mut r = d % TWO_PI;
    if r <= -PI {
        r += TWO_PI;
    } else if r > PI {
        r -= TWO_PI;
    }
    r
}

/// 对称矩阵（全 n×n 存储）的高斯消元求解，返回 None 表示奇异。
fn solve_normal(mut m: Vec<f64>, mut b: Vec<f64>, n: usize) -> Option<Vec<f64>> {
    for col in 0..n {
        // 选主元
        let mut piv = col;
        for r in (col + 1)..n {
            if m[r * n + col].abs() > m[piv * n + col].abs() {
                piv = r;
            }
        }
        if m[piv * n + col].abs() < 1e-15 {
            return None;
        }
        if piv != col {
            for k in 0..n {
                m.swap(col * n + k, piv * n + k);
            }
            b.swap(col, piv);
        }
        for r in (col + 1)..n {
            let f = m[r * n + col] / m[col * n + col];
            for k in col..n {
                m[r * n + k] -= f * m[col * n + k];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for k in (i + 1)..n {
            s -= m[i * n + k] * x[k];
        }
        x[i] = s / m[i * n + i];
    }
    Some(x)
}

/// 正规方程累加（全矩阵）：仅 AᵀA。
fn accumulate_ata(ata: &mut [f64], row: &[f64], n: usize) {
    for i in 0..n {
        for j in 0..n {
            ata[i * n + j] += row[i] * row[j];
        }
    }
}

/// 正规方程累加：仅 Aᵀb（可对多个右端分别调用）。
fn accumulate_atb(atb: &mut [f64], row: &[f64], y: f64, n: usize) {
    for i in 0..n {
        atb[i] += row[i] * y;
    }
}

/// 拟合结果：WCS + 残差统计。
pub struct FitResult {
    pub wcs: TanWcs,
    pub rms_arcsec: f64,
    pub max_arcsec: f64,
    pub used: usize,
    /// 与输入 pairs 等长的内点掩码。
    pub inliers: Vec<bool>,
}

/// 从对应点拟合 TAN（不带 SIP），自动做 2.5σ 鲁棒剔除。
///
/// 双层迭代：
/// * 内层——固定一组内点，切点（CRVAL）从内点球面质心起步，每轮线性 LS
///   解出自由 CRPIX/CD 后，用当前 WCS 把图像中心指向作为新切点重拟合，
///   直到切点移动量 < 0.1 毫角秒（通常 2~3 轮），消除 TAN 非线性残差；
/// * 外层——按收敛解的角残差做 2.5σ MAD 裁剪，最多 3 轮。
pub fn fit_tan(pairs: &[Correspondence], image_w: f64, image_h: f64) -> Option<FitResult> {
    if pairs.len() < 4 {
        return None;
    }
    let center = [image_w * 0.5, image_h * 0.5];
    let mut active: Vec<bool> = vec![true; pairs.len()];
    let mut wcs = None;
    let mut stats = (0.0, 0.0);

    for _outer in 0..3 {
        let sel: Vec<&Correspondence> = pairs
            .iter()
            .enumerate()
            .filter(|(i, _)| active[*i])
            .map(|(_, p)| p)
            .collect();
        if sel.len() < 4 {
            return None;
        }

        // 切点初值：内点单位向量球面平均（对 RA=0 跨界稳健）。
        let mean = sel
            .iter()
            .map(|p| radec_to_vec(p.ra, p.dec))
            .fold([0.0; 3], |acc, v| [acc[0] + v[0], acc[1] + v[1], acc[2] + v[2]]);
        let (mut ra0, mut dec0) = vec_to_radec(crate::angle::normalize(mean));

        // 内层切点收敛。
        let mut candidate = None;
        for _inner in 0..8 {
            // 无约束线性模型 ξ = a0 + a1·x + a2·y，η 同理；CRPIX 自由。
            // AᵀA 对两个右端是同一个，每点只能累加一次。
            let n = 3;
            let (mut ata, mut atb_xi, mut atb_eta) =
                (vec![0.0; n * n], vec![0.0; n], vec![0.0; n]);
            for p in &sel {
                let (xi, eta) = tan_project(ra0, dec0, p.ra, p.dec);
                let row = [1.0, p.px, p.py];
                accumulate_ata(&mut ata, &row, n);
                accumulate_atb(&mut atb_xi, &row, xi, n);
                accumulate_atb(&mut atb_eta, &row, eta, n);
            }
            let a = solve_normal(ata.clone(), atb_xi, n)?;
            let b = solve_normal(ata, atb_eta, n)?;
            let cd = [[a[1], a[2]], [b[1], b[2]]];
            // (a0,b0) = −CD·CRPIX → CRPIX = −CD⁻¹·(a0,b0)
            let det = cd[0][0] * cd[1][1] - cd[0][1] * cd[1][0];
            if det.abs() < 1e-25 {
                return None;
            }
            let inv_det = 1.0 / det;
            let crpix = [
                -inv_det * (cd[1][1] * a[0] - cd[0][1] * b[0]),
                -inv_det * (-cd[1][0] * a[0] + cd[0][0] * b[0]),
            ];
            let current = TanWcs {
                crpix,
                crval: [ra0, dec0],
                cd,
                sip: SipCoeffs::default(),
            };
            // 把切点移到当前模型的图像中心指向，再重拟合。
            let (nra, ndec) = current.pixel_to_world(center[0], center[1]);
            let shift = crate::angle::angular_sep(
                radec_to_vec(nra, ndec),
                radec_to_vec(ra0, dec0),
            );
            ra0 = nra;
            dec0 = ndec;
            candidate = Some(current);
            if shift < crate::angle::arcsec2rad(0.0001) {
                break;
            }
        }
        let candidate = candidate?;

        // 残差（角距离）与 2.5σ MAD 裁剪。
        let mut rms_list: Vec<f64> = sel
            .iter()
            .map(|p| {
                let (pra, pdec) = candidate.pixel_to_world(p.px, p.py);
                crate::angle::angular_sep(
                    radec_to_vec(pra, pdec),
                    radec_to_vec(p.ra, p.dec),
                )
            })
            .collect();
        rms_list.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let med = crate::angle::median(&rms_list);
        let cutoff = (med + 2.5 * 1.4826 * median_abs_dev(&rms_list, med))
            .max(crate::angle::arcsec2rad(0.5));
        let mut changed = false;
        for (i, p) in pairs.iter().enumerate() {
            if !active[i] {
                continue;
            }
            let (pra, pdec) = candidate.pixel_to_world(p.px, p.py);
            let r = crate::angle::angular_sep(
                radec_to_vec(pra, pdec),
                radec_to_vec(p.ra, p.dec),
            );
            if r > cutoff {
                active[i] = false;
                changed = true;
            }
        }
        let kept = active.iter().filter(|&&x| x).count();
        let (rms, mx) = residual_stats(&candidate, pairs, &active);
        stats = (rms, mx);
        wcs = Some(candidate);
        if !changed || kept < 4 {
            break;
        }
    }
    let wcs = wcs?;
    Some(FitResult {
        wcs,
        rms_arcsec: rad2arcsec(stats.0),
        max_arcsec: rad2arcsec(stats.1),
        used: active.iter().filter(|&&x| x).count(),
        inliers: active,
    })
}

fn median_abs_dev(sorted: &[f64], med: f64) -> f64 {
    let mut d: Vec<f64> = sorted.iter().map(|x| (x - med).abs()).collect();
    d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    crate::angle::median(&d)
}

fn residual_stats(wcs: &TanWcs, pairs: &[Correspondence], active: &[bool]) -> (f64, f64) {
    let mut sum2 = 0.0;
    let mut max: f64 = 0.0;
    let mut n = 0u32;
    for (i, p) in pairs.iter().enumerate() {
        if !active[i] {
            continue;
        }
        let (pra, pdec) = wcs.pixel_to_world(p.px, p.py);
        let r = crate::angle::angular_sep(
            radec_to_vec(pra, pdec),
            radec_to_vec(p.ra, p.dec),
        );
        sum2 += r * r;
        max = max.max(r);
        n += 1;
    }
    ((sum2 / n.max(1) as f64).sqrt(), max)
}

/// 在 TAN 解基础上追加二阶 SIP 拟合（对应点 ≥ 12 才有意义）。
/// 系数在零附近线性化，一轮正规方程即可。
pub fn fit_sip(wcs: &TanWcs, pairs: &[Correspondence], active: &[bool]) -> TanWcs {
    let n = 6;
    let mut ata = vec![0.0; n * n];
    let mut atb = vec![0.0; n];
    let mut used = 0usize;
    for (i, p) in pairs.iter().enumerate() {
        if !active[i] {
            continue;
        }
        let u = p.px - wcs.crpix[0];
        let v = p.py - wcs.crpix[1];
        let (xi_t, eta_t) = tan_project(wcs.crval[0], wcs.crval[1], p.ra, p.dec);
        let xi0 = wcs.cd[0][0] * u + wcs.cd[0][1] * v;
        let eta0 = wcs.cd[1][0] * u + wcs.cd[1][1] * v;
        // 残差关于 SIP 系数线性化：
        // dxi ≈ cd11·(A20 u² + A21 uv + A02 v²) + cd12·(B20 u² + B21 uv + B02 v²)
        let (u2, uv, v2) = (u * u, u * v, v * v);
        let basis_a = [u2, uv, v2];
        let basis_b = [u2, uv, v2];
        let row = [
            wcs.cd[0][0] * basis_a[0],
            wcs.cd[0][0] * basis_a[1],
            wcs.cd[0][0] * basis_a[2],
            wcs.cd[0][1] * basis_b[0],
            wcs.cd[0][1] * basis_b[1],
            wcs.cd[0][1] * basis_b[2],
        ];
        // 每颗星产生两行观测方程（ξ、η），各自贡献一次 AᵀA/Aᵀb，
        // 等价于 2N 行的线性系统。
        accumulate_ata(&mut ata, &row, n);
        accumulate_atb(&mut atb, &row, xi_t - xi0, n);
        let row2 = [
            wcs.cd[1][0] * basis_a[0],
            wcs.cd[1][0] * basis_a[1],
            wcs.cd[1][0] * basis_a[2],
            wcs.cd[1][1] * basis_b[0],
            wcs.cd[1][1] * basis_b[1],
            wcs.cd[1][1] * basis_b[2],
        ];
        accumulate_ata(&mut ata, &row2, n);
        accumulate_atb(&mut atb, &row2, eta_t - eta0, n);
        used += 1;
    }
    if used < 12 {
        return *wcs;
    }
    let Some(x) = solve_normal(ata, atb, n) else {
        return *wcs;
    };
    let mut out = *wcs;
    out.sip = SipCoeffs {
        a20: x[0],
        a21: x[1],
        a02: x[2],
        b20: x[3],
        b21: x[4],
        b02: x[5],
    };
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::angle::{arcsec2rad, deg2rad};

    fn synthetic_wcs() -> TanWcs {
        // 中心 M45 附近，比例尺 2 角秒/像元，滚动 12°，6000x4000。
        let scale = arcsec2rad(2.0);
        let roll = deg2rad(12.0);
        let cd = [
            [scale * roll.cos(), -scale * roll.sin()],
            [scale * roll.sin(), scale * roll.cos()],
        ];
        TanWcs {
            crpix: [3000.0, 2000.0],
            crval: [deg2rad(56.75), deg2rad(24.12)],
            cd,
            sip: SipCoeffs::default(),
        }
    }

    #[test]
    fn tan_roundtrip() {
        let (ra0, dec0) = (deg2rad(120.0), deg2rad(-30.0));
        // TAN 投影只在切点周围（角距 < 90°）有定义，测试点都取近切点。
        for (ra, dec) in [
            (120.5, -29.5),
            (119.2, -30.7),
            (120.01, -30.01),
            (118.0, -28.0),
        ] {
            let (xi, eta) = tan_project(ra0, dec0, deg2rad(ra), deg2rad(dec));
            let (r2, d2) = tan_deproject(ra0, dec0, xi, eta);
            assert!((r2 - deg2rad(ra)).abs() < 1e-12);
            assert!((d2 - deg2rad(dec)).abs() < 1e-12);
        }
    }

    #[test]
    fn tan_roundtrip_across_ra_zero() {
        // 切点靠近 RA=0，验证赤经环绕处的正反投影一致性。
        let (ra0, dec0) = (deg2rad(0.02), deg2rad(45.0));
        for (ra, dec) in [(359.9, 45.1), (0.05, 44.9), (0.0, 45.0)] {
            let (xi, eta) = tan_project(ra0, dec0, deg2rad(ra), deg2rad(dec));
            let (r2, d2) = tan_deproject(ra0, dec0, xi, eta);
            let expect_ra = if ra > 180.0 { deg2rad(ra) - crate::angle::TWO_PI } else { deg2rad(ra) };
            let diff = ((r2 - expect_ra + crate::angle::PI) % crate::angle::TWO_PI)
                - crate::angle::PI;
            assert!(diff.abs() < 1e-12, "ra wrap: {ra} -> {}", r2.to_degrees());
            assert!((d2 - deg2rad(dec)).abs() < 1e-12);
        }
    }

    #[test]
    fn pixel_world_roundtrip() {
        let w = synthetic_wcs();
        for (x, y) in [(0.0, 0.0), (3000.0, 2000.0), (5999.5, 3999.5), (123.4, 567.8)] {
            let (ra, dec) = w.pixel_to_world(x, y);
            let (x2, y2) = w.world_to_pixel(ra, dec);
            assert!((x2 - x).abs() < 1e-6, "x {x} -> {x2}");
            assert!((y2 - y).abs() < 1e-6, "y {y} -> {y2}");
        }
    }

    #[test]
    fn fit_recovers_known_wcs() {
        let truth = synthetic_wcs();
        // 生成 40 个网格对应点（带 0.05 像元噪声由调用侧负责，这里先无噪）。
        let mut pairs = Vec::new();
        for i in 0..8 {
            for j in 0..5 {
                let px = 300.0 + i as f64 * 700.0;
                let py = 250.0 + j as f64 * 750.0;
                let (ra, dec) = truth.pixel_to_world(px, py);
                pairs.push(Correspondence { px, py, ra, dec });
            }
        }
        let fit = fit_tan(&pairs, 6000.0, 4000.0).expect("fit ok");
        assert!(fit.rms_arcsec < 0.01, "rms={}", fit.rms_arcsec);
        assert!((fit.wcs.crval[0] - truth.crval[0]).abs() < 1e-9);
        assert!((fit.wcs.crval[1] - truth.crval[1]).abs() < 1e-9);
        let (sx, sy) = fit.wcs.pixel_scales_arcsec();
        assert!((sx - 2.0).abs() < 1e-6 && (sy - 2.0).abs() < 1e-6);
        assert!((fit.wcs.roll_deg() - 12.0).abs() < 1e-6);
    }

    #[test]
    fn fit_robust_to_outliers() {
        let truth = synthetic_wcs();
        let mut pairs = Vec::new();
        for i in 0..8 {
            for j in 0..5 {
                let px = 300.0 + i as f64 * 700.0;
                let py = 250.0 + j as f64 * 750.0;
                let (ra, dec) = truth.pixel_to_world(px, py);
                pairs.push(Correspondence { px, py, ra, dec });
            }
        }
        // 污染 6 个点：随机挪到别处
        for k in 0..6 {
            let p = &mut pairs[k * 5];
            p.ra = deg2rad(200.0 + k as f64);
            p.dec = deg2rad(5.0);
        }
        let fit = fit_tan(&pairs, 6000.0, 4000.0).expect("fit ok");
        assert!(fit.used >= 34, "应剔除污染点，实际保留 {}", fit.used);
        assert!(fit.rms_arcsec < 0.05);
    }
}

#[cfg(test)]
mod solve_normal_tests {
    use super::{accumulate_ata, accumulate_atb, solve_normal};

    #[test]
    fn known_3x3_system() {
        // 设计矩阵行 [1, x, y]，取 4 个点，手算 AᵀA 与 Aᵀb。
        let rows = [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0]];
        let y = [2.0, 5.0, 7.0, 10.0]; // y = 2 + 3x + 5y
        let n = 3;
        let mut ata = vec![0.0; n * n];
        let mut atb = vec![0.0; n];
        for (r, yv) in rows.iter().zip(y) {
            accumulate_ata(&mut ata, r, n);
            accumulate_atb(&mut atb, r, yv, n);
        }
        let x = solve_normal(ata, atb, n).unwrap();
        assert!((x[0] - 2.0).abs() < 1e-9, "got {:?}", x);
        assert!((x[1] - 3.0).abs() < 1e-9, "got {:?}", x);
        assert!((x[2] - 5.0).abs() < 1e-9, "got {:?}", x);
    }
}
