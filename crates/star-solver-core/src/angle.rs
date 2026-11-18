//! 角度/球面小工具。库内统一弧度制。

pub const TWO_PI: f64 = std::f64::consts::TAU;
pub const PI: f64 = std::f64::consts::PI;

/// 角度（度）→ 弧度。
#[inline]
pub fn deg2rad(d: f64) -> f64 {
    d.to_radians()
}

/// 弧度 → 角度（度）。
#[inline]
pub fn rad2deg(r: f64) -> f64 {
    r.to_degrees()
}

/// 弧度 → 角秒。
#[inline]
pub fn rad2arcsec(r: f64) -> f64 {
    r.to_degrees() * 3600.0
}

/// 角秒 → 弧度。
#[inline]
pub fn arcsec2rad(a: f64) -> f64 {
    (a / 3600.0).to_radians()
}

/// 赤经归一化到 [0, 2π)。
#[inline]
pub fn wrap_ra(ra: f64) -> f64 {
    let mut r = ra % TWO_PI;
    if r < 0.0 {
        r += TWO_PI;
    }
    r
}

/// 赤纬/角度归一化到 (-π, π]。
#[inline]
pub fn wrap_angle(a: f64) -> f64 {
    let mut r = a % TWO_PI;
    if r <= -PI {
        r += TWO_PI;
    } else if r > PI {
        r -= TWO_PI;
    }
    r
}

/// 单位球向量（J2000 笛卡尔：x 指向春分点，z 指向北天极）。
pub type Vec3 = [f64; 3];

/// 赤经赤纬（弧度）→ 单位向量。
#[inline]
pub fn radec_to_vec(ra: f64, dec: f64) -> Vec3 {
    let (sd, cd) = dec.sin_cos();
    let (sr, cr) = ra.sin_cos();
    [cd * cr, cd * sr, sd]
}

/// 单位向量 → 赤经赤纬（弧度）。
#[inline]
pub fn vec_to_radec(v: Vec3) -> (f64, f64) {
    let dec = v[2].asin().clamp(-PI / 2.0, PI / 2.0);
    let ra = wrap_ra(v[1].atan2(v[0]));
    (ra, dec)
}

#[inline]
pub fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
pub fn normalize(v: Vec3) -> Vec3 {
    let n = (dot(v, v)).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}

/// 两点角距离（弧度），输入单位向量。
#[inline]
pub fn angular_sep(a: Vec3, b: Vec3) -> f64 {
    dot(a, b).clamp(-1.0, 1.0).acos()
}

/// Rodrigues 旋转：把向量 v 绕单位轴 k 转 theta。
pub fn rotate_about(v: Vec3, k: Vec3, theta: f64) -> Vec3 {
    let (s, c) = theta.sin_cos();
    let kv = dot(k, v);
    let kxv = cross(k, v);
    [
        v[0] * c + kxv[0] * s + k[0] * kv * (1.0 - c),
        v[1] * c + kxv[1] * s + k[1] * kv * (1.0 - c),
        v[2] * c + kxv[2] * s + k[2] * kv * (1.0 - c),
    ]
}

/// 升序切片的中位数（f64, NaN 行为不做特殊处理——本项目数据不含 NaN）。
pub fn median(sorted: &[f64]) -> f64 {
    debug_assert!(sorted.windows(2).all(|w| w[0] <= w[1]));
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        0.5 * (sorted[n / 2 - 1] + sorted[n / 2])
    }
}
