//! 深空天体（星云、星团、星系）表与视场内筛选。
//!
//! 数据表本身（NGC/IC 等）由外部静态资源提供，这里只定义与 WCS 无关的
//! 纯逻辑：按解出的视场从全表筛出落在画幅内的天体，并给出像元坐标。
//!
//! 这样核心库不绑定某一份星表，前端/CLI 可各自加载精简的 DSO 数据。

use crate::wcs::TanWcs;

/// 深空天体类别。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DsoKind {
    /// 星系
    Galaxy,
    /// 发射/反射/暗星云
    Nebula,
    /// 行星状星云
    PlanetaryNebula,
    /// 疏散星团
    OpenCluster,
    /// 球状星团
    GlobularCluster,
    /// 超新星遗迹
    SupernovaRemnant,
    /// 其它/未分类
    Other,
}

impl DsoKind {
    /// 紧凑标识符（用于序列化）。
    pub fn code(self) -> &'static str {
        match self {
            DsoKind::Galaxy => "G",
            DsoKind::Nebula => "N",
            DsoKind::PlanetaryNebula => "PN",
            DsoKind::OpenCluster => "OC",
            DsoKind::GlobularCluster => "GC",
            DsoKind::SupernovaRemnant => "SNR",
            DsoKind::Other => "X",
        }
    }

    pub fn from_code(s: &str) -> Self {
        match s {
            "G" => DsoKind::Galaxy,
            "N" => DsoKind::Nebula,
            "PN" => DsoKind::PlanetaryNebula,
            "OC" => DsoKind::OpenCluster,
            "GC" => DsoKind::GlobularCluster,
            "SNR" => DsoKind::SupernovaRemnant,
            _ => DsoKind::Other,
        }
    }
}

/// 一个深空天体。角度在库边界用弧度，`size_*` 为视直径（弧度），未知填 0。
#[derive(Clone, Debug)]
pub struct Dso {
    /// 显示名（如 "M42 / Orion Nebula"）。
    pub name: String,
    pub ra: f64,
    pub dec: f64,
    pub kind: DsoKind,
    pub major_axis: f64,
    pub minor_axis: f64,
    /// 位置角（弧度，北转东为正，与 WCS roll 同约定），未知填 0。
    pub position_angle: f64,
    pub mag: f32,
}

/// 投影到画幅内的 DSO（像元坐标 + 原始数据）。
#[derive(Clone, Debug)]
pub struct PlottedDso {
    pub dso: Dso,
    pub x: f64,
    pub y: f64,
}

/// 解析紧凑 CSV：`name,ra_deg,dec_deg,kind,mag,major_arcmin,minor_arcmin,pa_deg`。
/// 带表头，`#` 开头为注释；尺寸/位置角可留空（按 0 处理）。
pub fn parse_dso_csv(text: &str) -> Result<Vec<Dso>, String> {
    let mut out = Vec::new();
    for (lineno, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let lower = line.to_ascii_lowercase();
        if lower.starts_with("name,") {
            continue;
        }
        let c: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
        if c.len() < 5 {
            return Err(format!("line {}: 至少需要 name,ra,dec,kind,mag", lineno + 1));
        }
        let ra = c[1]
            .parse::<f64>()
            .map_err(|e| format!("line {} ra: {e}", lineno + 1))?
            .to_radians();
        let dec = c[2]
            .parse::<f64>()
            .map_err(|e| format!("line {} dec: {e}", lineno + 1))?
            .to_radians();
        if !(0.0..360.0).contains(&c[1].parse::<f64>().unwrap_or(-1.0))
            || !(-90.0..=90.0).contains(&c[2].parse::<f64>().unwrap_or(-99.0))
        {
            return Err(format!("line {}: ra/dec 越界", lineno + 1));
        }
        let mag = c[4].parse::<f32>().unwrap_or(f32::NAN);
        let major_arcmin = c.get(5).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        let minor_arcmin = c.get(6).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        let pa_deg = c.get(7).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        out.push(Dso {
            name: c[0].to_string(),
            ra,
            dec,
            kind: DsoKind::from_code(c[3]),
            major_axis: crate::angle::arcsec2rad(major_arcmin * 60.0),
            minor_axis: crate::angle::arcsec2rad(minor_arcmin * 60.0),
            position_angle: pa_deg.to_radians(),
            mag,
        });
    }
    Ok(out)
}

/// 筛出落在画幅内（含 `margin` 像元容差）的 DSO，并投影到像元坐标。
pub fn dso_in_frame(
    wcs: &TanWcs,
    objects: &[Dso],
    width: f64,
    height: f64,
    margin: f64,
) -> Vec<PlottedDso> {
    let mut out = Vec::new();
    for dso in objects {
        // 先用角距离粗筛：视场半对角 + 天体半径，避免对全表做投影。
        let (cx_ra, cx_dec) = wcs.pixel_to_world(width * 0.5, height * 0.5);
        let sep = crate::angle::angular_sep(
            crate::angle::radec_to_vec(cx_ra, cx_dec),
            crate::angle::radec_to_vec(dso.ra, dso.dec),
        );
        let corner = wcs.pixel_to_world(0.0, 0.0);
        let half_diag = crate::angle::angular_sep(
            crate::angle::radec_to_vec(cx_ra, cx_dec),
            crate::angle::radec_to_vec(corner.0, corner.1),
        );
        if sep > half_diag + dso.major_axis + crate::angle::arcsec2rad(margin) {
            continue;
        }
        let (x, y) = wcs.world_to_pixel(dso.ra, dso.dec);
        if x >= -margin && x <= width + margin && y >= -margin && y <= height + margin {
            out.push(PlottedDso { dso: dso.clone(), x, y });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::angle::{arcsec2rad, deg2rad};
    use crate::wcs::SipCoeffs;

    fn wcs() -> TanWcs {
        let scale = arcsec2rad(9.0);
        let roll = deg2rad(10.0);
        TanWcs {
            crpix: [500.0, 400.0],
            crval: [deg2rad(180.0), deg2rad(0.0)],
            cd: [
                [scale * roll.cos(), -scale * roll.sin()],
                [scale * roll.sin(), scale * roll.cos()],
            ],
            sip: SipCoeffs::default(),
        }
    }

    #[test]
    fn csv_roundtrip_and_kind_codes() {
        let csv = "name,ra_deg,dec_deg,kind,mag,major_arcmin,minor_arcmin,pa_deg\n\
                   M42,83.82,-5.39,N,4.0,85,60,-40\n\
                   M31,10.68,41.27,G,3.4,178,63,35\n";
        let v = parse_dso_csv(csv).unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].kind, DsoKind::Nebula);
        assert_eq!(DsoKind::Galaxy.code(), "G");
        assert!(v[1].major_axis > 0.0);
    }

    #[test]
    fn filters_objects_inside_frame() {
        let w = wcs();
        // 中心处放一个天体，画幅外再放一个。
        let (cra, cdec) = w.pixel_to_world(500.0, 400.0);
        let inside = Dso {
            name: "center".into(), ra: cra, dec: cdec, kind: DsoKind::Nebula,
            major_axis: 0.0, minor_axis: 0.0, position_angle: 0.0, mag: 8.0,
        };
        let outside = Dso {
            name: "far".into(), ra: deg2rad(10.0), dec: deg2rad(80.0),
            kind: DsoKind::Galaxy, major_axis: 0.0, minor_axis: 0.0,
            position_angle: 0.0, mag: 9.0,
        };
        let plotted = dso_in_frame(&w, &[inside, outside], 1000.0, 800.0, 5.0);
        assert_eq!(plotted.len(), 1);
        assert_eq!(plotted[0].dso.name, "center");
        assert!((plotted[0].x - 500.0).abs() < 1e-4);
        assert!((plotted[0].y - 400.0).abs() < 1e-4);
    }
}
