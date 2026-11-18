//! Focused diagnostic: do image quads and catalog tangent-plane quads for the
//! same four physical stars hash alike?

use platersol_core::astro::{TanProjection, ARCSEC, DEG};
use platersol_core::quad::{make_quad, Pt};
use platersol_core::{render_synthetic_image, SyntheticCamera, SyntheticSky};

#[test]
fn same_four_stars_hash_consistency() {
    let sky = SyntheticSky::random(300_000, 3.0, 9.0, 20261118);
    let cam = SyntheticCamera {
        width: 2048,
        height: 1536,
        ra0: 10.68 * DEG,
        dec0: 41.27 * DEG,
        scale: 6.0 * ARCSEC,
        rotation: 57.0 * DEG,
        psf_sigma_px: 1.6,
        sky: 500.0,
        read_noise: 6.0,
        zp_flux: 4.0e7,
        k1: 0.0,
    };
    let (_img, injected) = render_synthetic_image(&sky, &cam, 9917);
    assert!(injected.len() >= 8, "only {} injected", injected.len());

    // Pick the four brightest injected stars with good frame separation.
    let mut stars: Vec<_> = injected.iter().collect();
    stars.sort_by(|a, b| b.0.mag.partial_cmp(&a.0.mag).unwrap()); // mag asc = brighter first
    let chosen = &stars[0..4];

    // Image-space quad in raw pixels.
    let img_pts: Vec<Pt> = chosen
        .iter()
        .enumerate()
        .map(|(i, (_, x, y))| Pt {
            x: *x,
            y: *y,
            id: i as u32,
        })
        .collect();
    let img_quad = make_quad([img_pts[0], img_pts[1], img_pts[2], img_pts[3]]).unwrap();

    // Catalog-space quad: project the same stars onto the camera-center
    // tangent plane (the gnomonic coordinates the index builder would use).
    let proj = TanProjection::new(cam.ra0, cam.dec0);
    let cat_pts: Vec<Pt> = chosen
        .iter()
        .enumerate()
        .map(|(i, (s, _, _))| {
            let (xi, eta) = proj.project(s.ra, s.dec).unwrap();
            Pt {
                x: xi,
                y: eta,
                id: i as u32,
            }
        })
        .collect();
    let cat_quad = make_quad([cat_pts[0], cat_pts[1], cat_pts[2], cat_pts[3]]).unwrap();

    eprintln!("img AB (px) = {:.3}  cat AB (rad) = {:.3e}", img_quad.ab, cat_quad.ab);
    eprintln!("implied plate scale (arcsec/px) = {:.4}",
        (cat_quad.ab / img_quad.ab) / ARCSEC);
    eprintln!("img c={:?} d={:?}", img_quad.c, img_quad.d);
    eprintln!("cat c={:?} d={:?}", cat_quad.c, cat_quad.d);
    eprintln!("img hash={:016x} mirror={:016x}", img_quad.hash, img_quad.hash_mirror);
    eprintln!("cat hash={:016x} mirror={:016x}", cat_quad.hash, cat_quad.hash_mirror);

    // The canonical inner coordinates must agree closely; the hash is a
    // quantized version of exactly these normalized coordinates.
    let dc = (img_quad.c.0 - cat_quad.c.0).hypot(img_quad.c.1 - cat_quad.c.1);
    let dd = (img_quad.d.0 - cat_quad.d.0).hypot(img_quad.d.1 - cat_quad.d.1);
    eprintln!("canonical coord diff c={:.3e} d={:.3e}", dc, dd);
    assert!(dc < 0.02, "canonical C mismatch {dc}");
    assert!(dd < 0.02, "canonical D mismatch {dd}");
}
