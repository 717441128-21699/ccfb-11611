/** 跨模块共享的类型定义。 */

export interface StarPoint {
  x: number;
  y: number;
  flux: number;
  area: number;
}

export interface SolveJson {
  center_ra_deg: number;
  center_dec_deg: number;
  pixel_scale_x: number;
  pixel_scale_y: number;
  roll_deg: number;
  fov_width_deg: number;
  fov_height_deg: number;
  fov_diag_deg: number;
  rms_arcsec: number;
  max_resid_arcsec: number;
  num_inliers: number;
  num_detected: number;
  elapsed_ms: number;
  fits_header: string;
  wcs: {
    crpix1: number;
    crpix2: number;
    crval1_deg: number;
    crval2_deg: number;
    cd11: number;
    cd12: number;
    cd21: number;
    cd22: number;
  };
  matched: Array<{
    cat_id: number;
    x: number;
    y: number;
    ra_deg: number;
    dec_deg: number;
    mag: number;
  }>;
}

export interface PlottedDso {
  name: string;
  kind: string;
  x: number;
  y: number;
  ra_deg: number;
  dec_deg: number;
  mag: number;
}
