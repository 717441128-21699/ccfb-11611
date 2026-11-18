// Overlay annotation renderer. Given a fitted WCS exposed by the WASM engine,
// projects reference sky positions (named stars, constellation lines, DSOs)
// back into image pixel space and draws labels/markers on a transparent
// canvas stacked over the photo.

import type { SolveResult } from "../wasm/loader";
import {
  CONSTELLATIONS,
  DEEP_SKY,
  NAMED_STARS,
  type Constellation,
  type DeepSkyObject,
  type NamedStar,
} from "./reference-data";

export interface OverlayOptions {
  showStars: boolean;
  showConstellations: boolean;
  showDso: boolean;
}

interface Projector {
  // Sky RA/Dec degrees -> image pixel, or null if outside the frame.
  radecToPixel(raDeg: number, decDeg: number): { x: number; y: number } | null;
}

const D2R = Math.PI / 180;
const R2D = 180 / Math.PI;

/**
 * Build a projector from the SolveResult by reconstructing the tangent-plane
 * affine (the same model the Rust WCS fitter uses). This mirrors the f64
 * coefficients returned via ps_wcs_coeffs but here we derive enough from the
 * reported center, scale, rotation and image dimensions for annotation.
 */
export function makeProjector(
  result: SolveResult,
  imgW: number,
  imgH: number,
): Projector {
  const ra0 = result.ra0Deg * D2R;
  const dec0 = result.dec0Deg * D2R;
  const scale = result.plateScaleArcsec * (Math.PI / (180 * 3600)); // rad/px
  const rot = result.rotationDeg * D2R;
  const sinD0 = Math.sin(dec0);
  const cosD0 = Math.cos(dec0);

  return {
    radecToPixel(raDeg, decDeg) {
      const ra = raDeg * D2R;
      const dec = decDeg * D2R;
      // Gnomonic projection onto the solution tangent plane.
      const sd = Math.sin(dec);
      const cd = Math.cos(dec);
      const dra = ra - ra0;
      const c = sinD0 * sd + cosD0 * cd * Math.cos(dra);
      if (c <= 1e-9) return null;
      const xi = (cd * Math.sin(dra)) / c;
      const eta = (cosD0 * sd - sinD0 * cd * Math.cos(dra)) / c;
      // Invert the rotated+scaled pixel affine: tangent coords -> pixels.
      const cs = Math.cos(rot);
      const sn = Math.sin(rot);
      const dx = (xi * cs + eta * sn) / scale;
      const dy = (-xi * sn + eta * cs) / scale;
      const x = dx + imgW * 0.5;
      const y = dy + imgH * 0.5;
      if (x < 0 || y < 0 || x >= imgW || y >= imgH) return null;
      return { x, y };
    },
  };
}

export function renderOverlay(
  ctx: CanvasRenderingContext2D,
  projector: Projector,
  imgW: number,
  imgH: number,
  opts: OverlayOptions,
): void {
  ctx.clearRect(0, 0, imgW, imgH);
  ctx.save();
  ctx.lineJoin = "round";

  if (opts.showConstellations) {
    for (const c of CONSTELLATIONS) drawConstellation(ctx, projector, c);
  }
  if (opts.showStars) {
    for (const s of NAMED_STARS) drawNamedStar(ctx, projector, s);
  }
  if (opts.showDso) {
    for (const d of DEEP_SKY) drawDso(ctx, projector, d);
  }
  ctx.restore();
}

function drawNamedStar(ctx: CanvasRenderingContext2D, p: Projector, s: NamedStar): void {
  const pt = p.radecToPixel(s.ra, s.dec);
  if (!pt) return;
  const r = Math.max(2.5, 6 - s.mag * 0.5);
  ctx.beginPath();
  ctx.arc(pt.x, pt.y, r, 0, Math.PI * 2);
  ctx.fillStyle = "rgba(255, 244, 200, 0.95)";
  ctx.shadowColor = "rgba(255, 230, 150, 0.9)";
  ctx.shadowBlur = 8;
  ctx.fill();
  ctx.shadowBlur = 0;
  label(ctx, pt.x + r + 4, pt.y - r, s.name, "#ffe9b0");
}

function drawConstellation(ctx: CanvasRenderingContext2D, p: Projector, c: Constellation): void {
  ctx.strokeStyle = "rgba(124, 200, 255, 0.55)";
  ctx.lineWidth = 1.4;
  ctx.setLineDash([6, 5]);
  for (const [[ra1, dec1], [ra2, dec2]] of c.lines) {
    const a = p.radecToPixel(ra1, dec1);
    const b = p.radecToPixel(ra2, dec2);
    if (!a || !b) continue;
    ctx.beginPath();
    ctx.moveTo(a.x, a.y);
    ctx.lineTo(b.x, b.y);
    ctx.stroke();
  }
  ctx.setLineDash([]);
}

function drawDso(ctx: CanvasRenderingContext2D, p: Projector, d: DeepSkyObject): void {
  const pt = p.radecToPixel(d.ra, d.dec);
  if (!pt) return;
  const color =
    d.type === "galaxy" ? "#9d8cff" : d.type === "nebula" ? "#7ce0c8" : "#ffb27a";
  ctx.beginPath();
  ctx.arc(pt.x, pt.y, 9, 0, Math.PI * 2);
  ctx.strokeStyle = color;
  ctx.lineWidth = 1.6;
  ctx.stroke();
  // Crosshair marker.
  ctx.beginPath();
  ctx.moveTo(pt.x - 13, pt.y); ctx.lineTo(pt.x + 13, pt.y);
  ctx.moveTo(pt.x, pt.y - 13); ctx.lineTo(pt.x, pt.y + 13);
  ctx.strokeStyle = color;
  ctx.globalAlpha = 0.7;
  ctx.stroke();
  ctx.globalAlpha = 1;
  label(ctx, pt.x + 14, pt.y + 4, `${d.name}`, color);
}

function label(ctx: CanvasRenderingContext2D, x: number, y: number, text: string, color: string): void {
  ctx.font = "13px 'Segoe UI', system-ui, sans-serif";
  const metrics = ctx.measureText(text);
  ctx.fillStyle = "rgba(5, 8, 16, 0.72)";
  ctx.fillRect(x - 3, y - 12, metrics.width + 6, 17);
  ctx.fillStyle = color;
  ctx.fillText(text, x, y);
}

export { R2D };
