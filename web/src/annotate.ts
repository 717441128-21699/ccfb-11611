/**
 * 解算结果的图像标注：星点圈选、亮星名称、星座连线。
 *
 * 星座连线这里用"亮星邻接图"做几何示意（最近邻连线），不依赖外部星座
 * 边界数据；接入真实星名/星座库后可替换 nameOf 与 constellationLines。
 */
import type { PlottedDso, SolveJson } from "./types";
export interface AnnotationOptions {
  showStars: boolean;
  showNames: boolean;
  showLines: boolean;
  showDso: boolean;
  /** 参与命名/连线的最亮星数量 */
  topBright: number;
  /** 星座连线的最近邻搜索半径（像元，降采样坐标） */
  linkRadius: number;
  nameOf?: (catId: number, ra: number, dec: number) => string | null;
}

/** DSO 类别 → 标注颜色。 */
const DSO_COLORS: Record<string, string> = {
  G: "#ff7eb6", // 星系：粉
  N: "#7CFFB2", // 星云：绿
  PN: "#9dff7c",
  OC: "#ffd24a", // 疏散星团：金
  GC: "#ffb86b", // 球状星团：橙
  SNR: "#ff6b6b",
  X: "#c9b8ff",
};

interface PlottedStar {
  x: number;
  y: number;
  mag: number;
  catId: number;
  ra: number;
  dec: number;
}

/** 在指定 canvas 上绘制标注。canvas 与降采样图同尺寸。 */
export function drawAnnotations(
  ctx: CanvasRenderingContext2D,
  width: number,
  height: number,
  solution: SolveJson,
  dsoList: PlottedDso[],
  opts: AnnotationOptions,
): void {
  ctx.clearRect(0, 0, width, height);

  // 仅标注匹配上的星（有星表身份），按星等升序（亮的在前）。
  const matched: PlottedStar[] = solution.matched
    .map((m) => ({
      x: m.x,
      y: m.y,
      mag: m.mag,
      catId: m.cat_id,
      ra: m.ra_deg,
      dec: m.dec_deg,
    }))
    .sort((a, b) => a.mag - b.mag);

  const bright = matched.slice(0, opts.topBright);

  if (opts.showLines) {
    drawConstellationLines(ctx, bright, opts.linkRadius);
  }
  if (opts.showDso) {
    drawDsoMarkers(ctx, dsoList);
  }
  if (opts.showStars) {
    drawStarMarkers(ctx, matched, bright.length);
  }
  if (opts.showNames) {
    drawNames(ctx, bright, opts.nameOf);
  }
}

function drawStarMarkers(
  ctx: CanvasRenderingContext2D,
  stars: PlottedStar[],
  brightCount: number,
): void {
  ctx.lineWidth = 1;
  stars.forEach((s, i) => {
    const isBright = i < brightCount;
    const r = isBright ? 6 : 3.5;
    ctx.beginPath();
    ctx.arc(s.x, s.y, r, 0, Math.PI * 2);
    ctx.strokeStyle = isBright ? "#ffd24a" : "#57d9ff";
    ctx.stroke();
    if (isBright) {
      ctx.beginPath();
      ctx.arc(s.x, s.y, 1.4, 0, Math.PI * 2);
      ctx.fillStyle = "#ffd24a";
      ctx.fill();
    }
  });
}

function drawNames(
  ctx: CanvasRenderingContext2D,
  stars: PlottedStar[],
  nameOf?: AnnotationOptions["nameOf"],
): void {
  ctx.font = "12px system-ui, sans-serif";
  ctx.textBaseline = "middle";
  stars.forEach((s) => {
    const label = nameOf ? nameOf(s.catId, s.ra, s.dec) : null;
    if (!label) return;
    // 文字描边底，保证在亮/暗背景上都可读
    ctx.lineWidth = 3;
    ctx.strokeStyle = "rgba(0,0,0,0.75)";
    ctx.strokeText(label, s.x + 8, s.y - 6);
    ctx.fillStyle = "#fff";
    ctx.fillText(label, s.x + 8, s.y - 6);
  });
}

/**
 * 深空天体标注：星系用菱形、星云/星团用方框，配类别配色与名称。
 */
function drawDsoMarkers(ctx: CanvasRenderingContext2D, objects: PlottedDso[]): void {
  const r = 8;
  ctx.lineWidth = 1.6;
  ctx.font = "bold 12px system-ui, sans-serif";
  ctx.textBaseline = "middle";
  objects.forEach((o) => {
    const color = DSO_COLORS[o.kind] ?? DSO_COLORS.X;
    ctx.strokeStyle = color;
    ctx.beginPath();
    if (o.kind === "G") {
      // 菱形
      ctx.moveTo(o.x, o.y - r);
      ctx.lineTo(o.x + r, o.y);
      ctx.lineTo(o.x, o.y + r);
      ctx.lineTo(o.x - r, o.y);
      ctx.closePath();
    } else {
      // 方框
      ctx.rect(o.x - r, o.y - r, r * 2, r * 2);
    }
    ctx.stroke();
    // 名称带描边底
    ctx.lineWidth = 3;
    ctx.strokeStyle = "rgba(0,0,0,0.8)";
    ctx.strokeText(o.name, o.x + r + 3, o.y);
    ctx.fillStyle = color;
    ctx.fillText(o.name, o.x + r + 3, o.y);
    ctx.lineWidth = 1.6;
  });
}

/**
 * 几何星座示意：每颗亮星连向半径内最近的两颗亮星，形成不交叉太多的网状结构。
 * 这是可视化占位，不是 IAU 官方星座连线。
 */
function drawConstellationLines(
  ctx: CanvasRenderingContext2D,
  stars: PlottedStar[],
  radius: number,
): void {
  const r2 = radius * radius;
  ctx.strokeStyle = "rgba(120, 200, 255, 0.55)";
  ctx.lineWidth = 1;
  const drawn = new Set<string>();
  stars.forEach((a, i) => {
    const near = stars
      .map((b, j) => ({ j, d2: (b.x - a.x) ** 2 + (b.y - a.y) ** 2 }))
      .filter((o) => o.j !== i && o.d2 <= r2)
      .sort((p, q) => p.d2 - q.d2)
      .slice(0, 2);
    near.forEach(({ j }) => {
      const key = i < j ? `${i}-${j}` : `${j}-${i}`;
      if (drawn.has(key)) return;
      drawn.add(key);
      const b = stars[j];
      ctx.beginPath();
      ctx.moveTo(a.x, a.y);
      ctx.lineTo(b.x, b.y);
      ctx.stroke();
    });
  });
}
