/**
 * 应用主逻辑：建库 → 导入图像 → 检测 → 盲解 → 标注与结果展示。
 * 所有重计算经 WASM 在本地完成。
 */
import { decodeToGray, type GrayBuffer } from "./image";
import {
  ensureWasm,
  buildIndex,
  detectStars,
  backgroundStats,
  solveImage,
  selfTest,
  projectDso,
  type SolveJson,
} from "./solver";
import { generateDemoCsv, DEFAULT_DEMO } from "./demoCatalog";
import { drawAnnotations } from "./annotate";
import dsoCsv from "./dso-sample.csv?raw";
import type { PlottedDso, StarPoint } from "./ui-types";
import type { WasmIndex } from "../pkg/star_solver_wasm.js";

const $ = <T extends HTMLElement>(id: string): T =>
  document.getElementById(id) as T;

const els = {
  btnDemo: $<HTMLButtonElement>("btn-demo-index"),
  fileCatalog: $<HTMLInputElement>("file-catalog"),
  magLimit: $<HTMLInputElement>("mag-limit"),
  indexStatus: $<HTMLParagraphElement>("index-status"),
  fileImage: $<HTMLInputElement>("file-image"),
  imageStatus: $<HTMLParagraphElement>("image-status"),
  scaleLo: $<HTMLInputElement>("scale-lo"),
  scaleHi: $<HTMLInputElement>("scale-hi"),
  kSigma: $<HTMLInputElement>("k-sigma"),
  maxStars: $<HTMLInputElement>("max-stars"),
  optStars: $<HTMLInputElement>("opt-stars"),
  optNames: $<HTMLInputElement>("opt-names"),
  optLines: $<HTMLInputElement>("opt-lines"),
  optDso: $<HTMLInputElement>("opt-dso"),
  btnSolve: $<HTMLButtonElement>("btn-solve"),
  btnSelftest: $<HTMLButtonElement>("btn-selftest"),
  runStatus: $<HTMLParagraphElement>("run-status"),
  base: $<HTMLCanvasElement>("base-canvas"),
  overlay: $<HTMLCanvasElement>("overlay-canvas"),
  result: $<HTMLDivElement>("result"),
};

let index: WasmIndex | null = null;
let gray: GrayBuffer | null = null;
let detected: StarPoint[] = [];
let lastSolution: SolveJson | null = null;
let plottedDso: PlottedDso[] = [];

function setStatus(el: HTMLParagraphElement, msg: string, kind: "" | "ok" | "err" = "") {
  el.textContent = msg;
  el.className = kind ? `status ${kind}` : "status";
}

function canSolve(): void {
  els.btnSolve.disabled = !(index && gray);
}

// ---- 星表建库 ----

async function buildFromCsv(csv: string, source: string): Promise<void> {
  setStatus(els.indexStatus, `正在为 ${source} 构建索引（首次较慢）……`);
  const t0 = performance.now();
  try {
    // 建库是 CPU 密集型，让出一帧让状态文字先渲染。
    await new Promise((r) => setTimeout(r, 30));
    const magLimit = parseFloat(els.magLimit.value) || 9.5;
    index = await buildIndex(csv, magLimit);
    const secs = ((performance.now() - t0) / 1000).toFixed(1);
    setStatus(
      els.indexStatus,
      `索引就绪：${index.star_count.toLocaleString()} 颗星、${index.band_count} 个尺度带，用时 ${secs}s`,
      "ok",
    );
    canSolve();
  } catch (e) {
    setStatus(els.indexStatus, `建库失败：${(e as Error).message}`, "err");
  }
}

els.btnDemo.addEventListener("click", () => {
  const csv = generateDemoCsv(DEFAULT_DEMO);
  void buildFromCsv(csv, `内置演示星表（${DEFAULT_DEMO.count} 颗）`);
});

els.fileCatalog.addEventListener("change", async () => {
  const file = els.fileCatalog.files?.[0];
  if (!file) return;
  const text = await file.text();
  await buildFromCsv(text, file.name);
});

// ---- 图像导入 ----

els.fileImage.addEventListener("change", async () => {
  const file = els.fileImage.files?.[0];
  if (!file) return;
  setStatus(els.imageStatus, "正在解码图像……");
  try {
    gray = await decodeToGray(file);
    paintBase(gray);
    const bg = await backgroundStats(
      Uint8Array.from(gray.data),
      gray.width,
      gray.height,
    );
    setStatus(
      els.imageStatus,
      `已载入 ${gray.originalWidth}×${gray.originalHeight}` +
        (gray.scale > 1 ? `（检测缩至 ${gray.width}×${gray.height}）` : "") +
        `，背景 ${bg.level.toFixed(1)} ± ${bg.sigma.toFixed(2)}`,
      "ok",
    );
    canSolve();
  } catch (e) {
    setStatus(els.imageStatus, `图像解码失败：${(e as Error).message}`, "err");
  }
});

function paintBase(g: GrayBuffer): void {
  els.base.width = g.width;
  els.base.height = g.height;
  els.overlay.width = g.width;
  els.overlay.height = g.height;
  const ctx = els.base.getContext("2d")!;
  // 复制到全新 ArrayBuffer 支撑的数组，满足 ImageData 的类型约束。
  const img = new ImageData(new Uint8ClampedArray(g.data), g.width, g.height);
  ctx.putImageData(img, 0, 0);
}

// ---- 解算 ----

els.btnSolve.addEventListener("click", async () => {
  if (!index || !gray) return;
  setStatus(els.runStatus, "星点检测中……");
  els.btnSolve.disabled = true;
  await new Promise((r) => setTimeout(r, 30));
  try {
    const g = gray;
    const kSigma = parseFloat(els.kSigma.value) || 5;
    const maxStars = parseInt(els.maxStars.value, 10) || 500;
    detected = await detectStars(
      Uint8Array.from(g.data),
      g.width,
      g.height,
      kSigma,
      maxStars,
    );
    if (detected.length < 4) {
      setStatus(els.runStatus, `仅检出 ${detected.length} 颗星，无法建 quad`, "err");
      els.btnSolve.disabled = false;
      return;
    }
    setStatus(els.runStatus, `检出 ${detected.length} 颗星，几何哈希盲解中……`);
    await new Promise((r) => setTimeout(r, 30));

    const res = await solveImage(
      index,
      Uint8Array.from(g.data),
      g.width,
      g.height,
      parseFloat(els.scaleLo.value) || 0.5,
      parseFloat(els.scaleHi.value) || 25,
    );
    if (!res.ok || !res.solution) {
      setStatus(els.runStatus, `解算失败：${res.error}`, "err");
      els.btnSolve.disabled = false;
      return;
    }
    lastSolution = res.solution;
    // 用解出的 WCS 把内置 DSO 表投影到当前画幅。
    try {
      plottedDso = await projectDso(
        JSON.stringify(lastSolution.wcs),
        dsoCsv,
        g.width,
        g.height,
      );
    } catch {
      plottedDso = [];
    }
    renderResult(lastSolution);
    redrawOverlay();
    setStatus(
      els.runStatus,
      `解算成功：内点 ${lastSolution.num_inliers}/${lastSolution.num_detected}，RMS ${lastSolution.rms_arcsec.toFixed(2)}″，用时 ${lastSolution.elapsed_ms}ms`,
      "ok",
    );
  } catch (e) {
    setStatus(els.runStatus, `出错：${(e as Error).message}`, "err");
  } finally {
    els.btnSolve.disabled = !(index && gray);
  }
});

// ---- 结果与标注 ----

function renderResult(s: SolveJson): void {
  els.result.classList.remove("hidden");
  els.result.innerHTML = `
    <h3>解算结果（ICRS / J2000）</h3>
    <dl class="kv">
      <dt>中心赤经 RA</dt><dd>${s.center_ra_deg.toFixed(6)}°</dd>
      <dt>中心赤纬 DEC</dt><dd>${s.center_dec_deg.toFixed(6)}°</dd>
      <dt>比例尺</dt><dd>${s.pixel_scale_x.toFixed(3)} × ${s.pixel_scale_y.toFixed(3)} 角秒/像元</dd>
      <dt>旋转角</dt><dd>${s.roll_deg.toFixed(3)}°（CROTA2，北转东为正）</dd>
      <dt>视场</dt><dd>${s.fov_width_deg.toFixed(3)}° × ${s.fov_height_deg.toFixed(3)}°（对角 ${s.fov_diag_deg.toFixed(3)}°）</dd>
      <dt>匹配内点</dt><dd>${s.num_inliers} / ${s.num_detected} 颗</dd>
      <dt>残差 RMS</dt><dd>${s.rms_arcsec.toFixed(3)}″（最大 ${s.max_resid_arcsec.toFixed(3)}″）</dd>
    </dl>
    <pre class="fits">${escapeHtml(s.fits_header)}</pre>`;
}

function redrawOverlay(): void {
  if (!gray || !lastSolution) return;
  const ctx = els.overlay.getContext("2d")!;
  ctx.clearRect(0, 0, gray.width, gray.height);
  drawAnnotations(ctx, gray.width, gray.height, lastSolution, plottedDso, {
    showStars: els.optStars.checked,
    showNames: els.optNames.checked,
    showLines: els.optLines.checked,
    showDso: els.optDso.checked,
    topBright: 24,
    linkRadius: Math.max(gray.width, gray.height) * 0.12,
    // 暂无外部星名库：用编号+星等占位，接入真实命名表后替换。
    nameOf: (_id, _ra, _dec) => null,
  });
}

[els.optStars, els.optNames, els.optLines, els.optDso].forEach((c) =>
  c.addEventListener("change", redrawOverlay),
);

// ---- 引擎自检 ----

els.btnSelftest.addEventListener("click", async () => {
  setStatus(els.runStatus, "正在运行合成星空端到端自检（约 1~2 分钟）……");
  await new Promise((r) => setTimeout(r, 30));
  const r = await selfTest();
  if (r.ok) {
    setStatus(els.runStatus, `引擎自检通过：${JSON.stringify(r.report)}`, "ok");
  } else {
    setStatus(els.runStatus, `引擎自检失败：${r.error}`, "err");
  }
});

function escapeHtml(s: string): string {
  return s
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;");
}

// 启动即初始化 WASM。
void ensureWasm().then(() => setStatus(els.runStatus, "WASM 引擎已就绪，全程离线运行"));
