// Main UI controller: image loading, WASM engine bootstrap, worker lifecycle,
// result display and overlay coordination. All heavy work runs in the worker.

import { decodeToGray, type GrayImage } from "./imaging/decode";
import { makeProjector, renderOverlay, type OverlayOptions } from "./sky/overlay";
import type { SolveResult, SourceC } from "./wasm/loader";
import type { CatStarRow } from "./wasm/chunk-store";
import SolverWorker from "./workers/solver.worker?worker";

const $ = <T extends HTMLElement>(id: string): T => document.getElementById(id) as T;

const els = {
  engine: $<HTMLDivElement>("engineStatus"),
  fileInput: $<HTMLInputElement>("fileInput"),
  drop: $<HTMLDivElement>("dropHint"),
  wrap: $<HTMLDivElement>("canvasWrap"),
  photo: $<HTMLCanvasElement>("photo"),
  overlay: $<HTMLCanvasElement>("overlay"),
  solveBtn: $<HTMLButtonElement>("solveBtn"),
  threshold: $<HTMLInputElement>("threshold"),
  hashTol: $<HTMLInputElement>("hashTol"),
  quorum: $<HTMLInputElement>("quorum"),
  indexUrl: $<HTMLInputElement>("indexUrl"),
  progressBlock: $<HTMLDivElement>("progressBlock"),
  stageText: $<HTMLDivElement>("stageText"),
  fill: $<HTMLDivElement>("progressFill"),
  resultBlock: $<HTMLDivElement>("resultBlock"),
  ra: $<HTMLElement>("ra"),
  dec: $<HTMLElement>("dec"),
  fov: $<HTMLElement>("fov"),
  rot: $<HTMLElement>("rot"),
  scale: $<HTMLElement>("scale"),
  rms: $<HTMLElement>("rms"),
  matches: $<HTMLElement>("matches"),
  elapsed: $<HTMLElement>("elapsed"),
  io: $<HTMLElement>("io"),
  toggleStars: $<HTMLInputElement>("toggleStars"),
  toggleConstel: $<HTMLInputElement>("toggleConstel"),
  toggleDso: $<HTMLInputElement>("toggleDso"),
};

let wasmBytes: ArrayBuffer | null = null;
let gray: GrayImage | null = null;
let result: SolveResult | null = null;
let worker: Worker | null = null;

function setEngine(state: string, text: string): void {
  els.engine.dataset.state = state;
  els.engine.textContent = text;
}

function setStage(stage: string, pct: number): void {
  els.progressBlock.hidden = false;
  els.stageText.textContent = stage;
  els.fill.style.width = `${Math.round(pct)}%`;
}

async function bootEngine(): Promise<void> {
  setEngine("working", "加载 WASM 引擎…");
  try {
    const res = await fetch(`${import.meta.env.BASE_URL}assets/platersol_wasm.wasm`);
    if (!res.ok) throw new Error(`WASM fetch HTTP ${res.status}`);
    wasmBytes = await res.arrayBuffer();
    setEngine("ready", "引擎就绪 · 等待照片");
  } catch (err) {
    setEngine("error", `引擎加载失败: ${msg(err)}`);
  }
}

function msg(err: unknown): string {
  return err instanceof Error ? err.message : String(err);
}

function drawPhoto(img: GrayImage): void {
  const canvas = els.photo;
  canvas.width = img.width;
  canvas.height = img.height;
  const ctx = canvas.getContext("2d")!;
  const out = ctx.createImageData(img.width, img.height);
  // Normalize to 0..255 via robust percentiles computed on a sample.
  const sample: number[] = [];
  const step = Math.max(1, Math.floor(img.pixels.length / 200000));
  for (let i = 0; i < img.pixels.length; i += step) sample.push(img.pixels[i]);
  sample.sort((a, b) => a - b);
  const lo = sample[Math.floor(sample.length * 0.02)];
  const hi = sample[Math.floor(sample.length * 0.995)];
  const span = Math.max(1e-6, hi - lo);
  for (let p = 0, i = 0; p < img.pixels.length; p++, i += 4) {
    const raw = ((img.pixels[p] - lo) / span) * 255;
    const v = Math.round(Math.min(255, Math.max(0, raw)));
    out.data[i] = v;
    out.data[i + 1] = v;
    out.data[i + 2] = v;
    out.data[i + 3] = 255;
  }
  ctx.putImageData(out, 0, 0);
  els.overlay.width = img.width;
  els.overlay.height = img.height;
  els.wrap.classList.add("has-image");
}

els.fileInput.addEventListener("change", async () => {
  const file = els.fileInput.files?.[0];
  if (!file) return;
  setEngine("working", "解码图像…");
  try {
    gray = await decodeToGray(file);
    drawPhoto(gray);
    els.solveBtn.disabled = wasmBytes === null;
    setEngine("ready", `已载入 ${gray.width}×${gray.height} · 可解算`);
  } catch (err) {
    setEngine("error", `图像解码失败: ${msg(err)}`);
  }
});

// Drag and drop.
["dragover", "drop"].forEach((evt) =>
  els.wrap.addEventListener(evt, (e) => e.preventDefault()),
);
els.wrap.addEventListener("drop", async (e: DragEvent) => {
  const file = e.dataTransfer?.files?.[0];
  if (!file) return;
  els.fileInput.files = e.dataTransfer!.files;
  els.fileInput.dispatchEvent(new Event("change"));
});

els.solveBtn.addEventListener("click", () => {
  if (!gray || !wasmBytes) return;
  if (worker) worker.terminate();
  worker = new SolverWorker();
  setEngine("working", "解算中…");
  els.resultBlock.hidden = true;
  setStage("启动…", 5);

  worker.onmessage = (ev: MessageEvent) => handleWorkerMessage(ev.data);
  worker.onerror = (e) => {
    setEngine("error", `Worker 错误: ${e.message}`);
    setStage("失败", 100);
  };

  // Transfer the pixel buffer (zero-copy) to the worker.
  const pixelCopy = new Float32Array(gray.pixels);
  worker.postMessage(
    {
      type: "solve",
      wasmBytes,
      indexBaseUrl: els.indexUrl.value.trim(),
      gray: { width: gray.width, height: gray.height, pixels: pixelCopy.buffer },
      thresholdSigma: parseFloat(els.threshold.value),
      hashTol: parseInt(els.hashTol.value, 10),
      quorum: parseInt(els.quorum.value, 10),
      brightest: 40,
      maxBucketReads: 600,
    },
    [pixelCopy.buffer],
  );
});

function handleWorkerMessage(m: {
  type: string;
  stage?: string;
  detail?: string;
  sources?: SourceC[];
  result?: SolveResult;
  matches?: Array<[number, number]>;
  catalog?: CatStarRow[];
  error?: string;
}): void {
  switch (m.type) {
    case "status":
      setStage(`${m.stage ?? ""} · ${m.detail ?? ""}`, stagePct(m.stage ?? ""));
      break;
    case "sources":
      setStage(`星点检测完成：${m.sources?.length ?? 0} 个源`, 45);
      break;
    case "solved":
      if (!m.result || !gray) return;
      result = m.result;
      presentResult(m.result);
      drawAnnotations();
      setStage("解算完成", 100);
      setEngine("solved", "解算成功");
      break;
    case "failed":
      setEngine("error", `解算失败: ${m.error}`);
      setStage(`失败: ${m.error}`, 100);
      break;
  }
}

function stagePct(stage: string): number {
  const map: Record<string, number> = {
    wasm: 10, index: 20, image: 25, extract: 40, plan: 55, chunks: 75, solve: 90,
  };
  return map[stage] ?? 30;
}

function presentResult(r: SolveResult): void {
  els.resultBlock.hidden = false;
  const raHms = toHms(r.ra0Deg / 15);
  const decDms = toDms(r.dec0Deg);
  els.ra.textContent = `${r.ra0Deg.toFixed(5)}° (${raHms})`;
  els.dec.textContent = `${r.dec0Deg.toFixed(5)}° (${decDms})`;
  els.fov.textContent = `${r.fovXDeg.toFixed(3)}° × ${r.fovYDeg.toFixed(3)}°`;
  els.rot.textContent = `${r.rotationDeg.toFixed(3)}°`;
  els.scale.textContent = `${r.plateScaleArcsec.toFixed(3)} 角秒/像素`;
  els.rms.textContent = `${r.rmsArcsec.toFixed(3)} 角秒`;
  els.matches.textContent = `${r.nInliers} 内点 / ${r.nMatches} 匹配`;
  els.elapsed.textContent = `${r.solveMs.toFixed(0)} ms`;
  els.io.textContent = `${r.bucketsRead} 桶 · ${(Number(r.indexBytes) / 1024).toFixed(0)} KiB`;
}

function drawAnnotations(): void {
  if (!result || !gray) return;
  const ctx = els.overlay.getContext("2d")!;
  const projector = makeProjector(result, gray.width, gray.height);
  const opts: OverlayOptions = {
    showStars: els.toggleStars.checked,
    showConstellations: els.toggleConstel.checked,
    showDso: els.toggleDso.checked,
  };
  renderOverlay(ctx, projector, gray.width, gray.height, opts);
}

[els.toggleStars, els.toggleConstel, els.toggleDso].forEach((t) =>
  t.addEventListener("change", drawAnnotations),
);

function toHms(hours: number): string {
  const h = Math.floor(hours);
  const mF = (hours - h) * 60;
  const m = Math.floor(mF);
  const s = (mF - m) * 60;
  return `${h}h ${m}m ${s.toFixed(1)}s`;
}
function toDms(deg: number): string {
  const sign = deg < 0 ? "-" : "+";
  const a = Math.abs(deg);
  const d = Math.floor(a);
  const mF = (a - d) * 60;
  const m = Math.floor(mF);
  const s = (mF - m) * 60;
  return `${sign}${d}° ${m}′ ${s.toFixed(1)}″`;
}

void bootEngine();
