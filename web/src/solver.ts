/**
 * WASM 模块的封装与类型映射。
 * 首次调用时惰性初始化，全部计算在本地完成，无任何网络请求。
 */
import init, {
  build_index,
  detect,
  background,
  solve_image,
  run_selftest,
  project_dso,
  type WasmIndex,
  type DetectedStar,
  type BackgroundStats,
  type SolveResult,
} from "../pkg/star_solver_wasm.js";
import wasmUrl from "../pkg/star_solver_wasm_bg.wasm?url";
import type { PlottedDso, SolveJson } from "./types";

export type { SolveJson } from "./types";

let ready: Promise<void> | null = null;

/** 确保 wasm 已实例化（幂等）。 */
export function ensureWasm(): Promise<void> {
  if (!ready) {
    ready = init(wasmUrl).then(() => undefined);
  }
  return ready;
}

export interface StarPoint {
  x: number;
  y: number;
  flux: number;
  area: number;
}

export interface BgStats {
  level: number;
  sigma: number;
}

export async function detectStars(
  gray: Uint8Array,
  width: number,
  height: number,
  kSigma = 5,
  maxStars = 500,
): Promise<StarPoint[]> {
  await ensureWasm();
  const out: DetectedStar[] = detect(gray, width, height, kSigma, maxStars);
  const pts = out.map((s) => ({ x: s.x, y: s.y, flux: s.flux, area: s.area }));
  out.forEach((s) => s.free());
  return pts;
}

export async function backgroundStats(
  gray: Uint8Array,
  width: number,
  height: number,
): Promise<BgStats> {
  await ensureWasm();
  const s: BackgroundStats = background(gray, width, height);
  const res = { level: s.level, sigma: s.sigma };
  s.free();
  return res;
}

/** 由 CSV 文本（id,ra_deg,dec_deg,mag）构建内存索引。 */
export async function buildIndex(csv: string, magLimit: number): Promise<WasmIndex> {
  await ensureWasm();
  return build_index(csv, magLimit);
}

export async function solveImage(
  index: WasmIndex,
  gray: Uint8Array,
  width: number,
  height: number,
  scaleLo: number,
  scaleHi: number,
): Promise<{ ok: boolean; error: string; solution?: SolveJson }> {
  await ensureWasm();
  const r: SolveResult = solve_image(index, gray, width, height, scaleLo, scaleHi);
  const ok = r.ok;
  const error = r.error;
  let solution: SolveJson | undefined;
  if (ok) solution = JSON.parse(r.json) as SolveJson;
  r.free();
  return { ok, error, solution };
}

export async function selfTest(): Promise<{ ok: boolean; error: string; report?: unknown }> {
  await ensureWasm();
  const r = run_selftest();
  const ok = r.ok;
  const error = r.error;
  const report = ok ? JSON.parse(r.json) : undefined;
  r.free();
  return { ok, error, report };
}

/** 用解出的 WCS 把 DSO CSV 投影到画幅，返回图内天体。 */
export async function projectDso(
  wcsJson: string,
  dsoCsv: string,
  width: number,
  height: number,
): Promise<PlottedDso[]> {
  await ensureWasm();
  const raw = project_dso(wcsJson, dsoCsv, width, height);
  return JSON.parse(raw) as PlottedDso[];
}
