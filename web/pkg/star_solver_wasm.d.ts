/* tslint:disable */
/* eslint-disable */
/**
 * 用内存索引对一张灰度图做盲解。
 *
 * `index` 为 [`build_index`] 返回的句柄；`scale_lo/hi` 是先验比例尺范围
 * （角秒/像元），未知焦距时给宽区间（如 0.3~60）。
 */
export function solve_image(index: WasmIndex, gray: Uint8Array, width: number, height: number, scale_lo: number, scale_hi: number): SolveResult;
/**
 * 内置合成星空自检（无需任何输入数据）。
 */
export function run_selftest(): SolveResult;
/**
 * 暴露镜像编码（JS 侧若要自行实现查询可用）。
 */
export function mirror_code_of(code: number): number;
/**
 * 用解出的 WCS 把 DSO 表投影到画幅，返回落在图内天体的 JSON 数组字符串。
 *
 * 入参 `wcs_json` 取自 solve_image 结果里的 `wcs` 字段；`dso_csv` 为
 * `name,ra_deg,dec_deg,kind,mag,major_arcmin,minor_arcmin,pa_deg` 格式。
 * 输出元素形如 `{"name","kind","x","y","ra_deg","dec_deg","mag"}`。
 */
export function project_dso(wcs_json: string, dso_csv: string, width: number, height: number): string;
/**
 * 星点检测入口：灰度字节（行优先，0..255）+ 宽高。
 */
export function detect(gray: Uint8Array, width: number, height: number, k_sigma: number, max_stars: number): DetectedStar[];
/**
 * 背景估计（供 UI 显示与手动调阈值）。
 */
export function background(gray: Uint8Array, width: number, height: number): BackgroundStats;
/**
 * 由 CSV 星表文本在内存中构建索引（小星表/演示用，返回各文件的字节表）。
 */
export function build_index(csv_text: string, mag_limit: number): WasmIndex;
export class BackgroundStats {
  private constructor();
  free(): void;
  level: number;
  sigma: number;
}
/**
 * 检测到的星点（JS 友好的扁平结构）。
 */
export class DetectedStar {
  private constructor();
  free(): void;
  x: number;
  y: number;
  flux: number;
  area: number;
}
/**
 * 完整解算的 JS 结果。
 */
export class SolveResult {
  private constructor();
  free(): void;
  readonly ok: boolean;
  /**
   * 成功时的完整 JSON（中心、比例尺、滚动、视场、内点、FITS 头、匹配星）。
   */
  readonly json: string;
  readonly error: string;
}
/**
 * 构建产物的持有者（不透明句柄，JS 不解码内部布局）。
 */
export class WasmIndex {
  private constructor();
  free(): void;
  /**
   * 第 band 个尺度带文件字节。
   */
  band_bytes(band: number): Uint8Array | undefined;
  /**
   * 合并后的星表块字节（含块头，可直接落盘为 stars.bin）。
   */
  star_bytes(): Uint8Array;
  /**
   * manifest 的 JSON 文本（供 UI 展示/下载）。
   */
  manifest_json(): string;
  readonly band_count: number;
  readonly star_count: number;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
  readonly memory: WebAssembly.Memory;
  readonly __wbg_backgroundstats_free: (a: number, b: number) => void;
  readonly __wbg_detectedstar_free: (a: number, b: number) => void;
  readonly __wbg_get_backgroundstats_level: (a: number) => number;
  readonly __wbg_get_backgroundstats_sigma: (a: number) => number;
  readonly __wbg_get_detectedstar_area: (a: number) => number;
  readonly __wbg_get_detectedstar_flux: (a: number) => number;
  readonly __wbg_set_backgroundstats_level: (a: number, b: number) => void;
  readonly __wbg_set_backgroundstats_sigma: (a: number, b: number) => void;
  readonly __wbg_set_detectedstar_area: (a: number, b: number) => void;
  readonly __wbg_set_detectedstar_flux: (a: number, b: number) => void;
  readonly __wbg_solveresult_free: (a: number, b: number) => void;
  readonly __wbg_wasmindex_free: (a: number, b: number) => void;
  readonly background: (a: number, b: number, c: number, d: number) => [number, number, number];
  readonly build_index: (a: number, b: number, c: number) => [number, number, number];
  readonly detect: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number, number, number];
  readonly mirror_code_of: (a: number) => number;
  readonly project_dso: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number, number, number];
  readonly run_selftest: () => number;
  readonly solve_image: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => number;
  readonly solveresult_error: (a: number) => [number, number];
  readonly solveresult_json: (a: number) => [number, number];
  readonly solveresult_ok: (a: number) => number;
  readonly wasmindex_band_bytes: (a: number, b: number) => [number, number];
  readonly wasmindex_band_count: (a: number) => number;
  readonly wasmindex_manifest_json: (a: number) => [number, number];
  readonly wasmindex_star_bytes: (a: number) => [number, number];
  readonly wasmindex_star_count: (a: number) => number;
  readonly __wbg_get_detectedstar_x: (a: number) => number;
  readonly __wbg_get_detectedstar_y: (a: number) => number;
  readonly __wbg_set_detectedstar_x: (a: number, b: number) => void;
  readonly __wbg_set_detectedstar_y: (a: number, b: number) => void;
  readonly __wbindgen_export_0: WebAssembly.Table;
  readonly __wbindgen_malloc: (a: number, b: number) => number;
  readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
  readonly __externref_table_dealloc: (a: number) => void;
  readonly __wbindgen_free: (a: number, b: number, c: number) => void;
  readonly __externref_drop_slice: (a: number, b: number) => void;
  readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;
/**
* Instantiates the given `module`, which can either be bytes or
* a precompiled `WebAssembly.Module`.
*
* @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
*
* @returns {InitOutput}
*/
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
* If `module_or_path` is {RequestInfo} or {URL}, makes a request and
* for everything else, calls `WebAssembly.instantiate` directly.
*
* @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
*
* @returns {Promise<InitOutput>}
*/
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
