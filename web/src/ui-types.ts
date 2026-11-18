/** UI 层本地类型别名（避免直接在 main 里到处引 wasm 包类型）。 */
export type { WasmIndex } from "../pkg/star_solver_wasm.js";

/** 结构兼容的别名，便于测试与解耦。 */
export type WasmIndexLike = {
  star_count: number;
  band_count: number;
};

export interface StarPoint {
  x: number;
  y: number;
  flux: number;
  area: number;
}

export type { PlottedDso } from "./types";
