/**
 * 内置演示星表：确定性 LCG 生成全天均匀星，导出为建库所需 CSV。
 *
 * 仅用于零数据体验与离线 demo；生产环境应导入真实星表（Gaia/Hipparcos 派生）。
 * 这里的随机算法与 Rust 端 synth/selftest 独立（JS 侧无需逐位一致，
 * 只要分布均匀、可复现即可）。
 */

function mulberry32(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a |= 0;
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export interface DemoCatalogOptions {
  count: number;
  magBright: number;
  magFaint: number;
  seed?: number;
}

/** 生成 CSV 行（含表头），列：id,ra_deg,dec_deg,mag。 */
export function generateDemoCsv(opts: DemoCatalogOptions): string {
  const rand = mulberry32(opts.seed ?? 20240517);
  const lines: string[] = ["id,ra_deg,dec_deg,mag"];
  const span = opts.magFaint - opts.magBright;
  for (let id = 0; id < opts.count; id++) {
    const u = rand();
    const v = rand();
    const raDeg = 360 * u;
    // 均匀球面：dec = asin(2v-1)
    const decDeg = (Math.asin(2 * v - 1) * 180) / Math.PI;
    const mag = opts.magBright + span * rand();
    lines.push(
      `${id},${raDeg.toFixed(6)},${decDeg.toFixed(6)},${mag.toFixed(3)}`,
    );
  }
  return lines.join("\n");
}

/** 默认演示规模：3 万颗，建库在浏览器内约几十秒（release wasm）。 */
export const DEFAULT_DEMO: DemoCatalogOptions = {
  count: 30000,
  magBright: 5,
  magFaint: 9.5,
};
