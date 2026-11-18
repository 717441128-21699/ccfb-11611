// Background solver worker. Owns the WASM session and ChunkStore, and drives
// the fetch-bucket → feed-bucket → solve loop off the main thread so the UI
// stays responsive during a multi-second blind solve.

import { PlateSolverWasm, type SolveResult, type SourceC } from "../wasm/loader";
import { ChunkStore, type CatStarRow } from "../wasm/chunk-store";
export interface SolveRequest {
  type: "solve";
  wasmBytes: ArrayBuffer;
  indexBaseUrl: string;
  gray: { width: number; height: number; pixels: ArrayBuffer };
  thresholdSigma: number;
  hashTol: number;
  quorum: number;
  brightest: number;
  maxBucketReads: number;
}

export type WorkerOutbound =
  | { type: "status"; stage: string; detail?: string }
  | { type: "sources"; sources: SourceC[] }
  | {
      type: "solved";
      result: SolveResult;
      matches: Array<[number, number]>;
      catalog: CatStarRow[];
    }
  | { type: "failed"; error: string };

const ctx = self as unknown as DedicatedWorkerGlobalScope;

ctx.onmessage = async (ev: MessageEvent<SolveRequest>) => {
  const req = ev.data;
  if (req.type !== "solve") return;
  const post = (m: WorkerOutbound) => ctx.postMessage(m);

  try {
    post({ type: "status", stage: "wasm", detail: "instantiating solver" });
    const wasm = await PlateSolverWasm.instantiate(new Uint8Array(req.wasmBytes));

    post({ type: "status", stage: "index", detail: "opening manifest + catalog" });
    const store = await ChunkStore.open(req.indexBaseUrl);
    wasm.init(JSON.stringify(store.manifest), JSON.stringify(store.catalog));

    post({ type: "status", stage: "image", detail: "loading grayscale pixels" });
    const pixels = new Float32Array(req.gray.pixels);
    wasm.setImageF32(req.gray.width, req.gray.height, pixels);

    post({ type: "status", stage: "extract", detail: "finding star sources" });
    const nSrc = wasm.extract({
      thresholdSigma: req.thresholdSigma,
      minPixels: 3,
      maxPixels: 400,
      border: 4,
      maxSources: 500,
    });
    const sources = wasm.getSources(nSrc);
    post({ type: "sources", sources });
    if (nSrc < 4) throw new Error(`only ${nSrc} sources detected; need at least 4`);

    // Planning pass: the engine reports exactly which (bin, bucket) chunks
    // its quad probes will touch, so only those are fetched over the network.
    post({ type: "status", stage: "plan", detail: "computing required index chunks" });
    const solveCfg = {
      hashTol: req.hashTol,
      quorum: req.quorum,
      brightest: req.brightest,
      maxBucketReads: req.maxBucketReads,
    };
    const needed = wasm.plan(solveCfg);
    post({
      type: "status",
      stage: "chunks",
      detail: `fetching ${needed.length} index buckets`,
    });

    // Fetch in bounded parallelism; empty buckets resolve to zero bytes and
    // are skipped when feeding.
    const CONCURRENCY = 8;
    for (let i = 0; i < needed.length; i += CONCURRENCY) {
      const batch = needed.slice(i, i + CONCURRENCY);
      const chunks = await Promise.all(
        batch.map((r) => store.loadBucket(r.bin, r.bucket).catch(() => new Uint8Array(0))),
      );
      for (let k = 0; k < batch.length; k++) {
        const bytes = chunks[k];
        if (bytes.byteLength > 0) wasm.feedBucket(batch[k].bin, batch[k].bucket, bytes);
      }
    }

    post({ type: "status", stage: "solve", detail: "matching quads and fitting WCS" });
    const result = wasm.solve(solveCfg);

    const matches = wasm.getMatches(result.nMatches + 8);
    post({
      type: "solved",
      result,
      matches,
      catalog: store.catalog,
    });
  } catch (err) {
    post({ type: "failed", error: err instanceof Error ? err.message : String(err) });
  }
};
