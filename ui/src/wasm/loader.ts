// Low-level bindings to the raw C-ABI platersol WASM module. No wasm-bindgen:
// we allocate in the module's linear memory and call the `ps_*` exports
// directly, mirroring crates/wasm/src/lib.rs.

/** Mirrors `SolveResultC` (all f64/u32/u64, little-endian, 15 fields). */
export interface SolveResult {
  ra0Deg: number;
  dec0Deg: number;
  fovXDeg: number;
  fovYDeg: number;
  fovDiagDeg: number;
  rotationDeg: number;
  plateScaleArcsec: number;
  rmsArcsec: number;
  nInliers: number;
  nMatches: number;
  nSources: number;
  bucketsRead: number;
  indexBytes: bigint;
  quadsTried: number;
  solveMs: number;
}

/** Mirrors `SourceC` (four f64: x, y, flux, peak = 32 bytes). */
export interface SourceC {
  x: number;
  y: number;
  flux: number;
  peak: number;
}

interface WasmExports {
  memory: WebAssembly.Memory;
  ps_alloc(len: number): number;
  ps_free(ptr: number, len: number): void;
  ps_last_error(buf: number, cap: number): number;
  ps_init(manifestPtr: number, manifestLen: number, catalogPtr: number, catalogLen: number): number;
  ps_feed_bucket(bin: number, bucket: number, ptr: number, len: number): number;
  ps_bucket_count(): number;
  ps_clear_buckets(): number;
  ps_plan(hashTol: number, quorum: number, brightest: number, maxBucketReads: number): number;
  ps_take_missing_buckets(out: number, capPairs: number): number;
  ps_missing_count(): number;
  ps_set_image_f32(width: number, height: number, ptr: number, len: number): number;
  ps_extract(thresholdSigma: number, minPixels: number, maxPixels: number, border: number, maxSources: number): number;
  ps_get_sources(out: number, cap: number): number;
  ps_solve(hashTol: number, quorum: number, brightest: number, maxBucketReads: number, out: number): number;
  ps_pixel_to_radec(x: number, y: number, ra: number, dec: number): number;
  ps_get_matches(out: number, capPairs: number): number;
  ps_nearest_catalog_star(raDeg: number, decDeg: number, searchRadiusDeg: number, outId: number, outSep: number): number;
  ps_wcs_coeffs(a: number, b: number, ar: number, br: number): number;
}

// SolveResultC field layout (bytes): 11×f64 (88) + 4×u32 (16) + u64 (8) = 112.
const SOLVE_RESULT_SIZE = 112;
const SOURCE_C_SIZE = 32;

export class PlateSolverWasm {
  private exports: WasmExports;
  private enc = new TextEncoder();
  private dec = new TextDecoder();

  private constructor(exports: WasmExports) {
    this.exports = exports;
  }

  static async instantiate(bytes: BufferSource): Promise<PlateSolverWasm> {
    const { instance } = await WebAssembly.instantiate(bytes, {});
    return new PlateSolverWasm(instance.exports as unknown as WasmExports);
  }

  private get memory(): WebAssembly.Memory {
    return this.exports.memory;
  }

  /** Copy a byte slice into freshly allocated module memory; returns ptr+len. */
  private writeBytes(data: Uint8Array): { ptr: number; len: number } {
    const len = data.byteLength;
    const ptr = this.exports.ps_alloc(len);
    if (ptr === 0) throw new Error("ps_alloc returned null");
    new Uint8Array(this.memory.buffer, ptr, len).set(data);
    return { ptr, len };
  }

  private free(ptr: number, len: number): void {
    this.exports.ps_free(ptr, len);
  }

  lastError(): string {
    const cap = 512;
    const ptr = this.exports.ps_alloc(cap);
    try {
      const n = this.exports.ps_last_error(ptr, cap);
      if (n < 0) return `<error buffer needed ${-n} bytes>`;
      const bytes = new Uint8Array(this.memory.buffer, ptr, Math.max(0, n));
      return this.dec.decode(bytes);
    } finally {
      this.free(ptr, cap);
    }
  }

  /** Initialize the session with manifest JSON and catalog JSON. */
  init(manifestJson: string, catalogJson: string): void {
    const m = this.writeBytes(this.enc.encode(manifestJson));
    const c = this.writeBytes(this.enc.encode(catalogJson));
    try {
      const rc = this.exports.ps_init(m.ptr, m.len, c.ptr, c.len);
      if (rc < 0) throw new Error(this.lastError());
    } finally {
      this.free(m.ptr, m.len);
      this.free(c.ptr, c.len);
    }
  }

  /** Feed one index bucket's raw bytes (as laid out by the builder). */
  feedBucket(bin: number, bucket: number, bytes: Uint8Array): void {
    const b = this.writeBytes(bytes);
    try {
      const rc = this.exports.ps_feed_bucket(bin, bucket, b.ptr, b.len);
      if (rc < 0) throw new Error(this.lastError());
    } finally {
      this.free(b.ptr, b.len);
    }
  }

  bucketCount(): number {
    return this.exports.ps_bucket_count();
  }

  clearBuckets(): void {
    this.exports.ps_clear_buckets();
  }

  /**
   * Planning pass: populate the engine's missing-bucket set and return the
   * exact (bin, bucket) pairs the host must fetch before a real solve.
   */
  plan(cfg: { hashTol: number; quorum: number; brightest: number; maxBucketReads: number }): Array<{ bin: number; bucket: number }> {
    const need = this.exports.ps_plan(cfg.hashTol, cfg.quorum, cfg.brightest, cfg.maxBucketReads);
    if (need < 0) throw new Error(this.lastError());
    const bytes = need * 8;
    const ptr = this.exports.ps_alloc(bytes);
    try {
      const rc = this.exports.ps_take_missing_buckets(ptr, need);
      if (rc < 0) throw new Error(this.lastError());
      const dv = new DataView(this.memory.buffer, ptr, bytes);
      const out: Array<{ bin: number; bucket: number }> = [];
      for (let i = 0; i < rc; i++) {
        out.push({ bin: dv.getUint32(i * 8, true), bucket: dv.getUint32(i * 8 + 4, true) });
      }
      return out;
    } finally {
      this.free(ptr, bytes);
    }
  }

  missingCount(): number {
    return this.exports.ps_missing_count();
  }

  /** Load a grayscale f32 image (row-major, height rows of width). */
  setImageF32(width: number, height: number, pixels: Float32Array): void {
    const bytes = new Uint8Array(pixels.buffer, pixels.byteOffset, pixels.byteLength);
    const b = this.writeBytes(bytes);
    try {
      const rc = this.exports.ps_set_image_f32(width, height, b.ptr, width * height);
      if (rc < 0) throw new Error(this.lastError());
    } finally {
      this.free(b.ptr, b.len);
    }
  }

  /** Run source extraction; returns the number of detected sources. */
  extract(cfg: {
    thresholdSigma: number;
    minPixels: number;
    maxPixels: number;
    border: number;
    maxSources: number;
  }): number {
    const n = this.exports.ps_extract(
      cfg.thresholdSigma,
      cfg.minPixels,
      cfg.maxPixels,
      cfg.border,
      cfg.maxSources,
    );
    if (n < 0) throw new Error(this.lastError());
    return n;
  }

  /** Read back the extracted sources. */
  getSources(count: number): SourceC[] {
    const bytes = count * SOURCE_C_SIZE;
    const ptr = this.exports.ps_alloc(bytes);
    try {
      const rc = this.exports.ps_get_sources(ptr, count);
      if (rc < 0) throw new Error(this.lastError());
      const dv = new DataView(this.memory.buffer, ptr, bytes);
      const out: SourceC[] = [];
      for (let i = 0; i < rc; i++) {
        const o = i * SOURCE_C_SIZE;
        out.push({ x: dv.getFloat64(o, true), y: dv.getFloat64(o + 8, true), flux: dv.getFloat64(o + 16, true), peak: dv.getFloat64(o + 24, true) });
      }
      return out;
    } finally {
      this.free(ptr, bytes);
    }
  }

  /** Run the blind solve; throws with the engine's error on failure. */
  solve(cfg: {
    hashTol: number;
    quorum: number;
    brightest: number;
    maxBucketReads: number;
  }): SolveResult {
    const ptr = this.exports.ps_alloc(SOLVE_RESULT_SIZE);
    try {
      const rc = this.exports.ps_solve(
        cfg.hashTol,
        cfg.quorum,
        cfg.brightest,
        cfg.maxBucketReads,
        ptr,
      );
      if (rc < 0) throw new Error(this.lastError());
      const dv = new DataView(this.memory.buffer, ptr, SOLVE_RESULT_SIZE);
      let o = 0;
      const f64 = (): number => {
        const v = dv.getFloat64(o, true);
        o += 8;
        return v;
      };
      const u32 = (): number => {
        const v = dv.getUint32(o, true);
        o += 4;
        return v;
      };
      const ra0Deg = f64();
      const dec0Deg = f64();
      const fovXDeg = f64();
      const fovYDeg = f64();
      const fovDiagDeg = f64();
      const rotationDeg = f64();
      const plateScaleArcsec = f64();
      const rmsArcsec = f64();
      const nInliers = u32();
      const nMatches = u32();
      const nSources = u32();
      const bucketsRead = u32();
      const indexBytes = dv.getBigUint64(o, true);
      o += 8;
      const quadsTried = u32();
      const solveMs = f64();
      return {
        ra0Deg, dec0Deg, fovXDeg, fovYDeg, fovDiagDeg, rotationDeg,
        plateScaleArcsec, rmsArcsec, nInliers, nMatches, nSources,
        bucketsRead, indexBytes, quadsTried, solveMs,
      };
    } finally {
      this.free(ptr, SOLVE_RESULT_SIZE);
    }
  }

  pixelToRaDec(x: number, y: number): { raDeg: number; decDeg: number } {
    const ptr = this.exports.ps_alloc(16);
    try {
      const rc = this.exports.ps_pixel_to_radec(x, y, ptr, ptr + 8);
      if (rc < 0) throw new Error(this.lastError());
      const dv = new DataView(this.memory.buffer, ptr, 16);
      return { raDeg: dv.getFloat64(0, true), decDeg: dv.getFloat64(8, true) };
    } finally {
      this.free(ptr, 16);
    }
  }

  getMatches(pairCap: number): Array<[number, number]> {
    const bytes = pairCap * 8;
    const ptr = this.exports.ps_alloc(bytes);
    try {
      const rc = this.exports.ps_get_matches(ptr, pairCap);
      if (rc < 0) throw new Error(this.lastError());
      const dv = new DataView(this.memory.buffer, ptr, bytes);
      const out: Array<[number, number]> = [];
      for (let i = 0; i < rc; i++) {
        out.push([dv.getUint32(i * 8, true), dv.getUint32(i * 8 + 4, true)]);
      }
      return out;
    } finally {
      this.free(ptr, bytes);
    }
  }

  nearestCatalogStar(raDeg: number, decDeg: number, searchRadiusDeg: number): { id: number; sepArcsec: number } | null {
    const ptr = this.exports.ps_alloc(16);
    try {
      const rc = this.exports.ps_nearest_catalog_star(raDeg, decDeg, searchRadiusDeg, ptr, ptr + 8);
      if (rc < 0) return null;
      const dv = new DataView(this.memory.buffer, ptr, 16);
      return { id: dv.getUint32(0, true), sepArcsec: dv.getFloat64(8, true) };
    } finally {
      this.free(ptr, 16);
    }
  }
}
