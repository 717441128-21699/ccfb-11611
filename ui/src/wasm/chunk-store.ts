// Chunked geometric-hash index access. The builder writes one `bin_XXXX.bin`
// per scale bin and an `index.json` manifest with per-bucket byte offsets.
// This store fetches exactly the bucket runs a solve needs via HTTP Range
// requests and caches them, so only a fraction of the multi-hundred-MB index
// is ever transferred or held in memory.

export interface ScaleBin {
  index: number;
  lo_rad: number;
  hi_rad: number;
}

export interface IndexManifest {
  magic: number;
  version: number;
  catalog: string;
  n_stars: number;
  n_quads: number;
  bins: ScaleBin[];
  /** bucket_offsets[bin][bucket] = byte offset; length N_BUCKETS+1. */
  bucket_offsets: number[][];
  build_tol_frac: number;
}

export interface CatStarRow {
  id: number;
  ra_deg: number;
  dec_deg: number;
  mag: number;
  name?: string | null;
}

export const N_BUCKETS = 1024;
const HEADER_SIZE = 36;

export interface FetchRange {
  (url: string, start: number, endExclusive: number): Promise<Uint8Array>;
}

/** Default Range fetch using the Fetch API with an explicit byte range. */
export const fetchRange: FetchRange = async (url, start, end) => {
  const res = await fetch(url, {
    headers: { Range: `bytes=${start}-${end - 1}` },
  });
  if (!res.ok && res.status !== 206) {
    throw new Error(`range fetch failed: HTTP ${res.status} for ${url}`);
  }
  return new Uint8Array(await res.arrayBuffer());
};

export class ChunkStore {
  readonly manifest: IndexManifest;
  readonly catalog: CatStarRow[];
  private baseUrl: string;
  private doFetch: FetchRange;
  private cache = new Map<string, Uint8Array>();
  bytesFetched = 0;
  reads = 0;

  private constructor(
    manifest: IndexManifest,
    catalog: CatStarRow[],
    baseUrl: string,
    doFetch: FetchRange,
  ) {
    this.manifest = manifest;
    this.catalog = catalog;
    this.baseUrl = baseUrl;
    this.doFetch = doFetch;
  }

  static async open(baseUrl: string, doFetch: FetchRange = fetchRange): Promise<ChunkStore> {
    const [manifestRes, catalogRes] = await Promise.all([
      fetch(`${baseUrl}/index.json`),
      fetch(`${baseUrl}/catalog.json`),
    ]);
    if (!manifestRes.ok || !catalogRes.ok) {
      throw new Error(`cannot open index at ${baseUrl}`);
    }
    const manifest = (await manifestRes.json()) as IndexManifest;
    const catalog = (await catalogRes.json()) as CatStarRow[];
    return new ChunkStore(manifest, catalog, baseUrl, doFetch);
  }

  private binUrl(bin: number): string {
    return `${this.baseUrl}/bin_${bin.toString().padStart(4, "0")}.bin`;
  }

  /** Raw entry bytes for one (bin, bucket), fetched on demand and cached. */
  async loadBucket(bin: number, bucket: number): Promise<Uint8Array> {
    const key = `${bin}/${bucket}`;
    const hit = this.cache.get(key);
    if (hit) return hit;

    const offsets = this.manifest.bucket_offsets[bin];
    const start = offsets[bucket];
    const end = offsets[bucket + 1];
    if (start === end) {
      const empty = new Uint8Array(0);
      this.cache.set(key, empty);
      return empty;
    }
    // Bucket offsets in the manifest are absolute within the bin file and
    // already include the 36-byte header.
    const bytes = await this.doFetch(this.binUrl(bin), start, end);
    this.reads += 1;
    this.bytesFetched += bytes.byteLength;
    this.cache.set(key, bytes);
    return bytes;
  }

  /** Preload a set of buckets (e.g. the solver's next probe batch). */
  async preload(requests: Array<{ bin: number; bucket: number }>): Promise<void> {
    await Promise.all(
      requests.map((r) => this.loadBucket(r.bin, r.bucket).catch(() => new Uint8Array(0))),
    );
  }

  get catalogById(): Map<number, CatStarRow> {
    const m = new Map<number, CatStarRow>();
    for (const row of this.catalog) m.set(row.id, row);
    return m;
  }

  get cacheMb(): number {
    let total = 0;
    for (const v of this.cache.values()) total += v.byteLength;
    return total / (1024 * 1024);
  }

  get headerSize(): number {
    return HEADER_SIZE;
  }
}
