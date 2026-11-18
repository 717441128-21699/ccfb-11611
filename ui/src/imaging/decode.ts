// Decode a browser image file into a grayscale f32 buffer suitable for the
// WASM source extractor. Large photos (up to 24 MP and beyond) are decoded at
// full resolution onto an OffscreenCanvas; an optional downscale keeps the
// longest side bounded when memory is constrained.

export interface GrayImage {
  width: number;
  height: number;
  pixels: Float32Array; // row-major, higher = brighter
}

export interface DecodeOptions {
  /** Cap the longest side; 0/undefined keeps native resolution. */
  maxLongSide?: number;
}

export async function decodeToGray(file: Blob, opts: DecodeOptions = {}): Promise<GrayImage> {
  const bitmap = await createImageBitmap(file);
  let { width, height } = bitmap;

  const cap = opts.maxLongSide ?? 0;
  let scale = 1;
  if (cap > 0) {
    const longSide = Math.max(width, height);
    if (longSide > cap) scale = cap / longSide;
  }
  if (scale !== 1) {
    width = Math.max(1, Math.round(width * scale));
    height = Math.max(1, Math.round(height * scale));
  }

  const canvas = new OffscreenCanvas(width, height);
  const ctx = canvas.getContext("2d", { willReadFrequently: true });
  if (!ctx) throw new Error("2d canvas context unavailable");
  ctx.drawImage(bitmap, 0, 0, width, height);
  bitmap.close();

  const imageData = ctx.getImageData(0, 0, width, height);
  const rgba = imageData.data;
  const pixels = new Float32Array(width * height);

  // ITU-R BT.601 luma, applied to straight (non-premultiplied) canvas pixels.
  for (let i = 0, p = 0; i < rgba.length; i += 4, p++) {
    pixels[p] = 0.299 * rgba[i] + 0.587 * rgba[i + 1] + 0.114 * rgba[i + 2];
  }
  return { width, height, pixels };
}
