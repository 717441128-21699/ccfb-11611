/**
 * 浏览器端图像解码与灰度化。
 *
 * 用 createImageBitmap 走原生快速解码，再绘到离屏 Canvas 取灰度字节。
 * 大图按长边等比缩到 maxDim，兼顾检测速度与亚像素精度。
 */

export interface GrayBuffer {
  width: number;
  height: number;
  /** 行优先，0..255 */
  data: Uint8ClampedArray;
  /** 原始（未缩放）尺寸，供结果标注时换算回原图坐标 */
  originalWidth: number;
  originalHeight: number;
  scale: number;
}

/** 解码任意浏览器支持的图片（JPEG/PNG/WebP/BMP）为灰度缓冲。 */
export async function decodeToGray(
  file: Blob,
  maxDim = 4000,
): Promise<GrayBuffer> {
  const bitmap = await createImageBitmap(file);
  const originalWidth = bitmap.width;
  const originalHeight = bitmap.height;

  const scale = Math.max(1, Math.ceil(Math.max(originalWidth, originalHeight) / maxDim));
  const width = Math.floor(originalWidth / scale);
  const height = Math.floor(originalHeight / scale);

  const canvas = document.createElement("canvas");
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext("2d", { willReadFrequently: true });
  if (!ctx) throw new Error("无法创建 2D 绘图上下文");
  ctx.drawImage(bitmap, 0, 0, width, height);
  bitmap.close();

  const rgba = ctx.getImageData(0, 0, width, height).data;
  // ITU-R BT.601 亮度加权
  const gray = new Uint8ClampedArray(width * height);
  for (let i = 0, j = 0; i < rgba.length; i += 4, j++) {
    gray[j] = Math.round(
      0.299 * rgba[i] + 0.587 * rgba[i + 1] + 0.114 * rgba[i + 2],
    );
  }
  return { width, height, data: gray, originalWidth, originalHeight, scale };
}

/** 把检测/解算坐标（降采样图）换回原图坐标。 */
export function toOriginal(coord: number, scale: number): number {
  return coord * scale;
}
