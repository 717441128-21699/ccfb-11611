import { defineConfig } from "vite";

// wasm 以 target=web 构建，Vite 用 ?url 或默认资源管线均可处理。
export default defineConfig({
  base: "./",
  server: {
    // 大文件上传与本地资源，关掉 fs 限制无关项；保持默认端口。
    port: 5173,
  },
  build: {
    target: "es2020",
    // wasm 由 Vite 作为资源 emit，不内联（体积超阈值）。
    assetsInlineLimit: 0,
  },
  optimizeDeps: {
    exclude: ["star-solver-wasm"],
  },
});
