import { defineConfig } from "vite";

// The WASM module and the generated data chunks are served as static assets.
export default defineConfig({
  base: "./",
  build: {
    target: "es2022",
    outDir: "dist",
    assetsInlineLimit: 0,
  },
  worker: {
    format: "es",
  },
  server: {
    headers: {
      // SharedArrayBuffer-free; cross-origin isolation not required since the
      // solver runs on the main module's linear memory inside the worker.
      "Cross-Origin-Opener-Policy": "same-origin",
    },
  },
});
