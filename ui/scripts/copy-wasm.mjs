// Copy the freshly built WASM module into Vite's static asset tree so the UI
// can fetch it as a local file (fully offline; no CDN at runtime).
import { copyFileSync, mkdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const root = resolve(here, "..");
const src = resolve(root, "../target/wasm32-unknown-unknown/release/platersol_wasm.wasm");
const destDir = resolve(root, "public/assets");
const dest = resolve(destDir, "platersol_wasm.wasm");

mkdirSync(destDir, { recursive: true });
copyFileSync(src, dest);
console.log(`copied ${src} -> ${dest}`);
