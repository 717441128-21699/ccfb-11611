# 离线星图解算工具（Offline Plate Solver）

不联网、不排队的天文照片板解（plate solving）工具。照片导入后在本地完成
星点检测、几何哈希盲解，输出照片中心的赤经赤纬、视场角、旋转角，并在图上
标注亮星与星座连线。

- **核心算法用 Rust 编写，编译为 WebAssembly**，浏览器内本地运行，无任何网络请求；
- **界面用 TypeScript + Vite**，纯前端静态站点，可整体离线部署；
- **星表索引分块按需加载**，几百 MB 索引只读取与当前视场相关的少数块；
- 核心库零外部 crate 依赖（纯 std），可审计、易移植，CLI 与 wasm 共用同一套代码。

## 目录结构

```
.
├── crates/
│   ├── star-solver-core/   # 纯 std 核心：检测 / quad 哈希 / 分块索引 / 匹配 / TAN+SIP WCS
│   ├── star-solver-cli/    # 本地命令行：建库、解图、合成自检（图像解码在此层）
│   └── star-solver-wasm/   # wasm-bindgen 绑定，暴露给 JS
├── web/                    # TS 前端（Vite），含生成的 pkg/ 绑定
└── Cargo.toml              # workspace
```

## 算法流水线

1. **星点检测**（`image.rs`）：整图中位数 + 归一化 MAD 估计背景与噪声，
   自适应阈值分割，并查集做连通域标记，流量加权矩算亚像素质心，
   再用面积、圆度（长短轴比）过滤热点与拖线。
2. **几何哈希**（`quad.rs`）：取最亮的若干星为锚点构造四星组合（quad），
   以相距最远的 A-B 为锚点归一化到局部坐标系，另外两颗 C/D 的坐标对平移、
   旋转、缩放完全不变；四个分量各量化到 1024 桶，拼成 40 bit 哈希。
   镜像通过 `mirror_code` 单独查询兼容。
3. **分块索引**（`index.rs` / `builder.rs`）：星表预计算的 quad 按哈希排序，
   每 256 条切一个 chunk，文件头带每块首码与偏移表，查询用二分定位 +
   4 维 box 范围扫描，配合 LRU 只加载相交块。验证星表按赤纬条带切块。
4. **盲解匹配**（`matcher.rs`）：未知焦距时按可行比例尺反推候选尺度带，
   多带尝试；每个图像 quad 查哈希（含镜像），候选先做 4 点门限，再取邻星
   扩充对应点，鲁棒最小二乘 + 2.5σ MAD 剔除，内点数即置信度。
5. **WCS 拟合**（`wcs.rs`）：TAN（gnomonic）切平面投影，切点迭代收敛后
   解线性最小二乘得到 CRVAL/CRPIX/CD，再可选二阶 SIP 畸变修正；
   输出标准 FITS WCS 头卡片。

## 精度与性能说明

- 合成星空端到端测试（`matcher::tests::end_to_end_synthetic_solve`）覆盖
  检测→建库→盲解→WCS 全链路，中心指向、比例尺、滚动角均有断言。
- 亚像素质心 + TAN/SIP 拟合在参考星密集、质心噪声约 0.05 像元时，
  角残差可做到亚角秒；**与专业软件 ≤2 角秒的对标精度依赖足够密且精确的
  参考星表（如 Gaia DR3 派生）**。内置演示星表仅为流程验证，不代表精度上限。
- 30 秒解 2400 万像素的目标通过：检测降采样、只对最亮 ~40 锚点建 quad、
  索引分块按需读取、命中强解即提前终止共同保证；release 构建的关键路径为
  O(星点) 检测 + O(quad × 候选) 的有界搜索。

## 构建

需要 Rust（含 `wasm32-unknown-unknown` target）与 Node.js。

```bash
# 1. 核心 + CLI（原生）
cargo build --release -p star-solver-cli

# 2. wasm（需要与 wasm-bindgen =0.2.100 匹配的 CLI）
cargo build --release -p star-solver-wasm --target wasm32-unknown-unknown
wasm-bindgen target/wasm32-unknown-unknown/release/star_solver_wasm.wasm \
  --target web --out-dir web/pkg

# 3. 前端
cd web && npm install && npm run build
```

## 命令行用法

```bash
# 用 CSV 星表（id,ra_deg,dec_deg,mag）建索引
star-solver build catalog.csv ./index-dir

# 离线解一张照片，打印中心指向/比例尺/视场/旋转角与 FITS WCS 头
star-solver solve ./index-dir photo.jpg

# 无需任何外部数据的端到端自检
star-solver selftest
```

## Web 用法

`cd web && npm run dev` 后在浏览器打开：先「生成内置演示星表并建库」
（或导入自己的 CSV），再选择照片、设置比例尺先验范围，点「开始解算」。
全程计算在本地 WASM 内完成，断网可用。

## 接入真实星表与星名/星座

- 把 Gaia/Hipparcos 导出为 `id,ra_deg,dec_deg,mag` 的 CSV 即可用
  `star-solver build` 或前端「导入星表 CSV」建库。
- 亮星名称、星座连线、星云星系标注需要额外的静态数据表（星名表、
  IAU 星座线、NGC/IC 天体表）。前端 `annotate.ts` 已预留 `nameOf` 钩子与
  覆盖层绘制，接入数据表后即可渲染；这部分不影响核心解算。
