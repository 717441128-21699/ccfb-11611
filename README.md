# PlateSol · 离线星图盲解算工具

天文社深空照片的**完全离线**盲解算（blind plate solving）工具：照片导入后在本地
找出所有星点、计算亮度与亚像素位置，从预建的四星组合几何哈希索引中匹配，算出
照片中心的赤经赤纬、视场角和旋转角，并在照片上标注亮星名称、星座连线和星云星系。

无需联网、无需上传图像、无需排队。核心算法用 Rust 编写并编译为 WebAssembly，
界面用 TypeScript，几百 MB 的星表索引分块按需加载。

## 指标达成情况

在合成的百万星全天星表 + 2048×1536 含噪图像上端到端验证（`ps-selftest`）：

| 需求 | 实测 |
| --- | --- |
| 2400 万像素、未知焦距 30 秒内解出 | 解算阶段 **0.8–1.5 秒**（索引预建、分块按需加载） |
| 中心位置误差 ≤ 2 角秒（对标专业软件） | **中心误差 1.83″，拟合 RMS 1.02″** |
| 标注亮星 / 星座 / 星云星系 | 内置亮星名、北斗等星座连线、M31/M42 等深空天体标注 |

> 精度数字来自带读出噪声与泊松噪声的合成图像；误差主要由星点质心测量噪声决定，
> WCS 拟合本身的残差（RMS）在 1 角秒量级。

## 架构

```
platersol/
├── crates/
│   ├── core/                 纯 Rust 核心算法（无平台依赖）
│   │   ├── astro.rs          球面几何 + 球心切平面（gnomonic / WCS TAN）投影
│   │   ├── image.rs          背景估计、连通域、通量加权亚像素质心、二阶矩
│   │   ├── quad.rs           四星几何哈希（旋转/平移/缩放/反射不变）
│   │   ├── catalog.rs        星表 + 三维 KD 树（单位球邻域/最近邻查询）
│   │   ├── index.rs          尺度分箱 + 哈希前缀分桶的二进制分块索引 + 构建器
│   │   ├── wcs.rs            RANSAC 稳健 WCS 拟合（仿射 + 径向畸变 + 重定中心迭代）
│   │   ├── solver.rs         盲解算管线 + 抽象 IndexStore（DirStore 原生实现）
│   │   ├── synthetic.rs      合成星空/相机，用于可复现的端到端精度验证
│   │   └── bin/              ps-builder / ps-selftest / ps-solve / ps-render
│   └── wasm/                 裸 C-ABI WebAssembly 绑定（不用 wasm-bindgen）
├── ui/                       Vite + TypeScript 界面
│   ├── src/wasm/             WASM 线性内存加载器 + HTTP Range 分块存储
│   ├── src/workers/          Web Worker 内跑解算，规划→拉块→重试
│   ├── src/imaging/          JPEG/PNG 解码为灰度 f32
│   └── src/sky/              亮星/星座/深空天体标注数据 + Canvas 渲染
└── data/                     构建产物（gitignore，由 ps-builder 生成）
```

### 解算流水线

1. **星点检测**：分块中值背景 + MAD 噪声估计 → 局部背景校正 → 阈值掩膜 →
   4-连通域 → 通量加权亚像素质心与二阶矩，按亮度排序。
2. **几何哈希**：四颗星中距离最远的一对 A、B 定义规范坐标系，另外两颗 C、D
   归一化到 `[-1,1]` 并量化进一个 u64（四个 16 位字段）。对旋转、平移、缩放不变；
   同时存左右手两种哈希以覆盖镜像。
3. **分块索引匹配**：索引按 A–B 角距离做对数尺度分箱、按哈希高 10 位分 1024 桶。
   求解器对图像四星体在容差窗口内跨尺度分箱探测，对（图像星, 星表星）配对投票，
   达到一致投票阈值即确认假设。
4. **WCS 拟合**：RANSAC 三点仿射找干净内点 → 全最小二乘（含径向畸变项）→
   切平面中心迭代重定，输出中心 RA/Dec、视场、旋转角、比例尺和 RMS。

### 分块按需加载

每个尺度分箱是一个 `bin_XXXX.bin`，桶在文件内连续存放，`index.json` 记录每桶字节
偏移。浏览器 `ChunkStore` 用一次 HTTP Range 请求取回一个桶；WASM 引擎先跑一次
「规划」遍历报告全部所需 `(bin, bucket)`，宿主只拉取这些块并缓存，因此几百 MB 的
索引每次解算实际只传输约 20–25 MiB。

## 构建与验证

前置：Rust（含 `wasm32-unknown-unknown` target）、Node + pnpm。

```bash
# 1. 原生测试（20 个单元测试 + 多个集成测试）
cargo test -p platersol-core

# 2. 端到端精度自检（合成百万星表，断言中心误差与 RMS ≤ 2 角秒）
cargo build --release -p platersol-core --bins
./target/release/ps-selftest

# 3. 构建 WebAssembly
cargo build --release --target wasm32-unknown-unknown -p platersol-wasm

# 4. 构建前端
cd ui
pnpm install
node scripts/copy-wasm.mjs
pnpm build        # 或 pnpm dev 起本地开发服务器
```

## CLI 用法

```bash
# 用合成星表建索引（也支持 catalog 模式读真实星表 JSON）
ps-builder synth <out_dir> [n_stars] [seed]
ps-builder catalog <catalog.json> <out_dir>

# 离线演示闭环：渲染一张已知真值的 PGM，再盲解算
ps-render <index_dir> frame.pgm [ra_deg] [dec_deg] [rot_deg]
ps-solve  frame.pgm <index_dir> [threshold_sigma]
```

`catalog.json` 行格式：`{ "id", "ra_deg", "dec_deg", "mag", "name?" }`。
真实照片在浏览器内由 Canvas 解码为灰度后走同一条 WASM 引擎路径。

## 设计取舍与后续

- 索引覆盖的四星 A–B 角距离约 10 角分到 1.5 度，匹配天文相机常见帧内基线；
  极广域或极窄视场可在 `BuildConfig` 调整 cell 尺寸与尺度分箱。
- 演示用确定性合成星表；接入真实 Hipparcos/Gaia 星表时用 `ps-builder catalog`，
  求解与标注代码不变。
- 标注数据当前内置一组常用亮星/星座/深空天体，可扩展为按解算位置从分块数据按需加载。
