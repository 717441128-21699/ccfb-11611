// Reference annotation data for overlay rendering after a successful solve.
// Positions are J2000 RA/Dec in degrees. This is a compact built-in set so the
// app works fully offline; a larger named-star/DSO catalog can be loaded from
// the same chunked data directory without changing the rendering path.

export interface NamedStar {
  ra: number;
  dec: number;
  mag: number;
  name: string;
}

export interface Constellation {
  name: string;
  // Pairs of RA/Dec points to connect (great-circle chords drawn on screen).
  lines: Array<[[number, number], [number, number]]>;
}

export interface DeepSkyObject {
  ra: number;
  dec: number;
  type: "nebula" | "galaxy" | "cluster";
  name: string;
  mag: number;
}

// A small, recognizable subset of the brightest named stars.
export const NAMED_STARS: NamedStar[] = [
  { ra: 310.358, dec: 38.784, mag: 0.03, name: "天津四 Deneb" },
  { ra: 304.316, dec: 40.257, mag: 1.25, name: "织女一 Vega" },
  { ra: 297.696, dec: 38.784, mag: 0.77, name: "河鼓二 Altair" },
  { ra: 83.794, dec: -5.391, mag: 0.5, name: "参宿七 Rigel" },
  { ra: 88.793, dec: 7.407, mag: 0.42, name: "参宿四 Betelgeuse" },
  { ra: 95.988, dec: -1.202, mag: -1.46, name: "天狼 Sirius" },
  { ra: 114.825, dec: -28.972, mag: -0.74, name: "老人 Canopus" },
  { ra: 213.915, dec: 19.182, mag: 0.0, name: "大角 Arcturus" },
  { ra: 152.093, dec: -69.028, mag: 0.46, name: "马腹一 Hadar" },
  { ra: 207.54, dec: -16.716, mag: 0.61, name: "角宿一 Spica" },
  { ra: 346.187, dec: -29.876, mag: 1.16, name: "北落师门 Fomalhaut" },
  { ra: 140.813, dec: -60.685, mag: -0.05, name: "南门二 Rigil Kent" },
  { ra: 243.321, dec: -26.432, mag: 0.91, name: "心宿二 Antares" },
  { ra: 54.908, dec: 49.861, mag: 1.79, name: "天船三 Mirfak" },
  { ra: 37.954, dec: 89.264, mag: 1.98, name: "北极星 Polaris" },
  { ra: 10.127, dec: 56.537, mag: 2.28, name: "王良四 Schedar" },
  { ra: 10.897, dec: 59.15, mag: 2.47, name: "策 Gamma Cas" },
  { ra: 42.296, dec: 49.304, mag: 2.86, name: "天船二 Algol" },
  { ra: 191.908, dec: -17.964, mag: 0.97, name: "五帝座一 Denebola" },
  { ra: 15.866, dec: 35.631, mag: 2.42, name: "天船一 Capulus" },
];

// The Big Dipper (part of Ursa Major) as a recognizable connected asterism.
export const CONSTELLATIONS: Constellation[] = [
  {
    name: "大熊座 北斗七星",
    lines: [
      [[165.932, 61.751], [166.449, 56.383]], // Dubhe-Merak
      [[166.449, 56.383], [178.458, 53.694]], // Merak-Phecda
      [[178.458, 53.694], [183.857, 57.033]], // Phecda-Megrez
      [[183.857, 57.033], [165.932, 61.751]], // Megrez-Dubhe (bowl)
      [[183.857, 57.033], [193.507, 55.96]], // Megrez-Alioth
      [[193.507, 55.96], [200.982, 54.925]], // Alioth-Mizar
      [[200.982, 54.925], [210.956, 49.313]], // Mizar-Alkaid
    ],
  },
  {
    name: "夏季大三角",
    lines: [
      [[310.358, 38.784], [304.316, 40.257]], // Deneb-Vega
      [[304.316, 40.257], [297.696, 38.784]], // Vega-Altair
      [[297.696, 38.784], [310.358, 38.784]], // Altair-Deneb
    ],
  },
  {
    name: "猎户座 Orion",
    lines: [
      [[83.794, -5.391], [88.793, 7.407]], // Rigel-Betelgeuse
      [[88.793, 7.407], [85.19, -1.202]], // Betelgeuse-Bellatrix
      [[83.794, -5.391], [81.283, -6.84]], // Rigel-Saiph
      [[85.19, -1.202], [81.283, -6.84]], // Bellatrix-Saiph (frame)
      [[84.053, -1.943], [86.939, -1.943]], // belt Alnitak-Mintaka
    ],
  },
];

// Well-known deep-sky objects (a representative offline subset).
export const DEEP_SKY: DeepSkyObject[] = [
  { ra: 83.822, dec: -5.391, type: "nebula", name: "猎户座大星云 M42", mag: 4.0 },
  { ra: 83.092, dec: -5.91, type: "nebula", name: "火焰星云 M43", mag: 9.0 },
  { ra: 10.685, dec: 41.269, type: "galaxy", name: "仙女座星系 M31", mag: 3.44 },
  { ra: 10.78, dec: 40.864, type: "galaxy", name: "M32", mag: 8.1 },
  { ra: 10.162, dec: 40.996, type: "galaxy", name: "M110", mag: 8.0 },
  { ra: 245.904, dec: -26.432, type: "cluster", name: "蝴蝶星团 M6", mag: 4.2 },
  { ra: 245.15, dec: -26.27, type: "cluster", name: "托勒密星团 M7", mag: 3.3 },
  { ra: 287.85, dec: 33.029, type: "nebula", name: "环状星云 M57", mag: 8.8 },
  { ra: 299.904, dec: 35.202, type: "nebula", name: "面纱星云 NGC6960", mag: 7.0 },
  { ra: 283.396, dec: 33.133, type: "cluster", name: "M27 哑铃星云", mag: 7.5 },
  { ra: 185.531, dec: 4.306, type: "cluster", name: "室女座星系团 M87", mag: 8.6 },
  { ra: 187.7, dec: 12.39, type: "galaxy", name: "黑眼星系 M64", mag: 8.5 },
];
