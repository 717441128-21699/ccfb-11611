//! # platersol-core
//!
//! Offline blind plate solving in pure Rust:
//!
//! - [`image`] — grayscale image + source extraction (background estimation,
//!   connected components, flux-weighted sub-pixel centroids, second moments)
//! - [`quad`] — scale/rotation/reflection-invariant four-star geometric hash
//! - [`catalog`] — star catalog with a static 3-D KD-tree
//! - [`index`] — scale-binned, bucket-chunked quad hash index on disk
//! - [`wcs`] — tangent-plane WCS fitting with optional radial distortion
//! - [`solver`] — the blind solve pipeline over an abstract [`IndexStore`]
//! - [`synthetic`] — synthetic sky/imager for end-to-end accuracy tests

pub mod astro;
pub mod catalog;
pub mod image;
pub mod index;
pub mod quad;
pub mod solver;
pub mod stats;
pub mod synthetic;
pub mod wcs;

pub use astro::{TanProjection, Vec3, ARCSEC, DEG};
pub use catalog::{CatStar, CatStarRow, KdTree};
pub use image::{extract_sources, Background, ExtractConfig, GrayImage, Source};
pub use index::{
    build_index, parse_manifest, BuildConfig, BuiltIndex, IndexManifest, LoadedBin,
};
pub use quad::{hash_distance, make_quad, Pt, Quad};
pub use solver::{
    solve, DirStore, IndexStore, SolveConfig, SolveError, SolveOutcome,
};
pub use synthetic::{render_synthetic_image, SyntheticCamera, SyntheticSky};
pub use wcs::{fit_wcs, Correspondence, FitConfig as WcsFitConfig, WcsSolution};
