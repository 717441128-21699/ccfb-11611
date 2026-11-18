//! Raw C-ABI WebAssembly surface for platersol.
//!
//! No `wasm-bindgen`: the TypeScript loader in `ui/src/wasm/loader.ts`
//! allocates linear-memory buffers and calls these exports directly, so the
//! toolchain is just `cargo build --target wasm32-unknown-unknown`.
//!
//! Error convention: functions returning `i32` return 0 on success and a
//! negative error code otherwise. Textual details are copied into the
//! thread-local last-error buffer, readable via `ps_last_error`.

// Raw pointers crossing the WASM/JS boundary are an intentional FFI contract;
// the JS loader is the sole caller and always supplies valid pointers.
#![allow(clippy::not_unsafe_ptr_arg_deref)]
#![allow(clippy::missing_safety_doc)]

use platersol_core::astro::{DEG, Vec3};
use platersol_core::catalog::{CatStar, CatStarRow, KdTree};
use platersol_core::image::{ExtractConfig, GrayImage, Source};
use platersol_core::index::{parse_manifest, IndexManifest};
use platersol_core::solver::{solve, IndexStore, SolveConfig, SolveOutcome};
use platersol_core::wcs::WcsSolution;
use std::cell::RefCell;
use std::collections::HashMap;
use std::sync::atomic::{AtomicU32, Ordering};

thread_local! {
    static ERROR: RefCell<String> = const { RefCell::new(String::new()) };
}

fn set_err(msg: impl Into<String>) -> i32 {
    ERROR.with(|e| *e.borrow_mut() = msg.into());
    -1
}

/// Streaming index store used inside the browser. The JS side owns the
/// chunk fetch/cache logic and pushes bucket bytes in via `ps_feed_bucket`.
pub struct WasmStore {
    manifest: IndexManifest,
    catalog: HashMap<u32, CatStar>,
    buckets: HashMap<(usize, usize), Vec<u8>>,
    bytes_fetched: u64,
    /// Ids the JS side has registered as available (resident catalog).
    known: Vec<u32>,
    /// Buckets the engine requested but JS has not fed yet. A planning pass
    /// populates this so the host fetches exactly the needed chunks.
    pub(crate) missing: std::collections::HashSet<(usize, usize)>,
    /// When true, `load_bucket` records a miss but returns an empty bucket
    /// instead of an error, so a planning solve walks its *full* probe path
    /// and records every chunk a real solve would touch.
    pub(crate) planning: bool,
}

impl IndexStore for WasmStore {
    fn manifest(&self) -> &IndexManifest {
        &self.manifest
    }

    fn load_bucket(&mut self, bin: usize, bucket: usize) -> std::io::Result<Vec<u8>> {
        if let Some(v) = self.buckets.get(&(bin, bucket)) {
            return Ok(v.clone());
        }
        self.missing.insert((bin, bucket));
        if self.planning {
            // Empty bucket: the solver treats this as "no candidates here"
            // and continues probing the remaining bins/buckets.
            return Ok(Vec::new());
        }
        Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("bucket {bin}/{bucket} not fed to solver"),
        ))
    }

    fn star_by_id(&self, id: u32) -> Option<CatStar> {
        self.catalog.get(&id).copied()
    }

    fn known_ids(&self) -> Vec<u32> {
        self.known.clone()
    }

    fn bytes_fetched(&self) -> u64 {
        self.bytes_fetched
    }
}

/// Global solver session (one at a time is plenty for a browser tab).
pub struct Session {
    store: WasmStore,
    image: Option<GrayImage>,
    sources: Vec<Source>,
    outcome: Option<SolveOutcome>,
}

static SESSION_ID: AtomicU32 = AtomicU32::new(1);
thread_local! {
    static SESSION: RefCell<Option<Session>> = const { RefCell::new(None) };
}

fn with_session<R>(f: impl FnOnce(&mut Session) -> R) -> Result<R, i32> {
    SESSION.with(|s| {
        if s.borrow().is_none() {
            return Err(set_err("session not initialized; call ps_init first"));
        }
        Ok(f(s.borrow_mut().as_mut().unwrap()))
    })
}

/// POD result struct mirrored exactly by the TypeScript typings.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct SolveResultC {
    pub ra0_deg: f64,
    pub dec0_deg: f64,
    pub fov_x_deg: f64,
    pub fov_y_deg: f64,
    pub fov_diag_deg: f64,
    pub rotation_deg: f64,
    pub plate_scale_as: f64,
    pub rms_arcsec: f64,
    pub n_inliers: u32,
    pub n_matches: u32,
    pub n_sources: u32,
    pub buckets_read: u32,
    pub index_bytes: u64,
    pub quads_tried: u32,
    pub solve_ms: f64,
}

#[no_mangle]
pub extern "C" fn ps_alloc(len: usize) -> *mut u8 {
    let mut v: Vec<u8> = Vec::with_capacity(len);
    let p = v.as_mut_ptr();
    std::mem::forget(v);
    p
}

#[no_mangle]
pub extern "C" fn ps_free(ptr: *mut u8, len: usize) {
    if !ptr.is_null() && len > 0 {
        unsafe {
            drop(Vec::from_raw_parts(ptr, len, len));
        }
    }
}

#[no_mangle]
pub extern "C" fn ps_last_error(buf: *mut u8, cap: usize) -> i32 {
    ERROR.with(|e| {
        let msg = e.borrow();
        let bytes = msg.as_bytes();
        let n = bytes.len().min(cap.saturating_sub(1));
        if buf.is_null() || cap == 0 {
            return -(bytes.len() as i32) - 1;
        }
        unsafe {
            std::ptr::copy_nonoverlapping(bytes.as_ptr(), buf, n);
            *buf.add(n) = 0;
        }
        n as i32
    })
}

/// Initialize a session with a manifest JSON blob and a catalog JSON blob
/// (array of CatStarRow). Takes ownership of the passed buffers.
#[no_mangle]
pub extern "C" fn ps_init(
    manifest_ptr: *const u8,
    manifest_len: usize,
    catalog_ptr: *const u8,
    catalog_len: usize,
) -> i32 {
    let manifest_bytes = unsafe { std::slice::from_raw_parts(manifest_ptr, manifest_len) };
    let catalog_bytes = unsafe { std::slice::from_raw_parts(catalog_ptr, catalog_len) };
    let manifest = match parse_manifest(manifest_bytes) {
        Ok(m) => m,
        Err(e) => return set_err(format!("bad manifest: {e}")),
    };
    let rows: Vec<CatStarRow> = match serde_json::from_slice(catalog_bytes) {
        Ok(r) => r,
        Err(e) => return set_err(format!("bad catalog: {e}")),
    };
    let catalog: HashMap<u32, CatStar> = rows
        .iter()
        .map(|r| {
            let ra = r.ra_deg * DEG;
            let dec = r.dec_deg * DEG;
            (
                r.id,
                CatStar {
                    id: r.id,
                    ra,
                    dec,
                    mag: r.mag,
                    v: Vec3::from_radec(ra, dec),
                },
            )
        })
        .collect();
    let known = catalog.keys().copied().collect();
    let id = SESSION_ID.fetch_add(1, Ordering::Relaxed);
    SESSION.with(|s| {
        *s.borrow_mut() = Some(Session {
            store: WasmStore {
                manifest,
                catalog,
                buckets: HashMap::new(),
                bytes_fetched: 0,
                known,
                missing: std::collections::HashSet::new(),
                planning: false,
            },
            image: None,
            sources: Vec::new(),
            outcome: None,
        });
    });
    id as i32
}

/// Feed one index bucket (`bin`/`bucket` ids from the manifest). The solver
/// copies the bytes, so the caller may free its buffer afterwards.
#[no_mangle]
pub extern "C" fn ps_feed_bucket(
    bin: u32,
    bucket: u32,
    ptr: *const u8,
    len: usize,
) -> i32 {
    let bytes = unsafe { std::slice::from_raw_parts(ptr, len) }.to_vec();
    SESSION.with(|s| {
        if let Some(sess) = s.borrow_mut().as_mut() {
            sess.store.bytes_fetched += len as u64;
            let key = (bin as usize, bucket as usize);
            sess.store.buckets.insert(key, bytes);
            sess.store.missing.remove(&key);
            0
        } else {
            set_err("ps_feed_bucket before ps_init")
        }
    })
}

/// How many buckets are currently resident (for the UI cache meter).
#[no_mangle]
pub extern "C" fn ps_bucket_count() -> u32 {
    SESSION.with(|s| {
        s.borrow()
            .as_ref()
            .map(|sess| sess.store.buckets.len() as u32)
            .unwrap_or(0)
    })
}

/// Load a grayscale image from an f32 pixel buffer (row-major, normalized or
/// raw flux — extraction only relies on relative values).
#[no_mangle]
pub extern "C" fn ps_set_image_f32(
    width: u32,
    height: u32,
    ptr: *const f32,
    len: usize,
) -> i32 {
    if (width as usize) * (height as usize) != len {
        return set_err("image buffer length does not match width*height");
    }
    let buf = unsafe { std::slice::from_raw_parts(ptr, len) }.to_vec();
    SESSION.with(|s| match s.borrow_mut().as_mut() {
        Some(sess) => {
            sess.image = Some(GrayImage {
                width: width as usize,
                height: height as usize,
                buf,
            });
            sess.sources.clear();
            sess.outcome = None;
            0
        }
        None => set_err("ps_set_image before ps_init"),
    })
}

/// Run source extraction on the loaded image. Returns the number of sources.
#[no_mangle]
pub extern "C" fn ps_extract(
    threshold_sigma: f32,
    min_pixels: u32,
    max_pixels: u32,
    border: u32,
    max_sources: u32,
) -> i32 {
    let cfg = ExtractConfig {
        threshold_sigma,
        min_pixels: min_pixels as usize,
        max_pixels: max_pixels as usize,
        border: border as usize,
        max_sources: if max_sources == 0 {
            None
        } else {
            Some(max_sources as usize)
        },
    };
    SESSION.with(|s| {
        let mut sess = s.borrow_mut();
        let Some(sess) = sess.as_mut() else {
            return set_err("session not initialized");
        };
        let Some(img) = sess.image.as_ref() else {
            return set_err("no image loaded");
        };
        sess.sources = platersol_core::extract_sources(img, &cfg);
        sess.sources.len() as i32
    })
}

/// Copy extracted sources into a caller-provided `SourcesC` array
/// (x, y, flux, peak f64s × n). Returns the required count when `cap` is too
/// small (negative = needed count).
#[repr(C)]
#[derive(Clone, Copy)]
pub struct SourceC {
    pub x: f64,
    pub y: f64,
    pub flux: f64,
    pub peak: f64,
}

#[no_mangle]
pub extern "C" fn ps_get_sources(out: *mut SourceC, cap: u32) -> i32 {
    SESSION.with(|s| {
        let sess = s.borrow();
        let Some(sess) = sess.as_ref() else {
            return set_err("session not initialized");
        };
        let n = sess.sources.len();
        if (cap as usize) < n {
            return -(n as i32);
        }
        for (i, src) in sess.sources.iter().enumerate() {
            unsafe {
                *out.add(i) = SourceC {
                    x: src.x,
                    y: src.y,
                    flux: src.flux,
                    peak: src.peak,
                };
            }
        }
        n as i32
    })
}

/// Run the blind solve. Bucket bytes must already have been fed for every
/// (bin, bucket) the solver may request; use `ps_requested_buckets` after a
/// failed solve to learn which ones to fetch, then retry.
#[no_mangle]
pub extern "C" fn ps_solve(
    hash_tol: u16,
    quorum: u32,
    brightest: u32,
    max_bucket_reads: u32,
    out: *mut SolveResultC,
) -> i32 {
    let t = std::time::Instant::now();
    SESSION.with(|s| {
        let mut sess = s.borrow_mut();
        let Some(sess) = sess.as_mut() else {
            return set_err("session not initialized");
        };
        if sess.sources.len() < 4 {
            return set_err("need at least 4 extracted sources");
        }
        let (w, h) = match sess.image.as_ref() {
            Some(img) => (img.width as f64, img.height as f64),
            None => return set_err("no image loaded"),
        };
        let cfg = SolveConfig {
            hash_tol,
            quorum: quorum.max(2) as usize,
            brightest: brightest.max(8) as usize,
            max_bucket_reads: max_bucket_reads.max(16) as usize,
            ..Default::default()
        };
        let result = solve(&sess.sources, w, h, &mut sess.store, &cfg);
        match result {
            Ok(o) => {
                let (fovx, fovy) = o.wcs.fov_xy_rad();
                let rc = SolveResultC {
                    ra0_deg: o.wcs.ra0 / DEG,
                    dec0_deg: o.wcs.dec0 / DEG,
                    fov_x_deg: fovx / DEG,
                    fov_y_deg: fovy / DEG,
                    fov_diag_deg: o.wcs.fov_diagonal_rad() / DEG,
                    rotation_deg: o.wcs.rotation_east_of_north() / DEG,
                    plate_scale_as: o.wcs.plate_scale_arcsec(),
                    rms_arcsec: o.wcs.rms_arcsec,
                    n_inliers: o.wcs.n_inliers as u32,
                    n_matches: o.matches.len() as u32,
                    n_sources: sess.sources.len() as u32,
                    buckets_read: o.buckets_read as u32,
                    index_bytes: o.index_bytes_read,
                    quads_tried: o.quads_tried as u32,
                    solve_ms: t.elapsed().as_secs_f64() * 1000.0,
                };
                unsafe {
                    out.write(rc);
                }
                sess.outcome = Some(o);
                0
            }
            Err(e) => set_err(format!("{e}")),
        }
    })
}

/// Map a pixel to sky coordinates after a successful solve.
#[no_mangle]
pub extern "C" fn ps_pixel_to_radec(x: f64, y: f64, ra: *mut f64, dec: *mut f64) -> i32 {
    with_session(|sess| match &sess.outcome {
        Some(o) => {
            let (r, d) = o.wcs.pixel_to_radec(x, y);
            unsafe {
                ra.write(r / DEG);
                dec.write(d / DEG);
            }
            0
        }
        None => set_err("no solution yet"),
    })
    .unwrap_or(-1)
}

/// Return matched `(source_index, catalog_id)` pairs (two u32 each).
#[no_mangle]
pub extern "C" fn ps_get_matches(out: *mut u32, cap_pairs: u32) -> i32 {
    SESSION.with(|s| {
        let sess = s.borrow();
        let Some(sess) = sess.as_ref() else {
            return set_err("session not initialized");
        };
        let Some(o) = &sess.outcome else {
            return set_err("no solution yet");
        };
        let need = o.matches.len() as u32;
        if cap_pairs < need {
            return -(need as i32);
        }
        for (i, (im, ca)) in o.matches.iter().enumerate() {
            unsafe {
                *out.add(2 * i) = *im as u32;
                *out.add(2 * i + 1) = *ca;
            }
        }
        need as i32
    })
}

/// Build (or rebuild) a KD-tree over the resident catalog and return the
/// catalog star id nearest to a sky position — used by the UI for labeling
/// arbitrary image points, not just matched ones.
#[no_mangle]
pub extern "C" fn ps_nearest_catalog_star(
    ra_deg: f64,
    dec_deg: f64,
    search_radius_deg: f64,
    out_id: *mut u32,
    out_sep_as: *mut f64,
) -> i32 {
    SESSION.with(|s| {
        let sess = s.borrow();
        let Some(sess) = sess.as_ref() else {
            return set_err("session not initialized");
        };
        let stars: Vec<CatStar> = sess.store.catalog.values().copied().collect();
        if stars.is_empty() {
            return set_err("empty catalog");
        }
        let tree = KdTree::build(stars);
        let q = Vec3::from_radec(ra_deg * DEG, dec_deg * DEG);
        let Some(near) = tree.nearest(q) else {
            return set_err("no nearest star");
        };
        let sep = near.v.separation(q);
        if sep > search_radius_deg * DEG {
            return set_err("no star within search radius");
        }
        unsafe {
            out_id.write(near.id);
            out_sep_as.write(sep / platersol_core::ARCSEC);
        }
        0
    })
}

/// Reset all fed buckets (e.g. when switching index sets).
#[no_mangle]
pub extern "C" fn ps_clear_buckets() -> i32 {
    SESSION.with(|s| match s.borrow_mut().as_mut() {
        Some(sess) => {
            sess.store.buckets.clear();
            sess.store.missing.clear();
            sess.store.bytes_fetched = 0;
            0
        }
        None => set_err("session not initialized"),
    })
}

/// Planning pass: run the solver purely to learn which index buckets it would
/// read, without requiring them resident. Missing buckets are recorded; call
/// `ps_take_missing_buckets` afterward to fetch and feed them, then re-solve.
/// Returns the number of distinct missing buckets, or a negative error.
#[no_mangle]
pub extern "C" fn ps_plan(hash_tol: u16, quorum: u32, brightest: u32, max_bucket_reads: u32) -> i32 {
    let result = SESSION.with(|s| {
        let mut sess = s.borrow_mut();
        let Some(sess) = sess.as_mut() else {
            return Err(set_err("session not initialized"));
        };
        if sess.sources.len() < 4 {
            return Err(set_err("need at least 4 extracted sources"));
        }
        sess.store.missing.clear();
        sess.store.planning = true;
        let (w, h) = match sess.image.as_ref() {
            Some(img) => (img.width as f64, img.height as f64),
            None => {
                sess.store.planning = false;
                return Err(set_err("no image loaded"));
            }
        };
        let cfg = SolveConfig {
            hash_tol,
            quorum: quorum.max(2) as usize,
            brightest: brightest.max(8) as usize,
            max_bucket_reads: max_bucket_reads.max(16) as usize,
            ..Default::default()
        };
        // In planning mode every missing bucket is recorded but returned
        // empty, so this run walks the full probe path and enumerates every
        // chunk a real solve would request. Its Ok/Err result is irrelevant.
        let _ = solve(&sess.sources, w, h, &mut sess.store, &cfg);
        sess.store.planning = false;
        Ok(sess.store.missing.len() as i32)
    });
    result.unwrap_or_else(|e| e)
}

/// Copy recorded missing-bucket ids into a caller array of `(bin u32,
/// bucket u32)` pairs. Returns the required pair count (negative if `cap` is
/// too small), so the host can size its buffer.
#[no_mangle]
pub extern "C" fn ps_take_missing_buckets(out: *mut u32, cap_pairs: u32) -> i32 {
    SESSION.with(|s| {
        let sess = s.borrow();
        let Some(sess) = sess.as_ref() else {
            return set_err("session not initialized");
        };
        let need = sess.store.missing.len() as u32;
        if cap_pairs < need {
            return -(need as i32);
        }
        for (i, (bin, bucket)) in sess.store.missing.iter().enumerate() {
            unsafe {
                *out.add(2 * i) = *bin as u32;
                *out.add(2 * i + 1) = *bucket as u32;
            }
        }
        need as i32
    })
}

/// Count of missing buckets currently recorded (for UI progress metering).
#[no_mangle]
pub extern "C" fn ps_missing_count() -> u32 {
    SESSION.with(|s| {
        s.borrow()
            .as_ref()
            .map(|sess| sess.store.missing.len() as u32)
            .unwrap_or(0)
    })
}

/// Expose the fitted WCS coefficients (for serialization / debugging).
#[no_mangle]
pub extern "C" fn ps_wcs_coeffs(a: *mut f64, b: *mut f64, ar: *mut f64, br: *mut f64) -> i32 {
    with_session(|sess| match &sess.outcome {
        Some(o) => {
            let w: &WcsSolution = &o.wcs;
            unsafe {
                for i in 0..3 {
                    a.add(i).write(w.a[i]);
                    b.add(i).write(w.b[i]);
                }
                ar.write(w.a_rad);
                br.write(w.b_rad);
            }
            0
        }
        None => set_err("no solution yet"),
    })
    .unwrap_or(-1)
}
