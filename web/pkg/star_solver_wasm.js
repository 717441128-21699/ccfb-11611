let wasm;

const cachedTextDecoder = (typeof TextDecoder !== 'undefined' ? new TextDecoder('utf-8', { ignoreBOM: true, fatal: true }) : { decode: () => { throw Error('TextDecoder not available') } } );

if (typeof TextDecoder !== 'undefined') { cachedTextDecoder.decode(); };

let cachedUint8ArrayMemory0 = null;

function getUint8ArrayMemory0() {
    if (cachedUint8ArrayMemory0 === null || cachedUint8ArrayMemory0.byteLength === 0) {
        cachedUint8ArrayMemory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8ArrayMemory0;
}

function getStringFromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return cachedTextDecoder.decode(getUint8ArrayMemory0().subarray(ptr, ptr + len));
}

function _assertClass(instance, klass) {
    if (!(instance instanceof klass)) {
        throw new Error(`expected instance of ${klass.name}`);
    }
}

let WASM_VECTOR_LEN = 0;

function passArray8ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 1, 1) >>> 0;
    getUint8ArrayMemory0().set(arg, ptr / 1);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}
/**
 * 用内存索引对一张灰度图做盲解。
 *
 * `index` 为 [`build_index`] 返回的句柄；`scale_lo/hi` 是先验比例尺范围
 * （角秒/像元），未知焦距时给宽区间（如 0.3~60）。
 * @param {WasmIndex} index
 * @param {Uint8Array} gray
 * @param {number} width
 * @param {number} height
 * @param {number} scale_lo
 * @param {number} scale_hi
 * @returns {SolveResult}
 */
export function solve_image(index, gray, width, height, scale_lo, scale_hi) {
    _assertClass(index, WasmIndex);
    const ptr0 = passArray8ToWasm0(gray, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.solve_image(index.__wbg_ptr, ptr0, len0, width, height, scale_lo, scale_hi);
    return SolveResult.__wrap(ret);
}

/**
 * 内置合成星空自检（无需任何输入数据）。
 * @returns {SolveResult}
 */
export function run_selftest() {
    const ret = wasm.run_selftest();
    return SolveResult.__wrap(ret);
}

/**
 * 暴露镜像编码（JS 侧若要自行实现查询可用）。
 * @param {number} code
 * @returns {number}
 */
export function mirror_code_of(code) {
    const ret = wasm.mirror_code_of(code);
    return ret;
}

const cachedTextEncoder = (typeof TextEncoder !== 'undefined' ? new TextEncoder('utf-8') : { encode: () => { throw Error('TextEncoder not available') } } );

const encodeString = (typeof cachedTextEncoder.encodeInto === 'function'
    ? function (arg, view) {
    return cachedTextEncoder.encodeInto(arg, view);
}
    : function (arg, view) {
    const buf = cachedTextEncoder.encode(arg);
    view.set(buf);
    return {
        read: arg.length,
        written: buf.length
    };
});

function passStringToWasm0(arg, malloc, realloc) {

    if (realloc === undefined) {
        const buf = cachedTextEncoder.encode(arg);
        const ptr = malloc(buf.length, 1) >>> 0;
        getUint8ArrayMemory0().subarray(ptr, ptr + buf.length).set(buf);
        WASM_VECTOR_LEN = buf.length;
        return ptr;
    }

    let len = arg.length;
    let ptr = malloc(len, 1) >>> 0;

    const mem = getUint8ArrayMemory0();

    let offset = 0;

    for (; offset < len; offset++) {
        const code = arg.charCodeAt(offset);
        if (code > 0x7F) break;
        mem[ptr + offset] = code;
    }

    if (offset !== len) {
        if (offset !== 0) {
            arg = arg.slice(offset);
        }
        ptr = realloc(ptr, len, len = offset + arg.length * 3, 1) >>> 0;
        const view = getUint8ArrayMemory0().subarray(ptr + offset, ptr + len);
        const ret = encodeString(arg, view);

        offset += ret.written;
        ptr = realloc(ptr, len, offset, 1) >>> 0;
    }

    WASM_VECTOR_LEN = offset;
    return ptr;
}

function takeFromExternrefTable0(idx) {
    const value = wasm.__wbindgen_export_0.get(idx);
    wasm.__externref_table_dealloc(idx);
    return value;
}
/**
 * 用解出的 WCS 把 DSO 表投影到画幅，返回落在图内天体的 JSON 数组字符串。
 *
 * 入参 `wcs_json` 取自 solve_image 结果里的 `wcs` 字段；`dso_csv` 为
 * `name,ra_deg,dec_deg,kind,mag,major_arcmin,minor_arcmin,pa_deg` 格式。
 * 输出元素形如 `{"name","kind","x","y","ra_deg","dec_deg","mag"}`。
 * @param {string} wcs_json
 * @param {string} dso_csv
 * @param {number} width
 * @param {number} height
 * @returns {string}
 */
export function project_dso(wcs_json, dso_csv, width, height) {
    let deferred4_0;
    let deferred4_1;
    try {
        const ptr0 = passStringToWasm0(wcs_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(dso_csv, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.project_dso(ptr0, len0, ptr1, len1, width, height);
        var ptr3 = ret[0];
        var len3 = ret[1];
        if (ret[3]) {
            ptr3 = 0; len3 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred4_0 = ptr3;
        deferred4_1 = len3;
        return getStringFromWasm0(ptr3, len3);
    } finally {
        wasm.__wbindgen_free(deferred4_0, deferred4_1, 1);
    }
}

let cachedDataViewMemory0 = null;

function getDataViewMemory0() {
    if (cachedDataViewMemory0 === null || cachedDataViewMemory0.buffer.detached === true || (cachedDataViewMemory0.buffer.detached === undefined && cachedDataViewMemory0.buffer !== wasm.memory.buffer)) {
        cachedDataViewMemory0 = new DataView(wasm.memory.buffer);
    }
    return cachedDataViewMemory0;
}

function getArrayJsValueFromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    const mem = getDataViewMemory0();
    const result = [];
    for (let i = ptr; i < ptr + 4 * len; i += 4) {
        result.push(wasm.__wbindgen_export_0.get(mem.getUint32(i, true)));
    }
    wasm.__externref_drop_slice(ptr, len);
    return result;
}
/**
 * 星点检测入口：灰度字节（行优先，0..255）+ 宽高。
 * @param {Uint8Array} gray
 * @param {number} width
 * @param {number} height
 * @param {number} k_sigma
 * @param {number} max_stars
 * @returns {DetectedStar[]}
 */
export function detect(gray, width, height, k_sigma, max_stars) {
    const ptr0 = passArray8ToWasm0(gray, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.detect(ptr0, len0, width, height, k_sigma, max_stars);
    if (ret[3]) {
        throw takeFromExternrefTable0(ret[2]);
    }
    var v2 = getArrayJsValueFromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
    return v2;
}

/**
 * 背景估计（供 UI 显示与手动调阈值）。
 * @param {Uint8Array} gray
 * @param {number} width
 * @param {number} height
 * @returns {BackgroundStats}
 */
export function background(gray, width, height) {
    const ptr0 = passArray8ToWasm0(gray, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.background(ptr0, len0, width, height);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return BackgroundStats.__wrap(ret[0]);
}

/**
 * 由 CSV 星表文本在内存中构建索引（小星表/演示用，返回各文件的字节表）。
 * @param {string} csv_text
 * @param {number} mag_limit
 * @returns {WasmIndex}
 */
export function build_index(csv_text, mag_limit) {
    const ptr0 = passStringToWasm0(csv_text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.build_index(ptr0, len0, mag_limit);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return WasmIndex.__wrap(ret[0]);
}

function getArrayU8FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint8ArrayMemory0().subarray(ptr / 1, ptr / 1 + len);
}

const BackgroundStatsFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_backgroundstats_free(ptr >>> 0, 1));

export class BackgroundStats {

    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(BackgroundStats.prototype);
        obj.__wbg_ptr = ptr;
        BackgroundStatsFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        BackgroundStatsFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_backgroundstats_free(ptr, 0);
    }
    /**
     * @returns {number}
     */
    get level() {
        const ret = wasm.__wbg_get_backgroundstats_level(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {number} arg0
     */
    set level(arg0) {
        wasm.__wbg_set_backgroundstats_level(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {number}
     */
    get sigma() {
        const ret = wasm.__wbg_get_backgroundstats_sigma(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {number} arg0
     */
    set sigma(arg0) {
        wasm.__wbg_set_backgroundstats_sigma(this.__wbg_ptr, arg0);
    }
}

const DetectedStarFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_detectedstar_free(ptr >>> 0, 1));
/**
 * 检测到的星点（JS 友好的扁平结构）。
 */
export class DetectedStar {

    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(DetectedStar.prototype);
        obj.__wbg_ptr = ptr;
        DetectedStarFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        DetectedStarFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_detectedstar_free(ptr, 0);
    }
    /**
     * @returns {number}
     */
    get x() {
        const ret = wasm.__wbg_get_backgroundstats_level(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {number} arg0
     */
    set x(arg0) {
        wasm.__wbg_set_backgroundstats_level(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {number}
     */
    get y() {
        const ret = wasm.__wbg_get_backgroundstats_sigma(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {number} arg0
     */
    set y(arg0) {
        wasm.__wbg_set_backgroundstats_sigma(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {number}
     */
    get flux() {
        const ret = wasm.__wbg_get_detectedstar_flux(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {number} arg0
     */
    set flux(arg0) {
        wasm.__wbg_set_detectedstar_flux(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {number}
     */
    get area() {
        const ret = wasm.__wbg_get_detectedstar_area(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @param {number} arg0
     */
    set area(arg0) {
        wasm.__wbg_set_detectedstar_area(this.__wbg_ptr, arg0);
    }
}

const SolveResultFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_solveresult_free(ptr >>> 0, 1));
/**
 * 完整解算的 JS 结果。
 */
export class SolveResult {

    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(SolveResult.prototype);
        obj.__wbg_ptr = ptr;
        SolveResultFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        SolveResultFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_solveresult_free(ptr, 0);
    }
    /**
     * @returns {boolean}
     */
    get ok() {
        const ret = wasm.solveresult_ok(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * 成功时的完整 JSON（中心、比例尺、滚动、视场、内点、FITS 头、匹配星）。
     * @returns {string}
     */
    get json() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.solveresult_json(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * @returns {string}
     */
    get error() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.solveresult_error(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
}

const WasmIndexFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_wasmindex_free(ptr >>> 0, 1));
/**
 * 构建产物的持有者（不透明句柄，JS 不解码内部布局）。
 */
export class WasmIndex {

    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(WasmIndex.prototype);
        obj.__wbg_ptr = ptr;
        WasmIndexFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        WasmIndexFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_wasmindex_free(ptr, 0);
    }
    /**
     * 第 band 个尺度带文件字节。
     * @param {number} band
     * @returns {Uint8Array | undefined}
     */
    band_bytes(band) {
        const ret = wasm.wasmindex_band_bytes(this.__wbg_ptr, band);
        let v1;
        if (ret[0] !== 0) {
            v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
            wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        }
        return v1;
    }
    /**
     * @returns {number}
     */
    get band_count() {
        const ret = wasm.wasmindex_band_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * 合并后的星表块字节（含块头，可直接落盘为 stars.bin）。
     * @returns {Uint8Array}
     */
    star_bytes() {
        const ret = wasm.wasmindex_star_bytes(this.__wbg_ptr);
        var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v1;
    }
    /**
     * @returns {number}
     */
    get star_count() {
        const ret = wasm.wasmindex_star_count(this.__wbg_ptr);
        return ret;
    }
    /**
     * manifest 的 JSON 文本（供 UI 展示/下载）。
     * @returns {string}
     */
    manifest_json() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.wasmindex_manifest_json(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
}

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);

            } catch (e) {
                if (module.headers.get('Content-Type') != 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else {
                    throw e;
                }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);

    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };

        } else {
            return instance;
        }
    }
}

function __wbg_get_imports() {
    const imports = {};
    imports.wbg = {};
    imports.wbg.__wbg_detectedstar_new = function(arg0) {
        const ret = DetectedStar.__wrap(arg0);
        return ret;
    };
    imports.wbg.__wbindgen_init_externref_table = function() {
        const table = wasm.__wbindgen_export_0;
        const offset = table.grow(4);
        table.set(0, undefined);
        table.set(offset + 0, undefined);
        table.set(offset + 1, null);
        table.set(offset + 2, true);
        table.set(offset + 3, false);
        ;
    };
    imports.wbg.__wbindgen_string_new = function(arg0, arg1) {
        const ret = getStringFromWasm0(arg0, arg1);
        return ret;
    };
    imports.wbg.__wbindgen_throw = function(arg0, arg1) {
        throw new Error(getStringFromWasm0(arg0, arg1));
    };

    return imports;
}

function __wbg_init_memory(imports, memory) {

}

function __wbg_finalize_init(instance, module) {
    wasm = instance.exports;
    __wbg_init.__wbindgen_wasm_module = module;
    cachedDataViewMemory0 = null;
    cachedUint8ArrayMemory0 = null;


    wasm.__wbindgen_start();
    return wasm;
}

function initSync(module) {
    if (wasm !== undefined) return wasm;


    if (typeof module !== 'undefined') {
        if (Object.getPrototypeOf(module) === Object.prototype) {
            ({module} = module)
        } else {
            console.warn('using deprecated parameters for `initSync()`; pass a single object instead')
        }
    }

    const imports = __wbg_get_imports();

    __wbg_init_memory(imports);

    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }

    const instance = new WebAssembly.Instance(module, imports);

    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(module_or_path) {
    if (wasm !== undefined) return wasm;


    if (typeof module_or_path !== 'undefined') {
        if (Object.getPrototypeOf(module_or_path) === Object.prototype) {
            ({module_or_path} = module_or_path)
        } else {
            console.warn('using deprecated parameters for the initialization function; pass a single object instead')
        }
    }

    if (typeof module_or_path === 'undefined') {
        module_or_path = new URL('star_solver_wasm_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    __wbg_init_memory(imports);

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync };
export default __wbg_init;
