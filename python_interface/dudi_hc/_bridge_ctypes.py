from __future__ import annotations
import ctypes as C
from pathlib import Path
import numpy as np

# Load the shared library placed by build_ctypes_bridge.sh
_lib_path = Path(__file__).with_name("libpy_dudihc_bridge.so")
if not _lib_path.exists():
    raise ImportError(
        f"Shared library not found: {_lib_path}. "
        f"Run: bash python_interface/fortran_bridge/build_ctypes_bridge.sh"
    )
_lib = C.CDLL(str(_lib_path))

# ---- method IDs (must match Fortran batching.f90) ----

METHOD_SIMPLE_EXPANSION = 1
METHOD_DELTA_EJECTION = 2
METHOD_V_INTEGRATION = 3

# ---- helpers / ctypes dtypes ----

Vec3d = np.ctypeslib.ndpointer(dtype=np.float64, shape=(3,), flags="C_CONTIGUOUS")
Vec1d = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
Int1d = np.ctypeslib.ndpointer(dtype=np.int32,  flags="C_CONTIGUOUS")


def _as_vec3(a) -> np.ndarray:
    arr = np.asarray(a, dtype=np.float64)
    if arr.shape != (3,):
        raise ValueError(f"expected shape (3,), got {arr.shape}")
    return np.ascontiguousarray(arr)


def _as_1d_f64(a, name: str = "array") -> np.ndarray:
    arr = np.asarray(a, dtype=np.float64)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1D, got shape {arr.shape}")
    return np.ascontiguousarray(arr)


def _as_2d_f64(a, name: str = "array") -> np.ndarray:
    arr = np.asarray(a, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be 2D, got shape {arr.shape}")
    return np.ascontiguousarray(arr)


def _as_1d_i32(a, name: str = "array") -> np.ndarray:
    arr = np.asarray(a, dtype=np.int32)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1D, got shape {arr.shape}")
    return np.ascontiguousarray(arr)


def _c_int(b: bool | int) -> int:
    return 1 if bool(b) else 0


# ======================================================================
#  scalar prototypes
# ======================================================================

# void py_hc_v_integration(double r, double alpha, double beta, double rvec[3],
#   double sr, double sAlphaM, double sBetaM, double sRR[3], double zeta, double eta,
#   double axis[3], int eject_distr, int ud_shape, double umin, double umax,
#   double Nparticles, double Tj, double dtau,
#   double comet_coords[3], double comet_vvec[3], double comet_v,
#   double muR, double tnow, double Rast_AU, int pericenter, double* density_out)
_lib.py_hc_v_integration.argtypes = [
    # point
    C.c_double, C.c_double, C.c_double, Vec3d,
    # source core
    C.c_double, C.c_double, C.c_double, Vec3d, C.c_double, C.c_double,
    Vec3d, C.c_int, C.c_int, C.c_double, C.c_double,
    # source extras
    C.c_double, C.c_double, C.c_double,    # Nparticles, Tj, dtau
    # comet
    Vec3d, Vec3d, C.c_double,
    # scalars
    C.c_double, C.c_double, C.c_double, C.c_int,
    # out
    C.POINTER(C.c_double),
]
_lib.py_hc_v_integration.restype = None

_lib.py_hc_delta_ejection.argtypes = [
    # point
    C.c_double, C.c_double, C.c_double, Vec3d,
    # source core
    C.c_double, C.c_double, C.c_double, Vec3d, C.c_double, C.c_double,
    Vec3d, C.c_int, C.c_int, C.c_double, C.c_double,
    # source extras
    C.c_double, C.c_double, C.c_double,    # Nparticles, Tj, dtau
    # comet
    Vec3d, Vec3d, C.c_double,
    # scalars
    C.c_double, C.c_double, C.c_double,
    # out
    C.POINTER(C.c_double),
]
_lib.py_hc_delta_ejection.restype = None

_lib.py_hc_simple_expansion.argtypes = [
    # point
    C.c_double, C.c_double, C.c_double, Vec3d,
    # source core
    C.c_double, C.c_double, C.c_double, Vec3d, C.c_double, C.c_double,
    Vec3d, C.c_int, C.c_int, C.c_double, C.c_double,
    # source extras
    C.c_double, C.c_double, C.c_double,    # Nparticles, Tj, dtau
    # cloud & scalar
    Vec3d, C.c_double,
    # out
    C.POINTER(C.c_double),
]
_lib.py_hc_simple_expansion.restype = None


# ======================================================================
#  scalar thin wrappers
# ======================================================================

def call_v_integration(
    *,
    point_r: float, point_alpha: float, point_beta: float, point_rvector,
    src_r: float, src_alphaM: float, src_betaM: float, src_rrM,
    src_zeta: float, src_eta: float, src_axis, src_eject_distr: int,
    src_ud_shape: int, src_umin: float, src_umax: float,
    src_Nparticles: float, src_Tj: float, src_dtau: float,
    comet_coords, comet_vastvec, comet_vast: float,
    muR: float, tnow: float, Rast_AU: float, pericenter: bool,
) -> float:
    out = C.c_double()
    _lib.py_hc_v_integration(
        float(point_r), float(point_alpha), float(point_beta), _as_vec3(point_rvector),
        float(src_r), float(src_alphaM), float(src_betaM), _as_vec3(src_rrM),
        float(src_zeta), float(src_eta), _as_vec3(src_axis),
        int(src_eject_distr), int(src_ud_shape), float(src_umin), float(src_umax),
        float(src_Nparticles), float(src_Tj), float(src_dtau),
        _as_vec3(comet_coords), _as_vec3(comet_vastvec), float(comet_vast),
        float(muR), float(tnow), float(Rast_AU), _c_int(pericenter),
        C.byref(out),
    )
    return float(out.value)


def call_delta_ejection(
    *,
    point_r: float, point_alpha: float, point_beta: float, point_rvector,
    src_r: float, src_alphaM: float, src_betaM: float, src_rrM,
    src_zeta: float, src_eta: float, src_axis, src_eject_distr: int,
    src_ud_shape: int, src_umin: float, src_umax: float,
    src_Nparticles: float, src_Tj: float, src_dtau: float,
    comet_coords, comet_vastvec, comet_vast: float,
    muR: float, dt: float, Rast_AU: float,
) -> float:
    out = C.c_double()
    _lib.py_hc_delta_ejection(
        float(point_r), float(point_alpha), float(point_beta), _as_vec3(point_rvector),
        float(src_r), float(src_alphaM), float(src_betaM), _as_vec3(src_rrM),
        float(src_zeta), float(src_eta), _as_vec3(src_axis),
        int(src_eject_distr), int(src_ud_shape), float(src_umin), float(src_umax),
        float(src_Nparticles), float(src_Tj), float(src_dtau),
        _as_vec3(comet_coords), _as_vec3(comet_vastvec), float(comet_vast),
        float(muR), float(dt), float(Rast_AU), C.byref(out),
    )
    return float(out.value)


def call_simple_expansion(
    *,
    point_r: float, point_alpha: float, point_beta: float, point_rvector,
    src_r: float, src_alphaM: float, src_betaM: float, src_rrM,
    src_zeta: float, src_eta: float, src_axis, src_eject_distr: int,
    src_ud_shape: int, src_umin: float, src_umax: float,
    src_Nparticles: float, src_Tj: float, src_dtau: float,
    cloudcentr, dt: float,
) -> float:
    out = C.c_double()
    _lib.py_hc_simple_expansion(
        float(point_r), float(point_alpha), float(point_beta), _as_vec3(point_rvector),
        float(src_r), float(src_alphaM), float(src_betaM), _as_vec3(src_rrM),
        float(src_zeta), float(src_eta), _as_vec3(src_axis),
        int(src_eject_distr), int(src_ud_shape), float(src_umin), float(src_umax),
        float(src_Nparticles), float(src_Tj), float(src_dtau),
        _as_vec3(cloudcentr), float(dt), C.byref(out),
    )
    return float(out.value)


# ======================================================================
#  batched prototypes
#  (Fortran side: py_hc_batch_points / py_hc_batch_sources must match)
# ======================================================================

# void py_hc_batch_points(
#   int n_points, double density_out[n_points],
#   double point_r[n_points], double point_alpha[n_points],
#   double point_beta[n_points], double point_rvec[3*n_points],
#   double sr, double sAlphaM, double sBetaM, double sRR[3],
#   double zeta, double eta, double axis[3], int eject_distr, int ud_shape,
#   double umin, double umax, double Nparticles, double Tj, double dtau,
#   double comet_coords[3], double comet_vvec[3], double comet_v,
#   double muR, double tnow, double dt, double Rast_AU, int pericenter,
#   double cloudcentr[3], int method_id)
_lib.py_hc_batch_points.argtypes = [
    C.c_int,      # n_points
    Vec1d,        # density_out
    Vec1d, Vec1d, Vec1d,  # point_r, alpha, beta
    Vec1d,        # point_rvector flat (len = 3*n_points)
    # source (same as scalar)
    C.c_double, C.c_double, C.c_double, Vec3d, C.c_double, C.c_double,
    Vec3d, C.c_int, C.c_int, C.c_double, C.c_double,
    C.c_double, C.c_double, C.c_double,
    # comet
    Vec3d, Vec3d, C.c_double,
    # scalars & extra
    C.c_double, C.c_double, C.c_double, C.c_double, C.c_int,
    Vec3d, C.c_int,
]
_lib.py_hc_batch_points.restype = None

# void py_hc_batch_sources(
#   int n_sources, double density_out[n_sources],
#   double point_r, double point_alpha, double point_beta, double point_rvec[3],
#   double src_r[n], double src_alphaM[n], double src_betaM[n],
#   double src_rrM[3*n], double zeta[n], double eta[n],
#   double axis[3*n], int eject_distr[n], int ud_shape[n],
#   double umin[n], double umax[n], double Nparticles[n],
#   double Tj[n], double dtau[n],
#   double comet_coords[3], double comet_vvec[3], double comet_v,
#   double muR, double tnow, double dt, double Rast_AU, int pericenter,
#   double cloudcentr[3], int method_id)
_lib.py_hc_batch_sources.argtypes = [
    C.c_int,  # n_sources
    Vec1d,    # density_out
    # point (single)
    C.c_double, C.c_double, C.c_double, Vec3d,
    # source arrays
    Vec1d, Vec1d, Vec1d,  # src_r, src_alphaM, src_betaM
    Vec1d, Vec1d, Vec1d,  # src_rrM_flat, zeta, eta
    Vec1d, Int1d, Int1d,  # axis_flat, eject_distr, ud_shape
    Vec1d, Vec1d,         # umin, umax
    Vec1d, Vec1d, Vec1d,  # Nparticles, Tj, dtau
    # comet
    Vec3d, Vec3d, C.c_double,
    # scalars
    C.c_double, C.c_double, C.c_double, C.c_double, C.c_int,
    Vec3d, C.c_int,
]
_lib.py_hc_batch_sources.restype = None


# ======================================================================
#  batched thin wrappers (return numpy arrays)
# ======================================================================

def call_batch_points(
    *,
    point_r,
    point_alpha,
    point_beta,
    point_rvector,
    src_r: float, src_alphaM: float, src_betaM: float, src_rrM,
    src_zeta: float, src_eta: float, src_axis, src_eject_distr: int,
    src_ud_shape: int, src_umin: float, src_umax: float,
    src_Nparticles: float, src_Tj: float, src_dtau: float,
    comet_coords, comet_vastvec, comet_vast: float,
    muR: float, tnow: float, dt: float, Rast_AU: float,
    pericenter: bool, cloudcentr, method_id: int,
) -> np.ndarray:
    """
    Low-level batched wrapper over hc_DUDI_batch_points.

    point_* are 1D arrays, point_rvector is (N, 3).
    Returns density array of length N.
    """
    r = _as_1d_f64(point_r, "point_r")
    alpha = _as_1d_f64(point_alpha, "point_alpha")
    beta = _as_1d_f64(point_beta, "point_beta")
    if not (r.size == alpha.size == beta.size):
        raise ValueError("point_r, point_alpha, point_beta must have same length")
    n = int(r.size)

    rvec2d = _as_2d_f64(point_rvector, "point_rvector")
    if rvec2d.shape[0] != n or rvec2d.shape[1] != 3:
        raise ValueError(
            f"point_rvector must have shape (N,3) with N={n}, got {rvec2d.shape}"
        )
    rvec_flat = np.ascontiguousarray(rvec2d.reshape(-1))

    density = np.empty(n, dtype=np.float64)

    _lib.py_hc_batch_points(
        n,
        density,
        r, alpha, beta,
        rvec_flat,
        float(src_r), float(src_alphaM), float(src_betaM), _as_vec3(src_rrM),
        float(src_zeta), float(src_eta), _as_vec3(src_axis),
        int(src_eject_distr), int(src_ud_shape),
        float(src_umin), float(src_umax),
        float(src_Nparticles), float(src_Tj), float(src_dtau),
        _as_vec3(comet_coords), _as_vec3(comet_vastvec), float(comet_vast),
        float(muR), float(tnow), float(dt), float(Rast_AU), _c_int(pericenter),
        _as_vec3(cloudcentr), int(method_id),
    )
    return density


def call_batch_sources(
    *,
    point_r: float,
    point_alpha: float,
    point_beta: float,
    point_rvector,
    src_r,
    src_alphaM,
    src_betaM,
    src_rrM,
    src_zeta,
    src_eta,
    src_axis,
    src_eject_distr,
    src_ud_shape,
    src_umin,
    src_umax,
    src_Nparticles,
    src_Tj,
    src_dtau,
    comet_coords,
    comet_vastvec,
    comet_vast: float,
    muR: float,
    tnow: float,
    dt: float,
    Rast_AU: float,
    pericenter: bool,
    cloudcentr,
    method_id: int,
) -> np.ndarray:
    """
    Low-level batched wrapper over hc_DUDI_batch_sources.

    src_* are 1D arrays (or (N,3) for vectors). Returns density array of length N.
    """
    r_arr = _as_1d_f64(src_r, "src_r")
    alphaM_arr = _as_1d_f64(src_alphaM, "src_alphaM")
    betaM_arr = _as_1d_f64(src_betaM, "src_betaM")
    if not (r_arr.size == alphaM_arr.size == betaM_arr.size):
        raise ValueError("src_r, src_alphaM, src_betaM must have same length")
    n = int(r_arr.size)

    rrM2d = _as_2d_f64(src_rrM, "src_rrM")
    axis2d = _as_2d_f64(src_axis, "src_axis")
    if rrM2d.shape != (n, 3):
        raise ValueError(f"src_rrM must have shape (N,3), got {rrM2d.shape}")
    if axis2d.shape != (n, 3):
        raise ValueError(f"src_axis must have shape (N,3), got {axis2d.shape}")
    rrM_flat = np.ascontiguousarray(rrM2d.reshape(-1))
    axis_flat = np.ascontiguousarray(axis2d.reshape(-1))

    zeta_arr = _as_1d_f64(src_zeta, "src_zeta")
    eta_arr = _as_1d_f64(src_eta, "src_eta")
    if not (zeta_arr.size == eta_arr.size == n):
        raise ValueError("src_zeta and src_eta must have length N")

    eject_arr = _as_1d_i32(src_eject_distr, "src_eject_distr")
    udshape_arr = _as_1d_i32(src_ud_shape, "src_ud_shape")
    if not (eject_arr.size == udshape_arr.size == n):
        raise ValueError("src_eject_distr and src_ud_shape must have length N")

    umin_arr = _as_1d_f64(src_umin, "src_umin")
    umax_arr = _as_1d_f64(src_umax, "src_umax")
    Np_arr = _as_1d_f64(src_Nparticles, "src_Nparticles")
    Tj_arr = _as_1d_f64(src_Tj, "src_Tj")
    dtau_arr = _as_1d_f64(src_dtau, "src_dtau")
    for name, arr in [
        ("src_umin", umin_arr),
        ("src_umax", umax_arr),
        ("src_Nparticles", Np_arr),
        ("src_Tj", Tj_arr),
        ("src_dtau", dtau_arr),
    ]:
        if arr.size != n:
            raise ValueError(f"{name} must have length N={n}")

    density = np.empty(n, dtype=np.float64)

    _lib.py_hc_batch_sources(
        n,
        density,
        float(point_r), float(point_alpha), float(point_beta), _as_vec3(point_rvector),
        r_arr, alphaM_arr, betaM_arr,
        rrM_flat, zeta_arr, eta_arr,
        axis_flat, eject_arr, udshape_arr,
        umin_arr, umax_arr,
        Np_arr, Tj_arr, dtau_arr,
        _as_vec3(comet_coords), _as_vec3(comet_vastvec), float(comet_vast),
        float(muR), float(tnow), float(dt), float(Rast_AU), _c_int(pericenter),
        _as_vec3(cloudcentr), int(method_id),
    )
    return density


# ======================================================================
#  sizes
# ======================================================================

_lib.py_get_nlats.restype = C.c_int
_lib.py_get_nlons.restype = C.c_int


def get_nlats() -> int:
    return int(_lib.py_get_nlats())


def get_nlons() -> int:
    return int(_lib.py_get_nlons())


# ======================================================================
#  lon bounds
# ======================================================================

_lib.py_set_lon_bounds.argtypes = [C.c_double, C.c_double]
_lib.py_get_lon_bounds.argtypes = [C.POINTER(C.c_double), C.POINTER(C.c_double)]


def set_lon_bounds(lonmin: float, lonmax: float) -> None:
    _lib.py_set_lon_bounds(float(lonmin), float(lonmax))


def get_lon_bounds() -> tuple[float, float]:
    a = C.c_double()
    b = C.c_double()
    _lib.py_get_lon_bounds(C.byref(a), C.byref(b))
    return float(a.value), float(b.value)


# ======================================================================
#  lats/lons (real*4)
# ======================================================================

_lib.py_set_lats.argtypes = [C.c_int, C.POINTER(C.c_float)]
_lib.py_set_lons.argtypes = [C.c_int, C.POINTER(C.c_float)]
_lib.py_get_lats.argtypes = [C.c_int, C.POINTER(C.c_float)]
_lib.py_get_lons.argtypes = [C.c_int, C.POINTER(C.c_float)]


def set_lats(arr) -> None:
    arr = np.asarray(arr, dtype=np.float32, order="C")
    n = arr.size
    _lib.py_set_lats(int(n), arr.ctypes.data_as(C.POINTER(C.c_float)))


def set_lons(arr) -> None:
    arr = np.asarray(arr, dtype=np.float32, order="C")
    n = arr.size
    _lib.py_set_lons(int(n), arr.ctypes.data_as(C.POINTER(C.c_float)))


def get_lats() -> np.ndarray:
    n = get_nlats()
    out = np.empty(n, dtype=np.float32, order="C")
    _lib.py_get_lats(int(n), out.ctypes.data_as(C.POINTER(C.c_float)))
    return out


def get_lons() -> np.ndarray:
    n = get_nlons()
    out = np.empty(n, dtype=np.float32, order="C")
    _lib.py_get_lons(int(n), out.ctypes.data_as(C.POINTER(C.c_float)))
    return out


# ======================================================================
#  maps (real*8), shape (nlons, nlats), Fortran-order preferred
# ======================================================================

_lib.py_set_rmap1.argtypes = [C.c_int, C.c_int, C.POINTER(C.c_double)]
_lib.py_set_rmap2.argtypes = [C.c_int, C.c_int, C.POINTER(C.c_double)]
_lib.py_set_ratemap.argtypes = [C.c_int, C.c_int, C.POINTER(C.c_double)]


def _as_f64_fortran_2d(a, nx, ny):
    a = np.asarray(a, dtype=np.float64, order="F")
    if a.shape != (nx, ny):
        raise ValueError(f"array must have shape {(nx, ny)} (got {a.shape})")
    return a


def set_rmap1(A) -> None:
    nx, ny = get_nlons(), get_nlats()
    A = _as_f64_fortran_2d(A, nx, ny)
    _lib.py_set_rmap1(nx, ny, A.ctypes.data_as(C.POINTER(C.c_double)))


def set_rmap2(A) -> None:
    nx, ny = get_nlons(), get_nlats()
    A = _as_f64_fortran_2d(A, nx, ny)
    _lib.py_set_rmap2(nx, ny, A.ctypes.data_as(C.POINTER(C.c_double)))


def set_ratemap(A) -> None:
    nx, ny = get_nlons(), get_nlats()
    A = _as_f64_fortran_2d(A, nx, ny)
    _lib.py_set_ratemap(nx, ny, A.ctypes.data_as(C.POINTER(C.c_double)))


# ======================================================================
#  rMtmp (real*8, len=3)
# ======================================================================

_lib.py_set_rmtmp.argtypes = [C.POINTER(C.c_double)]
_lib.py_get_rmtmp.argtypes = [C.POINTER(C.c_double)]


def set_rmtmp(v3) -> None:
    v = np.asarray(v3, dtype=np.float64, order="C")
    if v.size != 3:
        raise ValueError("rMtmp expects length-3 vector.")
    _lib.py_set_rmtmp(v.ctypes.data_as(C.POINTER(C.c_double)))


def get_rmtmp() -> np.ndarray:
    out = np.empty(3, dtype=np.float64)
    _lib.py_get_rmtmp(out.ctypes.data_as(C.POINTER(C.c_double)))
    return out
