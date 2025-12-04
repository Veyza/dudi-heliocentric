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
        float(muR), float(tnow), float(Rast_AU), C.c_int(pericenter),
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

# void py_hc_batch_sources_points(
#   int n_points, int Nt, int Ns,
#   double density[n_points],
#   double point_r[n_points], point_alpha[n_points], point_beta[n_points],
#   double point_rvector[3*n_points],
#   double src_r[Nt*Ns], src_alphaM[Nt*Ns], src_betaM[Nt*Ns],
#   double src_rrM[3*Nt*Ns], double src_zeta[Nt*Ns], double src_eta[Nt*Ns],
#   double src_axis[3*Nt*Ns],
#   int src_eject_distr[Nt*Ns], int src_ud_shape[Nt*Ns],
#   double src_umin[Nt*Ns], double src_umax[Nt*Ns],
#   double src_Nparticles[Nt*Ns], double src_Tj[Nt*Ns], double src_dtau[Nt*Ns],
#   double comet_coords[3*Nt], double comet_vvec[3*Nt], double comet_vast[Nt],
#   double muR, double tnow, double Rast_AU, int pericenter, int method_id)
_lib.py_hc_batch_sources_points.argtypes = [
    C.c_int,  # n_points
    C.c_int,  # Nt
    C.c_int,  # Ns
    Vec1d,    # density_out
    # points
    Vec1d, Vec1d, Vec1d,  # point_r, alpha, beta
    Vec1d,                # point_rvector_flat (3*n_points)
    # sources
    Vec1d, Vec1d, Vec1d,  # src_r, src_alphaM, src_betaM
    Vec1d, Vec1d, Vec1d,  # src_rrM_flat, src_zeta, src_eta
    Vec1d, Int1d, Int1d,  # src_axis_flat, eject_distr, ud_shape
    Vec1d, Vec1d,         # src_umin, src_umax
    Vec1d, Vec1d, Vec1d,  # Nparticles, Tj, dtau
    # comets
    Vec1d, Vec1d, Vec1d,  # comet_coords_flat, comet_vvec_flat, comet_vast
    # scalars
    C.c_double, C.c_double, C.c_double, C.c_int, C.c_int,
]
_lib.py_hc_batch_sources_points.restype = None




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
        float(muR), float(tnow), float(dt), float(Rast_AU), C.c_int(pericenter),
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
        float(muR), float(tnow), float(dt), float(Rast_AU), C.c_int(pericenter),
        _as_vec3(cloudcentr), int(method_id),
    )
    return density

def call_batch_sources_points(
    *,
    point_r,
    point_alpha,
    point_beta,
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
    comet_vast,
    muR: float,
    tnow: float,
    Rast_AU: float,
    pericenter: bool,
    method_id: int,
) -> np.ndarray:
    """
    Low-level wrapper for hc_DUDI_batch_sources_points (time × sources × points).

    Parameters
    ----------
    point_* : arrays over points (length n_points)
    src_*   : arrays over flattened (Nt, Ns) with C-order:
              idx = it * Ns + is
    src_rrM, src_axis : arrays of length 3*Nt*Ns, similarly flattened
    comet_coords, comet_vastvec : arrays of shape (Nt, 3)
    comet_vast : array of length Nt

    Returns
    -------
    density : ndarray, shape (n_points,)
        Total density per point (sum over times and sources).
    """
    # points
    r = _as_1d_f64(point_r, "point_r")
    alpha = _as_1d_f64(point_alpha, "point_alpha")
    beta = _as_1d_f64(point_beta, "point_beta")
    if not (r.size == alpha.size == beta.size):
        raise ValueError("point_r, point_alpha, point_beta must have same length")
    n_points = int(r.size)

    rvec2d = _as_2d_f64(point_rvector, "point_rvector")
    if rvec2d.shape != (n_points, 3):
        raise ValueError(f"point_rvector must have shape (N_points,3), got {rvec2d.shape}")
    rvec_flat = np.ascontiguousarray(rvec2d.reshape(-1))

    # sources: src_* expected as 2D (Nt, Ns) or 3D (Nt, Ns, 3)
    src_r_2d      = _as_2d_f64(src_r, "src_r")
    src_alpha_2d  = _as_2d_f64(src_alphaM, "src_alphaM")
    src_beta_2d   = _as_2d_f64(src_betaM, "src_betaM")
    Nt, Ns = src_r_2d.shape
    if src_alpha_2d.shape != (Nt, Ns) or src_beta_2d.shape != (Nt, Ns):
        raise ValueError("src_r, src_alphaM, src_betaM must all have shape (Nt,Ns)")

    src_rrM_3d = np.asarray(src_rrM, dtype=np.float64)
    if src_rrM_3d.ndim != 3 or src_rrM_3d.shape != (Nt, Ns, 3):
        raise ValueError(f"src_rrM must have shape (Nt,Ns,3), got {src_rrM_3d.shape}")
    src_axis_3d = np.asarray(src_axis, dtype=np.float64)
    if src_axis_3d.ndim != 3 or src_axis_3d.shape != (Nt, Ns, 3):
        raise ValueError(f"src_axis must have shape (Nt,Ns,3), got {src_axis_3d.shape}")

    src_zeta_2d = _as_2d_f64(src_zeta, "src_zeta")
    src_eta_2d  = _as_2d_f64(src_eta, "src_eta")
    if src_zeta_2d.shape != (Nt, Ns) or src_eta_2d.shape != (Nt, Ns):
        raise ValueError("src_zeta and src_eta must have shape (Nt,Ns)")

    eject = _as_1d_i32(src_eject_distr, "src_eject_distr")
    udsh  = _as_1d_i32(src_ud_shape, "src_ud_shape")
    if eject.size != Nt * Ns or udsh.size != Nt * Ns:
        raise ValueError("src_eject_distr and src_ud_shape must have length Nt*Ns")

    src_umin_2d  = _as_2d_f64(src_umin, "src_umin")
    src_umax_2d  = _as_2d_f64(src_umax, "src_umax")
    src_Np_2d    = _as_2d_f64(src_Nparticles, "src_Nparticles")
    src_Tj_2d    = _as_2d_f64(src_Tj, "src_Tj")
    src_dtau_2d  = _as_2d_f64(src_dtau, "src_dtau")
    for name, arr2d in [
        ("src_umin", src_umin_2d),
        ("src_umax", src_umax_2d),
        ("src_Nparticles", src_Np_2d),
        ("src_Tj", src_Tj_2d),
        ("src_dtau", src_dtau_2d),
    ]:
        if arr2d.shape != (Nt, Ns):
            raise ValueError(f"{name} must have shape (Nt,Ns)")

    # flatten with C-order -> idx = it*Ns + is
    src_r_flat      = np.ascontiguousarray(src_r_2d.reshape(-1))
    src_alpha_flat  = np.ascontiguousarray(src_alpha_2d.reshape(-1))
    src_beta_flat   = np.ascontiguousarray(src_beta_2d.reshape(-1))
    src_rrM_flat    = np.ascontiguousarray(src_rrM_3d.reshape(-1))
    src_zeta_flat   = np.ascontiguousarray(src_zeta_2d.reshape(-1))
    src_eta_flat    = np.ascontiguousarray(src_eta_2d.reshape(-1))
    src_axis_flat   = np.ascontiguousarray(src_axis_3d.reshape(-1))
    src_umin_flat   = np.ascontiguousarray(src_umin_2d.reshape(-1))
    src_umax_flat   = np.ascontiguousarray(src_umax_2d.reshape(-1))
    src_Np_flat     = np.ascontiguousarray(src_Np_2d.reshape(-1))
    src_Tj_flat     = np.ascontiguousarray(src_Tj_2d.reshape(-1))
    src_dtau_flat   = np.ascontiguousarray(src_dtau_2d.reshape(-1))

    # comets: coords, vvec: (Nt,3), vast: (Nt,)
    comet_coords_2d = _as_2d_f64(comet_coords, "comet_coords")
    comet_vvec_2d   = _as_2d_f64(comet_vastvec, "comet_vastvec")
    if comet_coords_2d.shape != (Nt, 3) or comet_vvec_2d.shape != (Nt, 3):
        raise ValueError("comet_coords and comet_vastvec must have shape (Nt,3)")
    comet_coords_flat = np.ascontiguousarray(comet_coords_2d.reshape(-1))
    comet_vvec_flat   = np.ascontiguousarray(comet_vvec_2d.reshape(-1))

    comet_vast_arr = _as_1d_f64(comet_vast, "comet_vast")
    if comet_vast_arr.size != Nt:
        raise ValueError("comet_vast must have length Nt")

    density = np.empty(n_points, dtype=np.float64)

    _lib.py_hc_batch_sources_points(
        int(n_points),
        int(Nt),
        int(Ns),
        density,
        r, alpha, beta,
        rvec_flat,
        src_r_flat, src_alpha_flat, src_beta_flat,
        src_rrM_flat, src_zeta_flat, src_eta_flat,
        src_axis_flat, eject, udsh,
        src_umin_flat, src_umax_flat,
        src_Np_flat, src_Tj_flat, src_dtau_flat,
        comet_coords_flat, comet_vvec_flat, comet_vast_arr,
        float(muR), float(tnow), float(Rast_AU),
        C.c_int(pericenter),
        int(method_id),
    )
    return density



# ======================================================================
# Ratemap / impact map interface to Fortran (low-level ctypes layer)
# ======================================================================

# --- file-reading wrappers --------------------------------------------

# void py_read_ratemap_get_rhel(const char *fname, double *rhel);
_lib.py_read_ratemap_get_rhel.argtypes = [C.c_char_p, C.POINTER(C.c_double)]
_lib.py_read_ratemap_get_rhel.restype = None


def call_read_ratemap_get_rhel(filename: str) -> float:
    """
    Call Fortran py_read_ratemap_get_rhel(fname, rhel) and return rhel.
    """
    rhel = C.c_double()
    fname_bytes = filename.encode("utf-8")
    _lib.py_read_ratemap_get_rhel(fname_bytes, C.byref(rhel))
    return float(rhel.value)


# void py_ratematr_interpolate(double rhel, double rhel1, double rhel2);
_lib.py_ratematr_interpolate.argtypes = [C.c_double, C.c_double, C.c_double]
_lib.py_ratematr_interpolate.restype = None


def call_ratematr_interpolate(rhel: float, rhel1: float, rhel2: float) -> None:
    """
    Call Fortran py_ratematr_interpolate(rhel, rhel1, rhel2).
    """
    _lib.py_ratematr_interpolate(
        C.c_double(rhel), C.c_double(rhel1), C.c_double(rhel2)
    )


# --- dimensions & longitude limits -------------------------------------

# void py_get_ratemap_dims(int *nlats, int *nlons);
_lib.py_get_ratemap_dims.argtypes = [C.POINTER(C.c_int), C.POINTER(C.c_int)]
_lib.py_get_ratemap_dims.restype = None


def call_get_ratemap_dims() -> tuple[int, int]:
    """
    Return (nlats, nlons) from Fortran.
    """
    nlats = C.c_int()
    nlons = C.c_int()
    _lib.py_get_ratemap_dims(C.byref(nlats), C.byref(nlons))
    return int(nlats.value), int(nlons.value)

# void py_get_lon_limits(double *lonmin, double *lonmax);
_lib.py_get_lon_limits.argtypes = [C.POINTER(C.c_double), C.POINTER(C.c_double)]
_lib.py_get_lon_limits.restype = None


def call_get_lon_limits() -> tuple[float, float]:
    lonmin = C.c_double()
    lonmax = C.c_double()
    _lib.py_get_lon_limits(C.byref(lonmin), C.byref(lonmax))
    return float(lonmin.value), float(lonmax.value)


# void py_set_lon_limits(const double *lonmin, const double *lonmax);
_lib.py_set_lon_limits.argtypes = [C.POINTER(C.c_double), C.POINTER(C.c_double)]
_lib.py_set_lon_limits.restype = None


def call_set_lon_limits(lonmin: float, lonmax: float) -> None:
    lonmin_c = C.c_double(lonmin)
    lonmax_c = C.c_double(lonmax)
    _lib.py_set_lon_limits(C.byref(lonmin_c), C.byref(lonmax_c))


# --- 1D grids: lats, lons ----------------------------------------------

# void py_get_lats(double *lats_out);  ! expects length nlats
_lib.py_get_lats.argtypes = [C.POINTER(C.c_float)]
_lib.py_get_lats.restype = None

# void py_set_lats(const double *lats_in);  ! expects length nlats
_lib.py_set_lats.argtypes = [C.POINTER(C.c_float)]
_lib.py_set_lats.restype = None


def call_get_lats(nlats: int) -> np.ndarray:
    """
    Return lats as a 1D NumPy array of shape (nlats,).
    """
    arr = np.empty(nlats, dtype=np.float32)
    ptr = arr.ctypes.data_as(C.POINTER(C.c_float))
    _lib.py_get_lats(ptr)
    return arr.astype(np.float64)


def call_set_lats(lats: np.ndarray) -> None:
    """
    Copy a 1D NumPy array into Fortran lats.
    """
    lats = np.asarray(lats, dtype=np.float32)
    ptr = lats.ctypes.data_as(C.POINTER(C.c_float))
    _lib.py_set_lats(ptr)


# void py_get_lons(double *lons_out);  ! expects length nlons
_lib.py_get_lons.argtypes = [C.POINTER(C.c_float)]
_lib.py_get_lons.restype = None

# void py_set_lons(const double *lons_in);  ! expects length nlons
_lib.py_set_lons.argtypes = [C.POINTER(C.c_float)]
_lib.py_set_lons.restype = None


def call_get_lons(nlons: int) -> np.ndarray:
    """
    Return lons as a 1D NumPy array of shape (nlons,).
    """
    arr = np.empty(nlons, dtype=np.float32)
    ptr = arr.ctypes.data_as(C.POINTER(C.c_float))
    _lib.py_get_lons(ptr)
    return arr.astype(np.float64)


def call_set_lons(lons: np.ndarray) -> None:
    """
    Copy a 1D NumPy array into Fortran lons.
    """
    lons = np.asarray(lons, dtype=np.float32)
    ptr = lons.ctypes.data_as(C.POINTER(C.c_float))
    _lib.py_set_lons(ptr)


# --- 2D maps: ratemap, rmap1, rmap2 ------------------------------------

# void py_get_ratemap(double *ratemap_out);  ! length nlats*nlons
_lib.py_get_ratemap.argtypes = [C.POINTER(C.c_double)]
_lib.py_get_ratemap.restype = None

# void py_set_ratemap(const double *ratemap_in);
_lib.py_set_ratemap.argtypes = [C.POINTER(C.c_double)]
_lib.py_set_ratemap.restype = None


def call_get_ratemap_flat(nlats: int, nlons: int) -> np.ndarray:
    """
    Return ratemap as a flat 1D array of length nlats*nlons.
    """
    arr = np.empty(nlats * nlons, dtype=np.float64)
    ptr = arr.ctypes.data_as(C.POINTER(C.c_double))
    _lib.py_get_ratemap(ptr)
    return arr


def call_set_ratemap_from_flat(ratemap_flat: np.ndarray) -> None:
    ratemap_flat = np.asarray(ratemap_flat, dtype=np.float64)
    ptr = ratemap_flat.ctypes.data_as(C.POINTER(C.c_double))
    _lib.py_set_ratemap(ptr)


# rmap1
_lib.py_get_rmap1.argtypes = [C.POINTER(C.c_double)]
_lib.py_get_rmap1.restype = None

_lib.py_set_rmap1.argtypes = [C.POINTER(C.c_double)]
_lib.py_set_rmap1.restype = None


def call_get_rmap1_flat(nlats: int, nlons: int) -> np.ndarray:
    arr = np.empty(nlats * nlons, dtype=np.float64)
    ptr = arr.ctypes.data_as(C.POINTER(C.c_double))
    _lib.py_get_rmap1(ptr)
    return arr


def call_set_rmap1_from_flat(rmap1_flat: np.ndarray) -> None:
    rmap1_flat = np.asarray(rmap1_flat, dtype=np.float64)
    ptr = rmap1_flat.ctypes.data_as(C.POINTER(C.c_double))
    _lib.py_set_rmap1(ptr)


# rmap2
_lib.py_get_rmap2.argtypes = [C.POINTER(C.c_double)]
_lib.py_get_rmap2.restype = None

_lib.py_set_rmap2.argtypes = [C.POINTER(C.c_double)]
_lib.py_set_rmap2.restype = None


def call_get_rmap2_flat(nlats: int, nlons: int) -> np.ndarray:
    arr = np.empty(nlats * nlons, dtype=np.float64)
    ptr = arr.ctypes.data_as(C.POINTER(C.c_double))
    _lib.py_get_rmap2(ptr)
    return arr


def call_set_rmap2_from_flat(rmap2_flat: np.ndarray) -> None:
    rmap2_flat = np.asarray(rmap2_flat, dtype=np.float64)
    ptr = rmap2_flat.ctypes.data_as(C.POINTER(C.c_double))
    _lib.py_set_rmap2(ptr)


# --- temporary vector rMtmp (size 3) -----------------------------------

_lib.py_get_rMtmp.argtypes = [C.POINTER(C.c_double)]
_lib.py_get_rMtmp.restype = None

_lib.py_set_rMtmp.argtypes = [C.POINTER(C.c_double)]
_lib.py_set_rMtmp.restype = None


def call_get_rMtmp() -> np.ndarray:
    arr = np.empty(3, dtype=np.float64)
    ptr = arr.ctypes.data_as(C.POINTER(C.c_double))
    _lib.py_get_rMtmp(ptr)
    return arr


def call_set_rMtmp(rMtmp: np.ndarray) -> None:
    rMtmp = np.asarray(rMtmp, dtype=np.float64)
    if rMtmp.size != 3:
        raise ValueError("rMtmp must have length 3")
    ptr = rMtmp.ctypes.data_as(C.POINTER(C.c_double))
    _lib.py_set_rMtmp(ptr)


# ======================================================================
# End of ratemap / impact map ctypes layer
# ======================================================================
