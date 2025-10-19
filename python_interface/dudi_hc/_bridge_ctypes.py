# python_interface/dudi_hc/_bridge_ctypes.py
from __future__ import annotations
import ctypes as C
from pathlib import Path
import numpy as np

# Load the shared library placed by build_ctypes_bridge.sh
_lib_path = Path(__file__).with_name("libpy_dudihc_bridge.so")
if not _lib_path.exists():
    raise ImportError(f"Shared library not found: {_lib_path}. "
                      f"Run: bash python_interface/fortran_bridge/build_ctypes_bridge.sh")
_lib = C.CDLL(str(_lib_path))

# Helpers
Vec3d = np.ctypeslib.ndpointer(dtype=np.float64, shape=(3,), flags="C_CONTIGUOUS")

def _as_vec3(a) -> np.ndarray:
    arr = np.asarray(a, dtype=np.float64)
    if arr.shape != (3,):
        raise ValueError(f"expected shape (3,), got {arr.shape}")
    return np.ascontiguousarray(arr)

def _c_int(b: bool | int) -> int:
    return 1 if bool(b) else 0

# ---- prototypes ----

# void py_hc_v_integration(double r, double alpha, double beta, double rvec[3],
#   double sr, double sAlphaM, double sBetaM, double sRR[3], double zeta, double eta,
#   double axis[3], int eject_distr, int ud_shape, double umin, double umax,
#   double comet_coords[3], double comet_vvec[3], double comet_v,
#   double muR, double tnow, double Rast_AU, int pericenter, double* density_out)
_lib.py_hc_v_integration.argtypes = [
    C.c_double, C.c_double, C.c_double, Vec3d,
    C.c_double, C.c_double, C.c_double, Vec3d, C.c_double, C.c_double,
    Vec3d, C.c_int, C.c_int, C.c_double, C.c_double,
    Vec3d, Vec3d, C.c_double,
    C.c_double, C.c_double, C.c_double, C.c_int, C.POINTER(C.c_double)
]
_lib.py_hc_v_integration.restype = None

_lib.py_hc_delta_ejection.argtypes = [
    C.c_double, C.c_double, C.c_double, Vec3d,
    C.c_double, C.c_double, C.c_double, Vec3d, C.c_double, C.c_double,
    Vec3d, C.c_int, C.c_int, C.c_double, C.c_double,
    Vec3d, Vec3d, C.c_double,
    C.c_double, C.c_double, C.c_double, C.POINTER(C.c_double)
]
_lib.py_hc_delta_ejection.restype = None

_lib.py_hc_simple_expansion.argtypes = [
    C.c_double, C.c_double, C.c_double, Vec3d,
    C.c_double, C.c_double, C.c_double, Vec3d, C.c_double, C.c_double,
    Vec3d, C.c_int, C.c_int, C.c_double, C.c_double,
    Vec3d, C.c_double, C.POINTER(C.c_double)
]
_lib.py_hc_simple_expansion.restype = None

# ---- thin Python functions returning float ----

def call_v_integration(
    *,
    point_r: float, point_alpha: float, point_beta: float, point_rvector,
    src_r: float, src_alphaM: float, src_betaM: float, src_rrM,
    src_zeta: float, src_eta: float, src_axis, src_eject_distr: int,
    src_ud_shape: int, src_umin: float, src_umax: float,
    comet_coords, comet_vastvec, comet_vast: float,
    muR: float, tnow: float, Rast_AU: float, pericenter: bool
) -> float:
    out = C.c_double()
    _lib.py_hc_v_integration(
        float(point_r), float(point_alpha), float(point_beta), _as_vec3(point_rvector),
        float(src_r), float(src_alphaM), float(src_betaM), _as_vec3(src_rrM),
        float(src_zeta), float(src_eta), _as_vec3(src_axis),
        int(src_eject_distr), int(src_ud_shape), float(src_umin), float(src_umax),
        _as_vec3(comet_coords), _as_vec3(comet_vastvec), float(comet_vast),
        float(muR), float(tnow), float(Rast_AU), _c_int(pericenter),
        C.byref(out)
    )
    return float(out.value)

def call_delta_ejection(
    *,
    point_r: float, point_alpha: float, point_beta: float, point_rvector,
    src_r: float, src_alphaM: float, src_betaM: float, src_rrM,
    src_zeta: float, src_eta: float, src_axis, src_eject_distr: int,
    src_ud_shape: int, src_umin: float, src_umax: float,
    comet_coords, comet_vastvec, comet_vast: float,
    muR: float, dt: float, Rast_AU: float
) -> float:
    out = C.c_double()
    _lib.py_hc_delta_ejection(
        float(point_r), float(point_alpha), float(point_beta), _as_vec3(point_rvector),
        float(src_r), float(src_alphaM), float(src_betaM), _as_vec3(src_rrM),
        float(src_zeta), float(src_eta), _as_vec3(src_axis),
        int(src_eject_distr), int(src_ud_shape), float(src_umin), float(src_umax),
        _as_vec3(comet_coords), _as_vec3(comet_vastvec), float(comet_vast),
        float(muR), float(dt), float(Rast_AU), C.byref(out)
    )
    return float(out.value)

def call_simple_expansion(
    *,
    point_r: float, point_alpha: float, point_beta: float, point_rvector,
    src_r: float, src_alphaM: float, src_betaM: float, src_rrM,
    src_zeta: float, src_eta: float, src_axis, src_eject_distr: int,
    src_ud_shape: int, src_umin: float, src_umax: float,
    cloudcentr, dt: float
) -> float:
    out = C.c_double()
    _lib.py_hc_simple_expansion(
        float(point_r), float(point_alpha), float(point_beta), _as_vec3(point_rvector),
        float(src_r), float(src_alphaM), float(src_betaM), _as_vec3(src_rrM),
        float(src_zeta), float(src_eta), _as_vec3(src_axis),
        int(src_eject_distr), int(src_ud_shape), float(src_umin), float(src_umax),
        _as_vec3(cloudcentr), float(dt), C.byref(out)
    )
    return float(out.value)
