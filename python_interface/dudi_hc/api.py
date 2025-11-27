from __future__ import annotations

"""
Thin, stable Python API for DUDI-hc.

This module maps your typed Python dataclasses (Point, Source, Comet)
to the flat, C-compatible arguments expected by the ctypes bridge.

Scalar functions
----------------
v_integration(point, source, comet, muR, tnow, Rast_AU, pericenter) -> float
delta_ejection(point, source, comet, muR, dt, Rast_AU) -> float
simple_expansion(point, source, cloudcentr, dt) -> float

Batched functions
-----------------
batch_over_points(points, source, comet, muR, tnow, dt, Rast_AU,
                  pericenter, cloudcentr, method) -> ndarray
batch_over_sources(point, sources, comet, muR, tnow, dt, Rast_AU,
                   pericenter, cloudcentr, method) -> ndarray
"""

from typing import Iterable, Sequence
import numpy as np

from .models import Point, Source, Comet
from ._bridge_ctypes import (
    call_v_integration as _call_v_integration,
    call_delta_ejection as _call_delta_ejection,
    call_simple_expansion as _call_simple_expansion,
    call_batch_points as _call_batch_points,
    call_batch_sources as _call_batch_sources,
    METHOD_SIMPLE_EXPANSION,
    METHOD_DELTA_EJECTION,
    METHOD_V_INTEGRATION,
)
from ._bridge_ctypes import (
    get_nlats,
    get_nlons,
    set_lon_bounds,
    get_lon_bounds,
    set_lats,
    set_lons,
    get_lats,
    get_lons,
    set_rmap1,
    set_rmap2,
    set_ratemap,
    set_rmtmp,
    get_rmtmp,
)


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------

def _vec3(x: Iterable[float]) -> np.ndarray:
    """Return x as contiguous float64 vector of shape (3,)."""
    arr = np.asarray(x, dtype=np.float64)
    if arr.shape != (3,):
        raise ValueError(f"Expected a 3-vector (shape (3,)), got shape {arr.shape}")
    return np.ascontiguousarray(arr)


def _method_to_id(method: int | str) -> int:
    """
    Map a user-facing method specifier to the integer ID used in Fortran.
    Accepts either an int (1/2/3) or a string:
      'simple_expansion', 'delta_ejection', 'v_integration'.
    """
    if isinstance(method, int):
        if method in (METHOD_SIMPLE_EXPANSION, METHOD_DELTA_EJECTION, METHOD_V_INTEGRATION):
            return method
        raise ValueError(f"Unknown method id: {method}")

    name = str(method).strip().lower()
    mapping = {
        "simple_expansion": METHOD_SIMPLE_EXPANSION,
        "simple": METHOD_SIMPLE_EXPANSION,
        "delta_ejection": METHOD_DELTA_EJECTION,
        "delta": METHOD_DELTA_EJECTION,
        "v_integration": METHOD_V_INTEGRATION,
        "vintegration": METHOD_V_INTEGRATION,
        "v-int": METHOD_V_INTEGRATION,
    }
    try:
        return mapping[name]
    except KeyError:
        raise ValueError(
            f"Unknown method '{method}'. "
            f"Expected one of: {', '.join(sorted(mapping.keys()))}"
        ) from None


# ----------------------------------------------------------------------
# scalar API
# ----------------------------------------------------------------------

def v_integration(
    point: Point,
    source: Source,
    comet: Comet,
    muR: float,
    tnow: float,
    Rast_AU: float,
    pericenter: bool,
) -> float:
    """
    Compute dust density via velocity integration at a single point.
    """
    # Validation: pericenter must be a real boolean (avoid auto-casting ints)
    if not isinstance(pericenter, bool):
        raise ValueError("pericenter must be a bool.")
    return _call_v_integration(
        point_r=float(point.r),
        point_alpha=float(point.alpha),
        point_beta=float(point.beta),
        point_rvector=_vec3(point.rvector),
        src_r=float(source.r),
        src_alphaM=float(source.alphaM),
        src_betaM=float(source.betaM),
        src_rrM=_vec3(source.rrM),
        src_zeta=float(source.zeta),
        src_eta=float(source.eta),
        src_axis=_vec3(source.symmetry_axis),
        src_eject_distr=int(source.ejection_angle_distr),
        src_ud_shape=int(source.ud.ud_shape),
        src_umin=float(source.ud.umin),
        src_umax=float(source.ud.umax),
        src_Nparticles=float(source.Nparticles),
        src_Tj=float(source.Tj),
        src_dtau=float(source.dtau),
        comet_coords=_vec3(comet.coords),
        comet_vastvec=_vec3(comet.Vastvec),
        comet_vast=float(comet.Vast),
        muR=float(muR),
        tnow=float(tnow),
        Rast_AU=float(Rast_AU),
        pericenter=bool(pericenter),
    )


def delta_ejection(
    point: Point,
    source: Source,
    comet: Comet,
    muR: float,
    dt: float,
    Rast_AU: float,
) -> float:
    """
    Compute density for a delta-function ejection at time offset dt, at one point.
    """
    # Validation: dt must be non-negative
    if float(dt) < 0.0:
        raise ValueError("dt must be >= 0.")
    return _call_delta_ejection(
        point_r=float(point.r),
        point_alpha=float(point.alpha),
        point_beta=float(point.beta),
        point_rvector=_vec3(point.rvector),
        src_r=float(source.r),
        src_alphaM=float(source.alphaM),
        src_betaM=float(source.betaM),
        src_rrM=_vec3(source.rrM),
        src_zeta=float(source.zeta),
        src_eta=float(source.eta),
        src_axis=_vec3(source.symmetry_axis),
        src_eject_distr=int(source.ejection_angle_distr),
        src_ud_shape=int(source.ud.ud_shape),
        src_umin=float(source.ud.umin),
        src_umax=float(source.ud.umax),
        src_Nparticles=float(source.Nparticles),
        src_Tj=float(source.Tj),
        src_dtau=float(source.dtau),
        comet_coords=_vec3(comet.coords),
        comet_vastvec=_vec3(comet.Vastvec),
        comet_vast=float(comet.Vast),
        muR=float(muR),
        dt=float(dt),
        Rast_AU=float(Rast_AU),
    )


def simple_expansion(
    point: Point,
    source: Source,
    cloudcentr: Iterable[float],
    dt: float,
) -> float:
    """
    Compute density in a simple expanding cloud centered at `cloudcentr`, at one point.
    """
    return _call_simple_expansion(
        point_r=float(point.r),
        point_alpha=float(point.alpha),
        point_beta=float(point.beta),
        point_rvector=_vec3(point.rvector),
        src_r=float(source.r),
        src_alphaM=float(source.alphaM),
        src_betaM=float(source.betaM),
        src_rrM=_vec3(source.rrM),
        src_zeta=float(source.zeta),
        src_eta=float(source.eta),
        src_axis=_vec3(source.symmetry_axis),
        src_eject_distr=int(source.ejection_angle_distr),
        src_ud_shape=int(source.ud.ud_shape),
        src_umin=float(source.ud.umin),
        src_umax=float(source.ud.umax),
        src_Nparticles=float(source.Nparticles),
        src_Tj=float(source.Tj),
        src_dtau=float(source.dtau),
        cloudcentr=_vec3(cloudcentr),
        dt=float(dt),
    )


# ----------------------------------------------------------------------
# batched API: over points / over sources
# ----------------------------------------------------------------------

def batch_over_points(
    points: Sequence[Point],
    source: Source,
    comet: Comet,
    muR: float,
    tnow: float,
    dt: float,
    Rast_AU: float,
    pericenter: bool,
    cloudcentr: Iterable[float],
    method: int | str = "v_integration",
) -> np.ndarray:
    """
    Compute densities for many points and a single source.

    Parameters
    ----------
    points : sequence of Point
        Points at which density is evaluated.
    source : Source
        Single dust source.
    comet : Comet
        Ephemeris snapshot.
    muR, tnow, dt, Rast_AU, pericenter, cloudcentr :
        Same meaning as in the scalar Fortran routines.
    method : {'v_integration','delta_ejection','simple_expansion'} or int
        Which DUDI method to use.

    Returns
    -------
    densities : ndarray, shape (N,)
        Density at each point.
    """
    pts = list(points)
    if len(pts) == 0:
        return np.empty(0, dtype=np.float64)

    # build point arrays
    point_r = np.array([float(p.r) for p in pts], dtype=np.float64)
    point_alpha = np.array([float(p.alpha) for p in pts], dtype=np.float64)
    point_beta = np.array([float(p.beta) for p in pts], dtype=np.float64)
    point_rvec = np.vstack([_vec3(p.rvector) for p in pts])  # (N, 3)

    method_id = _method_to_id(method)
    if not isinstance(pericenter, bool):
        raise ValueError("pericenter must be a bool.")

    densities = _call_batch_points(
        point_r=point_r,
        point_alpha=point_alpha,
        point_beta=point_beta,
        point_rvector=point_rvec,
        src_r=float(source.r),
        src_alphaM=float(source.alphaM),
        src_betaM=float(source.betaM),
        src_rrM=_vec3(source.rrM),
        src_zeta=float(source.zeta),
        src_eta=float(source.eta),
        src_axis=_vec3(source.symmetry_axis),
        src_eject_distr=int(source.ejection_angle_distr),
        src_ud_shape=int(source.ud.ud_shape),
        src_umin=float(source.ud.umin),
        src_umax=float(source.ud.umax),
        src_Nparticles=float(source.Nparticles),
        src_Tj=float(source.Tj),
        src_dtau=float(source.dtau),
        comet_coords=_vec3(comet.coords),
        comet_vastvec=_vec3(comet.Vastvec),
        comet_vast=float(comet.Vast),
        muR=float(muR),
        tnow=float(tnow),
        dt=float(dt),
        Rast_AU=float(Rast_AU),
        pericenter=bool(pericenter),
        cloudcentr=_vec3(cloudcentr),
        method_id=method_id,
    )
    return densities


def batch_over_sources(
    point: Point,
    sources: Sequence[Source],
    comet: Comet,
    muR: float,
    tnow: float,
    dt: float,
    Rast_AU: float,
    pericenter: bool,
    cloudcentr: Iterable[float],
    method: int | str = "v_integration",
) -> np.ndarray:
    """
    Compute densities for many sources and a single point.

    Parameters
    ----------
    point : Point
        Observation point.
    sources : sequence of Source
        Multiple dust sources.
    comet : Comet
        Ephemeris snapshot.
    muR, tnow, dt, Rast_AU, pericenter, cloudcentr :
        Same meaning as in the scalar Fortran routines.
    method : {'v_integration','delta_ejection','simple_expansion'} or int
        Which DUDI method to use.

    Returns
    -------
    densities : ndarray, shape (N,)
        Density contribution from each source.
    """
    srcs = list(sources)
    if len(srcs) == 0:
        return np.empty(0, dtype=np.float64)

    n = len(srcs)

    # scalar point
    point_r = float(point.r)
    point_alpha = float(point.alpha)
    point_beta = float(point.beta)
    point_rvec = _vec3(point.rvector)

    # build source arrays
    src_r = np.array([float(s.r) for s in srcs], dtype=np.float64)
    src_alphaM = np.array([float(s.alphaM) for s in srcs], dtype=np.float64)
    src_betaM = np.array([float(s.betaM) for s in srcs], dtype=np.float64)
    src_rrM = np.vstack([_vec3(s.rrM) for s in srcs])          # (N, 3)
    src_zeta = np.array([float(s.zeta) for s in srcs], dtype=np.float64)
    src_eta = np.array([float(s.eta) for s in srcs], dtype=np.float64)
    src_axis = np.vstack([_vec3(s.symmetry_axis) for s in srcs])  # (N, 3)
    src_eject_distr = np.array(
        [int(s.ejection_angle_distr) for s in srcs], dtype=np.int32
    )
    src_ud_shape = np.array(
        [int(s.ud.ud_shape) for s in srcs], dtype=np.int32
    )
    src_umin = np.array([float(s.ud.umin) for s in srcs], dtype=np.float64)
    src_umax = np.array([float(s.ud.umax) for s in srcs], dtype=np.float64)
    src_Nparticles = np.array(
        [float(s.Nparticles) for s in srcs], dtype=np.float64
    )
    src_Tj = np.array([float(s.Tj) for s in srcs], dtype=np.float64)
    src_dtau = np.array([float(s.dtau) for s in srcs], dtype=np.float64)

    method_id = _method_to_id(method)
    if not isinstance(pericenter, bool):
        raise ValueError("pericenter must be a bool.")

    densities = _call_batch_sources(
        point_r=point_r,
        point_alpha=point_alpha,
        point_beta=point_beta,
        point_rvector=point_rvec,
        src_r=src_r,
        src_alphaM=src_alphaM,
        src_betaM=src_betaM,
        src_rrM=src_rrM,
        src_zeta=src_zeta,
        src_eta=src_eta,
        src_axis=src_axis,
        src_eject_distr=src_eject_distr,
        src_ud_shape=src_ud_shape,
        src_umin=src_umin,
        src_umax=src_umax,
        src_Nparticles=src_Nparticles,
        src_Tj=src_Tj,
        src_dtau=src_dtau,
        comet_coords=_vec3(comet.coords),
        comet_vastvec=_vec3(comet.Vastvec),
        comet_vast=float(comet.Vast),
        muR=float(muR),
        tnow=float(tnow),
        dt=float(dt),
        Rast_AU=float(Rast_AU),
        pericenter=bool(pericenter),
        cloudcentr=_vec3(cloudcentr),
        method_id=method_id,
    )
    return densities
