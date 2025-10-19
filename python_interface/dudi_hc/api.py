# python_interface/dudi_hc/api.py
from __future__ import annotations

"""
Thin, stable Python API for DUDI-hc.

This module maps your typed Python dataclasses (Point, Source, Comet)
to the flat, C-compatible arguments expected by the ctypes bridge,
and returns plain Python floats.

Functions
---------
v_integration(point, source, comet, muR, tnow, Rast_AU, pericenter) -> float
delta_ejection(point, source, comet, muR, dt, Rast_AU) -> float
simple_expansion(point, source, cloudcentr, dt) -> float
"""

from typing import Iterable
import numpy as np

from .models import Point, Source, Comet
from ._bridge_ctypes import (
    call_v_integration as _call_v_integration,
    call_delta_ejection as _call_delta_ejection,
    call_simple_expansion as _call_simple_expansion,
)


def _vec3(x: Iterable[float]) -> np.ndarray:
    """Return x as contiguous float64 vector of shape (3,)."""
    arr = np.asarray(x, dtype=np.float64)
    if arr.shape != (3,):
        raise ValueError(f"Expected a 3-vector (shape (3,)), got shape {arr.shape}")
    return np.ascontiguousarray(arr)


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
    Compute dust density via velocity integration at a point.

    Parameters
    ----------
    point : Point
        Observer point (spherical r, alpha, beta + cartesian rvector).
    source : Source
        Dust source parameters (geometry, ejection distribution & speeds).
    comet : Comet
        Ephemeris snapshot (coords, Vastvec, Vast).
    muR : float
        Reduced gravitational parameter (code units, matches Fortran).
    tnow : float
        Epoch offset used by the HC kernel (matches Fortran).
    Rast_AU : float
        Scaling (AU) used by the kernel.
    pericenter : bool
        Whether pericenter branch is used inside HC kernel.

    Returns
    -------
    float
        Density at the point (returned as Python float).
    """
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
    Compute density for a delta-function ejection at time offset dt.

    Returns
    -------
    float
        Density at the point (Python float).
    """
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
    Compute density in a simple expanding cloud centered at `cloudcentr`.

    Parameters
    ----------
    cloudcentr : (3,) array-like
        Cloud center vector (float64, shape (3,)).

    Returns
    -------
    float
        Density at the point (Python float).
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
