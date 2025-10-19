from __future__ import annotations

"""
Public API for DUDI-heliocentric (thin, stable surface).

This module exposes three functions that mirror the core Fortran entry points.
They validate inputs and define a stable interface for downstream users. The
actual numerical work will be wired in a later step via an f2py bridge.

Mappings
--------
v_integration   -> Fortran: hc_DUDI_v_integration
delta_ejection  -> Fortran: hc_DUDI_delta_ejection
simple_expansion-> Fortran: hc_DUDI_simple_expansion
"""

from typing import Final
import math
import numpy as np

from .typing import Vec3
from .models import Point, Source, Comet, as_vec3

__all__ = ["v_integration", "delta_ejection", "simple_expansion"]


def _check_finite_scalar(x: float, name: str) -> None:
    if not (isinstance(x, (int, float)) and math.isfinite(float(x))):
        raise ValueError(f"{name} must be a finite float.")
def _check_nonneg_scalar(x: float, name: str) -> None:
    _check_finite_scalar(x, name)
    if float(x) < 0.0:
        raise ValueError(f"{name} must be >= 0.")
def _check_bool(b: bool, name: str) -> None:
    if not isinstance(b, (bool, np.bool_)):
        raise ValueError(f"{name} must be a boolean.")


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
    Compute dust number density using the v-integration method.

    Parameters
    ----------
    point : Point
        Position where density is evaluated. Requires `rvector` (AU).
    source : Source
        Dust source definition. Requires `rrM` (AU) and `symmetry_axis` (unit 3-vector).
    comet : Comet
        State of the dust-emitting body at ejection.
    muR : float
        Reduced gravitational parameter [AU^3/day^2].
    tnow : float
        Absolute time at which density is evaluated [day].
    Rast_AU : float
        Radius of the source body [AU]. May be 0 if re-collisions are neglected.
    pericenter : bool
        Whether particles have passed perihelion in their orbit from source to point.

    Returns
    -------
    float
        Number density at `point` (units per the model; typically 1/AU^3).

    Notes
    -----
    This is a thin wrapper for Fortran `hc_DUDI_v_integration`. In this step it
    only validates inputs and raises NotImplementedError. The numerical bridge
    will be added later.
    """
    # --- lightweight validation (cheap & early) ---
    _check_finite_scalar(muR, "muR")
    _check_finite_scalar(tnow, "tnow")
    _check_nonneg_scalar(Rast_AU, "Rast_AU")
    _check_bool(pericenter, "pericenter")

    # Ensure stored vectors are correct shape; do not modify them
    as_vec3(point.rvector, name="point.rvector")
    as_vec3(source.rrM, name="source.rrM")
    as_vec3(source.symmetry_axis, name="source.symmetry_axis")
    as_vec3(comet.coords, name="comet.coords")
    as_vec3(comet.Vastvec, name="comet.Vastvec")

    raise NotImplementedError("Fortran bridge not wired yet: hc_DUDI_v_integration")


def delta_ejection(
    point: Point,
    source: Source,
    comet: Comet,
    muR: float,
    dt: float,
    Rast_AU: float,
) -> float:
    """
    Compute dust number density using the delta-ejection method.

    Parameters
    ----------
    point : Point
        Position where density is evaluated. Requires `rvector` (AU).
    source : Source
        Dust source definition. Requires `rrM` (AU) and `symmetry_axis` (unit 3-vector).
    comet : Comet
        State of the dust-emitting body at ejection.
    muR : float
        Reduced gravitational parameter [AU^3/day^2].
    dt : float
        Time elapsed since dust ejection [day], dt >= 0.
    Rast_AU : float
        Radius of the source body [AU]. May be 0 if re-collisions are neglected.

    Returns
    -------
    float
        Number density at `point` (units per the model; typically 1/AU^3).

    Notes
    -----
    Thin wrapper for Fortran `hc_DUDI_delta_ejection`. Stub for now.
    """
    _check_finite_scalar(muR, "muR")
    _check_nonneg_scalar(dt, "dt")
    _check_nonneg_scalar(Rast_AU, "Rast_AU")

    as_vec3(point.rvector, name="point.rvector")
    as_vec3(source.rrM, name="source.rrM")
    as_vec3(source.symmetry_axis, name="source.symmetry_axis")
    as_vec3(comet.coords, name="comet.coords")
    as_vec3(comet.Vastvec, name="comet.Vastvec")

    raise NotImplementedError("Fortran bridge not wired yet: hc_DUDI_delta_ejection")


def simple_expansion(
    point: Point,
    source: Source,
    cloudcentr: Vec3,
    dt: float,
) -> float:
    """
    Compute dust number density using the simple expansion method.

    Parameters
    ----------
    point : Point
        Position where density is evaluated. Requires `rvector` (AU).
    source : Source
        Dust source definition. Requires `rrM` (AU) and `symmetry_axis` (unit 3-vector).
    cloudcentr : Vec3
        Heliocentric position of the prime cloud center [AU], shape (3,).
    dt : float
        Time elapsed since dust ejection [day], dt >= 0.

    Returns
    -------
    float
        Number density at `point` (units per the model; typically 1/AU^3).

    Notes
    -----
    Thin wrapper for Fortran `hc_DUDI_simple_expansion`. Stub for now.
    """
    _check_nonneg_scalar(dt, "dt")
    as_vec3(point.rvector, name="point.rvector")
    as_vec3(source.rrM, name="source.rrM")
    as_vec3(source.symmetry_axis, name="source.symmetry_axis")
    as_vec3(cloudcentr, name="cloudcentr")

    raise NotImplementedError("Fortran bridge not wired yet: hc_DUDI_simple_expansion")
