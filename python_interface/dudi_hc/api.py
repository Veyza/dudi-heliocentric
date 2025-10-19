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
    Examples
    --------
    >>> import numpy as np
    >>> from python_interface.dudi_hc.models import Point, Source, Comet, EjectionSpeedProperties, spherical_to_cartesian, normalize
    >>> from python_interface.dudi_hc import api
    >>> p = Point(1.0, 1.0, 0.5, spherical_to_cartesian(1.0, 1.0, 0.5))
    >>> s = Source(
    ...     r=1.0, alphaM=1.0, betaM=0.0,
    ...     rrM=spherical_to_cartesian(1.0, 1.0, 0.0),
    ...     zeta=0.3, eta=1.2,
    ...     symmetry_axis=normalize(np.array([0.1, 0.2, 0.97], float)),
    ...     ejection_angle_distr=3,
    ...     ud=EjectionSpeedProperties(ud_shape=1, umin=0.0, umax=0.01),
    ... )
    >>> Vastvec = np.array([0.0001, 0.0004, 0.00003], float)
    >>> c = Comet(coords=np.array([1.0, 0.0, 0.0], float), Vastvec=Vastvec, Vast=float(np.linalg.norm(Vastvec)))
    >>> api.v_integration(p, s, c, muR=0.6, tnow=0.0, Rast_AU=0.0, pericenter=False)
    Traceback (most recent call last):
        ...
    NotImplementedError: Fortran bridge not wired yet: hc_DUDI_v_integration
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
