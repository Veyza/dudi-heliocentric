"""
Data models for the DUDI-heliocentric Python interface.

These mirror the Fortran derived types:
- Point:    (r, alpha, beta, rvector)
- Source:   (r, alphaM, betaM, rrM, zeta, eta, symmetry_axis, ejection_angle_distr, ud)
- Comet:    (coords, Vastvec, Vast)
Units:
- Distances in AU, times in day, angles in rad, velocities/speeds in AU/day.
- 3-vectors are numpy float64 arrays with shape (3,).
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from .typing import Vec3

def as_vec3(x, name="vector") -> Vec3:
    v = np.asarray(x, dtype=np.float64)
    if v.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {v.shape}.")
    return v

def normalize(v: Vec3, atol: float = 1e-15) -> Vec3:
    v = as_vec3(v, "normalize(v)")
    n = float(np.linalg.norm(v))
    if not math.isfinite(n) or n <= atol:
        raise ValueError("Cannot normalize near-zero vector.")
    return (v / n).astype(np.float64)

def spherical_to_cartesian(r: float, alpha: float, beta: float) -> Vec3:
    import math
    sa, ca = math.sin(alpha), math.cos(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    return np.array([r*sa*cb, r*sa*sb, r*ca], dtype=np.float64)

@dataclass
class EjectionSpeedProperties:
    """Matches Fortran type(ejection_speed_properties).

    Fields
    ------
    ud_shape : int
        Selector for ejection speed PDF.
    umin : float
        Minimum ejection speed [AU/day], umin >= 0.
    umax : float
        Maximum ejection speed [AU/day], umax >= umin.
    """

    """Matches Fortran type(ejection_speed_properties)."""
    ud_shape: int
    umin: float   # AU/day
    umax: float   # AU/day
    def __post_init__(self):
        if not (isinstance(self.ud_shape, int)):
            raise ValueError("ud_shape must be int.")
        if not (math.isfinite(self.umin) and self.umin >= 0.0):
            raise ValueError("umin must be finite and >= 0.")
        if not (math.isfinite(self.umax) and self.umax >= self.umin):
            raise ValueError("umax must be finite and >= umin.")

@dataclass
class Point:
    """
    Location where the dust number density is evaluated.

    Fields
    ------
    r : float
        Heliocentric radial distance [AU], r >= 0.
    alpha : float
        Polar angle [rad].
    beta : float
        Eastern longitude [rad].
    rvector : Vec3
        Cartesian coordinates [AU], numpy array of shape (3,).

    Validation
    ----------
    - r is finite and >= 0
    - alpha, beta are finite
    - rvector has shape (3,) and dtype float64
    Matches Fortran 'point' structure."""
    r: float           # AU
    alpha: float       # rad
    beta: float        # rad
    rvector: Vec3      # AU, stored 3-vector
    def __post_init__(self):
        if not (math.isfinite(self.r) and self.r >= 0.0):
            raise ValueError("Point.r must be finite and >= 0.")
        if not (math.isfinite(self.alpha) and math.isfinite(self.beta)):
            raise ValueError("Point angles must be finite.")
        _ = as_vec3(self.rvector, "Point.rvector")

@dataclass
class Source:
    """
    Dust source definition at ejection.

    Fields
    ------
    r, alphaM, betaM : float
        Source heliocentric spherical coordinates [AU, rad, rad].
    rrM : Vec3
        Source Cartesian position [AU], shape (3,).
    zeta, eta : float
        Orientation angles [rad]. Interpretation depends on your local frame.
        They do not define `symmetry_axis` here; `symmetry_axis` must be provided.
    symmetry_axis : Vec3
        Unit 3-vector (shape (3,)) along the ejection symmetry axis (local frame).
    ejection_angle_distr : int
        Selector for the ejection direction distribution.
    ud : EjectionSpeedProperties
        Ejection speed PDF parameters (shape id, umin, umax) in AU/day.

    Validation
    ----------
    - r >= 0; all angles finite
    - rrM shape (3,)
    - symmetry_axis shape (3,) and norm ~ 1 (within 1e-9)
    - ejection_angle_distr is int
    Matches Fortran 'source' structure (heliocentric)."""
    r: float           # AU
    alphaM: float      # rad
    betaM: float       # rad
    rrM: Vec3          # AU, stored 3-vector
    zeta: float        # rad
    eta: float         # rad
    symmetry_axis: Vec3
    ejection_angle_distr : int
    ud : EjectionSpeedProperties
    Nparticles : float  # particles
    Tj : float          # days
    dtau : float        # days
    def __post_init__(self):
        if not (math.isfinite(self.r) and self.r >= 0.0):
            raise ValueError("Source.r must be finite and >= 0.")
        for name, val in (("alphaM", self.alphaM), ("betaM", self.betaM),
                          ("zeta", self.zeta), ("eta", self.eta)):
            if not math.isfinite(val):
                raise ValueError(f"Source.{name} must be finite.")
        _ = as_vec3(self.rrM, "Source.rrM")
        ax = as_vec3(self.symmetry_axis, "Source.symmetry_axis")
        n = float(np.linalg.norm(ax))
        if not (math.isfinite(n) and abs(n - 1.0) <= 1e-9):
            raise ValueError("Source.symmetry_axis must be unit length.")
        if not isinstance(self.ejection_angle_distr, int):
            raise ValueError("ejection_angle_distr must be int.")

@dataclass
class Comet:
    """
    State of the dust-emitting body at ejection.

    Fields
    ------
    coords : Vec3
        Heliocentric coordinates [AU], shape (3,).
    Vastvec : Vec3
        Heliocentric velocity vector [AU/day], shape (3,).
    Vast : float
        Speed magnitude [AU/day]. Typically equals ||Vastvec||.

    Validation
    ----------
    - coords, Vastvec shape (3,)
    - Vast finite and >= 0
    - We do not enforce Vast == ||Vastvec|| (documented expectation only).
    Matches Fortran 'comet' structure."""
    coords: Vec3     # AU
    Vastvec: Vec3    # AU/day
    Vast: float      # AU/day
    def __post_init__(self):
        _ = as_vec3(self.coords, "Comet.coords")
        v = as_vec3(self.Vastvec, "Comet.Vastvec")
        if not (math.isfinite(self.Vast) and self.Vast >= 0.0):
            raise ValueError("Comet.Vast must be finite and >= 0.")
        # Optional: don't enforce Vast == ||Vastvec||, but document expectation.

def __post_init__(self):
        if not (math.isfinite(self.Nparticles) and self.Nparticles > 0.0):
            raise ValueError("Source.Nparticles must be finite and > 0.")
        if not math.isfinite(self.Tj):
            raise ValueError("Source.Tj must be finite.")
        if not (math.isfinite(self.dtau) and self.dtau >= 0.0):
            raise ValueError("Source.dtau must be finite and >= 0.")