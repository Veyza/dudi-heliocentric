# python_interface/dudi_hc/models.py (essence)

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

@dataclass(frozen=True)
class EjectionSpeedProperties:
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

@dataclass(frozen=True)
class Point:
    """Matches Fortran 'point' structure."""
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

@dataclass(frozen=True)
class Source:
    """Matches Fortran 'source' structure (heliocentric)."""
    r: float           # AU
    alphaM: float      # rad
    betaM: float       # rad
    rrM: Vec3          # AU, stored 3-vector
    zeta: float        # rad
    eta: float         # rad
    symmetry_axis: Vec3
    ejection_angle_distr: int
    ud: EjectionSpeedProperties
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

@dataclass(frozen=True)
class Comet:
    """Matches Fortran 'comet' structure."""
    coords: Vec3     # AU
    Vastvec: Vec3    # AU/day
    Vast: float      # AU/day
    def __post_init__(self):
        _ = as_vec3(self.coords, "Comet.coords")
        v = as_vec3(self.Vastvec, "Comet.Vastvec")
        if not (math.isfinite(self.Vast) and self.Vast >= 0.0):
            raise ValueError("Comet.Vast must be finite and >= 0.")
        # Optional: don't enforce Vast == ||Vastvec||, but document expectation.
