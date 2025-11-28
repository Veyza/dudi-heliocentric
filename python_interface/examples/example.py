#!/usr/bin/env python3
"""
python_interface/examples/example.py
Python analogue of Fortran examples/example.f90

- Builds sources on a sphere around the comet at multiple ephemeris points
- Evaluates dust number density on an orbital-plane grid
- Uses the delta-ejection solution (as in example.f90)
- Writes results/result.dat (matrix)

NOTE: This is a faithful port; without OpenMP it will be slow in pure Python.
"""

from __future__ import annotations

from pathlib import Path
import math
import numpy as np

from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import time

# ---- import your thin API and dataclasses ----
from python_interface.dudi_hc.api import delta_ejection, batch_over_points, batch_over_sources, batch_over_points_sources
from python_interface.dudi_hc.models import Point, Source, Comet, EjectionSpeedProperties

# ---- constants (match const.f90) ----
PI = math.pi
AU_M = 1.495978707e11          # m
DAY_S = 86400.0                # s
MPS_to_AU_PER_DAY = DAY_S / AU_M
G = 6.674e-11                  # m^3 / kg / s^2
c = 2.99792458e8               # m/s
Msun = 1.9891e30               # kg
Lsun = 3.828e26                # W
rho = 2.5e3                    # kg/m^3

# ---- reuse the same orbital-plane grid you implemented for select_method.py ----
def orbital_plane_grid(nt1: int, nt2: int, resolution_m: tuple[float, float],
                       comet, center: np.ndarray):
    AU_M = 1.495978707e11

    R = np.asarray(comet.coords, dtype=float)
    V = np.asarray(comet.Vastvec, dtype=float)

    zvec = np.cross(V, R)
    nz = np.linalg.norm(zvec)
    zvec = zvec / nz if (np.isfinite(nz) and nz > 1e-15) else np.array([0.0, 0.0, 1.0])

    xvec = R.copy()
    xvec[0] *= 0.95  # avoid bad geometry (exactly like your Fortran)
    nx = np.linalg.norm(xvec)
    xvec = xvec / nx if (np.isfinite(nx) and nx > 1e-15) else np.array([1.0, 0.0, 0.0])

    yvec = np.cross(zvec, xvec)

    xstep = xvec * (resolution_m[0] / AU_M)
    ystep = yvec * (resolution_m[1] / AU_M)
    origin = center - nt1 * xstep * 0.5 - nt2 * ystep * 0.5

    def cart_to_spherical(rv):
        x, y, z = map(float, rv)
        r = math.sqrt(x*x + y*y + z*z)
        if r == 0.0:
            return 0.0, 0.0, 0.0
        alpha = math.acos(z / r)
        beta = math.atan2(y, x)
        return r, alpha, beta

    pts = np.empty((nt1, nt2), dtype=object)
    for j in range(nt2):
        for i in range(nt1):
            rvec = origin + (i + 1) * xstep + (j + 1) * ystep
            r, alpha, beta = cart_to_spherical(rvec)
            pts[i, j] = Point(r=r, alpha=alpha, beta=beta, rvector=rvec.astype(float))
    return pts

# ---- helpers that mirror your Fortran subroutines ----
def reduced_gravitational_parameter(Rg_m: float, Qpr: float) -> float:
    """
    Fortran:
      radiation = 3/(16π) * Lsun * Qpr / (Rg * rho * c)
      muR = G*Msun - radiation
      muR in AU^3/day^2
    """
    radiation = (3.0 / (16.0 * PI)) * Lsun * Qpr / (Rg_m * rho * c)   # m^3/s^2
    mu_SI = G * Msun - radiation                                      # m^3/s^2
    mu_AU_day2 = mu_SI / (AU_M**3) * (DAY_S**2)                       # AU^3/day^2
    return mu_AU_day2

def _uniform_points_on_unit_sphere(N: int) -> np.ndarray:
    # Marsaglia (1972)
    u = np.random.uniform(-1.0, 1.0, size=N)
    theta = np.random.uniform(0.0, 2.0 * PI, size=N)
    s = np.sqrt(1.0 - u**2)
    return np.stack([s * np.cos(theta), s * np.sin(theta), u], axis=1)

def _clamp(x, lo=-1.0, hi=1.0):
    return max(lo, min(hi, x))

def _eta_azimuth(sym_axis: np.ndarray, rrM_vec: np.ndarray) -> float:
    """
    A consistent azimuth around r̂ = rrM/|rrM|.
    """
    rhat = rrM_vec / np.linalg.norm(rrM_vec)
    # choose a stable "north" to build a local basis
    z = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(z, rhat)) > 0.99:
        z = np.array([1.0, 0.0, 0.0])
    e_theta = z - np.dot(z, rhat) * rhat
    e_theta /= np.linalg.norm(e_theta)
    e_phi = np.cross(rhat, e_theta)
    comp_theta = float(np.dot(sym_axis, e_theta))
    comp_phi = float(np.dot(sym_axis, e_phi))
    return math.atan2(comp_phi, comp_theta)

def get_sources(ephem_path: Path, Np: int, Ns: int, Rast_AU: float):
    """
    Python analogue of Fortran get_sources.
    Reads ephemeridae.dat: 7 columns per row -> t [day], coords[3] [AU], Vastvec[3] [AU/day]
    Returns:
      sources: list of list [Np][Ns] of Source
      comet:   list [Np] of Comet
    """
    data = np.loadtxt(ephem_path)
    if data.ndim != 2 or data.shape[1] < 7:
        raise ValueError("Ephemeridae must have at least 7 columns (t, x, y, z, vx, vy, vz).")
    if data.shape[0] < Np:
        raise ValueError(f"Ephemeridae has {data.shape[0]} rows; expected at least {Np}.")
    data = data[:Np, :7]

    comet = []
    for i in range(Np):
        t, x, y, z, vx, vy, vz = data[i]
        coords = np.array([x, y, z], float)
        Vastvec = np.array([vx, vy, vz], float)
        Vast = float(np.linalg.norm(Vastvec))
        comet.append(Comet(coords=coords, Vastvec=Vastvec, Vast=Vast))

    # build sources
    sources = [[None for _ in range(Ns)] for _ in range(Np)]
    for i in range(Np):
        # uniform points on sphere (unit vectors)
        xyz = _uniform_points_on_unit_sphere(Ns)
        for ii in range(Ns):
            axis = xyz[Ns - 1 - ii]  # mirror the Fortran Ns+1-ii order
            rrM = comet[i].coords + axis * Rast_AU
            r = float(np.linalg.norm(rrM))
            alphaM = 0.0 if r == 0.0 else math.acos(float(rrM[2]) / r)
            betaM = 0.0 if r == 0.0 else math.atan2(float(rrM[1]), float(rrM[0]))
            tmp = float(np.dot(axis, rrM / r)) if r > 0.0 else 1.0
            tmp = _clamp(tmp, -1.0, 1.0)
            zeta = math.acos(tmp)
            eta = _eta_azimuth(axis, rrM) if zeta > 1e-8 else 0.0

            ud = EjectionSpeedProperties(
                ud_shape=0,
                umin=5.0 * MPS_to_AU_PER_DAY,      # 5 m/s
                umax=100.0 * MPS_to_AU_PER_DAY,    # 100 m/s
            )
            # per Fortran:
            Tj = float(data[i, 0])                 # moment [day]
            dtau = 1e-2 / DAY_S                    # days
            Nparticles = 1e5
            if comet[i].Vast > 0.0:
                Nparticles += 1e5 * float(np.dot(comet[i].Vastvec, axis)) / comet[i].Vast

            sources[i][ii] = Source(
                r=r, alphaM=alphaM, betaM=betaM,
                rrM=rrM, zeta=zeta, eta=eta,
                symmetry_axis=axis, ejection_angle_distr=2,
                ud=ud, Nparticles=Nparticles, Tj=Tj, dtau=dtau,
            )
    return sources, comet

def matrix_out(path: Path, image: np.ndarray):
    path.parent.mkdir(parents=True, exist_ok=True)
    # ES12.4E2-like scientific format; numpy doesn't support 2-digit exponent width, but this is close.
    np.savetxt(path, image, fmt="%.4e")


def _cart_to_spherical_vec(rvec):
    x, y, z = float(rvec[0]), float(rvec[1]), float(rvec[2])
    r = math.sqrt(x*x + y*y + z*z)
    if r == 0.0:
        return 0.0, 0.0, 0.0
    alpha = math.acos(z / r)
    beta = math.atan2(y, x)
    return r, alpha, beta


def main() -> int:
    # repo root = ../../ from this file
    here = Path(__file__).resolve()
    repo_root = here.parents[2]
    ephem_path = repo_root / "input_data_files" / "ephemeridae.dat"
    out_path = repo_root / "results" / "result.dat"

    # parameters from the Fortran example
    Nt = 41          # number of time steps (was Np)
    Ns = 50          # number of sources per time
    Rast_m = 5e3
    Rast_AU = Rast_m / AU_M
    Qpr = 0.5
    Rg_m = 0.29e-6
    nt1, nt2 = 200, 200
    resolution_m = (2e3, 2e3)

    # μ_R (AU^3/day^2)
    muR = reduced_gravitational_parameter(Rg_m, Qpr)

    # ephemeris & sources
    sources, comet = get_sources(ephem_path, Nt, Ns, Rast_AU)

    # time now = moment at last ephemeris row
    tnow = sources[-1][0].Tj

    # Reconstruct the same plane frame used by orbital_plane_grid
    R = np.asarray(comet[-1].coords, dtype=float)
    V = np.asarray(comet[-1].Vastvec, dtype=float)

    # z = V × R (normal to orbital plane)
    zvec = np.cross(V, R)
    nz = np.linalg.norm(zvec)
    zvec = zvec / nz if (np.isfinite(nz) and nz > 1e-15) else np.array([0.0, 0.0, 1.0])

    # x = R, tweak x[0]*=0.95 to avoid degeneracy (same as Fortran)
    xvec = R.copy()
    xvec[0] *= 0.95
    nx = np.linalg.norm(xvec)
    xvec = xvec / nx if (np.isfinite(nx) and nx > 1e-15) else np.array([1.0, 0.0, 0.0])

    # y = z × x
    yvec = np.cross(zvec, xvec)

    # Step vectors (meters -> AU)
    xstep = xvec * (resolution_m[0] / AU_M)
    ystep = yvec * (resolution_m[1] / AU_M)

    # Lower-left corner of the grid (same as in orbital_plane_grid)
    origin = comet[-1].coords - nt1 * xstep * 0.5 - nt2 * ystep * 0.5

    # ------------------------------------------------------------------
    # Build the full grid of Points once
    # ------------------------------------------------------------------
    points: list[Point] = []
    for j in range(nt2):
        for i in range(nt1):
            rvec = origin + (i + 1) * xstep + (j + 1) * ystep
            r, alpha, beta = _cart_to_spherical_vec(rvec)
            pt = Point(r=r, alpha=alpha, beta=beta, rvector=rvec.astype(float))
            points.append(pt)

    # ------------------------------------------------------------------
    # Use the new time×sources×points batching (delta-ejection)
    # We mimic the old behaviour: i_t = 0..Nt-2 (exclude last row where dt=0)
    # ------------------------------------------------------------------
    sources_by_time = sources[: Nt - 1]   # shape (Nt-1, Ns)
    comets_by_time = comet[: Nt - 1]      # length Nt-1

    dens_flat = batch_over_points_sources(
        points=points,
        sources_by_time=sources_by_time,
        comets_by_time=comets_by_time,
        muR=muR,
        tnow=tnow,
        Rast_AU=Rast_AU,
        pericenter=False,          # not used by delta_ejection
        method="delta_ejection",
    )

    # reshape into (nt1, nt2)
    density = dens_flat.reshape(nt1, nt2)

    matrix_out(out_path, density)
    print(f"Wrote {out_path}")
    return 0



if __name__ == "__main__":
    raise SystemExit(main())
