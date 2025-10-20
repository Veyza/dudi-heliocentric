#!/usr/bin/env python3
"""
select_method.py — Python equivalent of Fortran examples/select_method.f90

Location: python_interface/examples/select_method.py

It reads input parameters from ../../input_data_files/orbit_and_time_test.dat,
computes dust number density on a 200x200 grid using three DUDI-hc methods
(v-integration, delta-ejection, simple expansion), saves matrices and
discrepancy maps into ../../results, and prints method recommendations.

Prereqs:
  1) Build the ctypes bridge (libpy_dudihc_bridge.so)
  2) `pip install -e .` (optional) so that `python_interface` is importable
Run:
  python python_interface/examples/select_method.py
"""
from __future__ import annotations

import math
from pathlib import Path
import numpy as np

# Import the thin Python API and datamodels
from python_interface.dudi_hc.api import v_integration, delta_ejection, simple_expansion
from python_interface.dudi_hc.models import Point, Source, Comet, EjectionSpeedProperties

# ---- physical constants (SI unless noted) ----
AU_M = 1.495978707e11                 # meters
DAY_S = 86400.0                       # seconds
MPS_to_AU_PER_DAY = DAY_S / AU_M      # convert m/s -> AU/day
GMSUN_AU3_PER_DAY2 = (0.01720209895 ** 2)  # Gaussian gravitational constant squared

ACCURACY_PERCENT = 5.0                # Fortran: real, parameter :: accuracy = 5.0
NT1 = 200
NT2 = 200
PERICENTER = False                    # Fortran: pericenter = .FALSE.


import math
import numpy as np

def propagate_two_body(r0, v0, mu, time, *, prefer_scipy=True):
    """
    Propagate a point mass under central gravity using RK4 (Fortran-equivalent).

    Parameters
    ----------
    r0 : array-like (3,)
        Initial heliocentric position [AU].
    v0 : array-like (3,)
        Initial heliocentric velocity [AU/day].
    mu : float
        Gravitational parameter in [AU^3/day^2] (e.g., G*M_sun).
    time : float
        Propagation duration [days], can be positive or negative.
    prefer_scipy : bool
        If True and SciPy is available, use solve_ivp; otherwise use RK4.

    Returns
    -------
    rT, vT : (3,), (3,)
        Propagated position and velocity at t = time, in AU and AU/day.
    """
    r0 = np.asarray(r0, dtype=float)
    v0 = np.asarray(v0, dtype=float)

    # --- optional SciPy path ---
    if prefer_scipy:
        try:
            from scipy.integrate import solve_ivp  # type: ignore
            def rhs(t, y):
                r = y[:3]
                v = y[3:]
                r2 = float(np.dot(r, r))
                if r2 == 0.0:
                    a = np.zeros(3)
                else:
                    a = -mu * r / (r2 * math.sqrt(r2))
                return np.hstack([v, a])

            y0 = np.hstack([r0, v0])
            # Tight tolerances to match fixed-step RK4 accuracy
            # DOP853 is excellent; fall back to RK45 if unavailable
            try:
                method = "DOP853"
            except Exception:
                method = "RK45"

            sol = solve_ivp(rhs, (0.0, float(time)), y0,
                            rtol=1e-10, atol=1e-12, method=method)
            rT = sol.y[:3, -1]
            vT = sol.y[3:, -1]
            return rT, vT
        except Exception:
            # Fall back to manual RK4 below
            pass

    # --- manual RK4 (mirrors your Fortran subroutine) ---
    r = r0.copy()
    v = v0.copy()

    # Match Fortran's adaptive Nstep logic: start at 200, grow by 20% until dt <= 3e-4
    Nstep = 200
    if time == 0.0:
        return r, v

    total_time = float(time)
    # Keep the sign of time; use positive step size for loop and apply direction in dt
    sign = 1.0 if total_time >= 0.0 else -1.0
    T = abs(total_time)

    dt = T / Nstep
    while dt > 3e-4:
        Nstep = int(Nstep * 1.2) + 1
        dt = T / Nstep

    dt *= sign  # bring back the direction

    def acc(rr):
        r2 = float(np.dot(rr, rr))
        if r2 == 0.0:
            return np.zeros(3)
        return -mu * rr / (r2 * math.sqrt(r2))

    for _ in range(Nstep):
        k1 = acc(r)
        l1 = v

        k2 = acc(r + l1 * (dt * 0.5))
        l2 = v + k1 * (dt * 0.5)

        k3 = acc(r + l2 * (dt * 0.5))
        l3 = v + k2 * (dt * 0.5)

        k4 = acc(r + l3 * dt)
        l4 = v + k3 * dt

        v += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        r += (dt / 6.0) * (l1 + 2.0 * l2 + 2.0 * l3 + l4)

    return r, v



def _read_test_input(path: Path):
    """
    orbit_and_time_test.dat with inline comments:
      line 1: 3 floats  -> comet coords (AU)                 # comment...
      line 2: 3 floats  -> comet Vastvec (AU/day)            # comment...
      line 3: 1 float   -> tnow (days)                       # comment...
      line 4: 1 float   -> dtau (seconds)                    # comment...

    Any text after '#' is ignored. Blank lines are skipped.
    """
    import numpy as np

    needed = [3, 3, 1, 1]
    parsed = []

    with path.open("r") as f:
        for need in needed:
            # Advance until we parse enough numbers from a non-empty, non-comment line
            while True:
                line = f.readline()
                if line == "":  # EOF
                    raise ValueError(
                        "Unexpected end of file while reading orbit_and_time_test.dat."
                    )
                # Strip inline comments and surrounding whitespace
                clean = line.split("#", 1)[0].strip()
                if not clean:
                    continue  # skip empty/comment-only lines
                nums = np.fromstring(clean, sep=" ", dtype=float)
                if nums.size < need:
                    raise ValueError(
                        f"Expected at least {need} numeric value(s) on a data line, got {nums.size}."
                    )
                parsed.append(nums[:need].astype(float))
                break

    coords = parsed[0]
    vvec = parsed[1]
    tnow = float(parsed[2][0])
    dtau_s = float(parsed[3][0])
    return coords, vvec, tnow, dtau_s



def _orthonormal_basis(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return two unit vectors spanning the plane orthogonal to `normal`.
    If `normal` is near zero, fall back to the global z-axis as normal."""
    n = np.asarray(normal, dtype=float)
    n_norm = float(np.linalg.norm(n))
    if not math.isfinite(n_norm) or n_norm < 1e-15:
        n = np.array([0.0, 0.0, 1.0], float)
        n_norm = 1.0
    n = n / n_norm
    helper = np.array([1.0, 0.0, 0.0], float) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0], float)
    u = np.cross(n, helper); u /= np.linalg.norm(u)
    v = np.cross(n, u);      v /= np.linalg.norm(v)
    return u, v


def _cart_to_spherical(rvec: np.ndarray) -> tuple[float, float, float]:
    """Return (r, alpha, beta) for vector rvec (AU).
    alpha: polar angle (0..pi), beta: eastern longitude (-pi..pi)."""
    x, y, z = map(float, rvec)
    r = math.sqrt(x*x + y*y + z*z)
    if r == 0.0:
        alpha = 0.0; beta = 0.0
    else:
        alpha = math.acos(z / r)
        beta  = math.atan2(y, x)
    return r, alpha, beta


def orbital_plane_grid(nt1: int, nt2: int, resolution_m: tuple[float, float],
                       comet, center: np.ndarray):
    """
    Fortran-equivalent grid with the '0.95 on x-component' tweak:
      zvec = cross(V, R); zvec = zvec / |zvec|
      xvec = R; xvec(1) = xvec(1) * 0.95; tmpvec = xvec; xvec = xvec / |tmpvec|
      yvec = cross(zvec, xvec)
    """
    AU_M = 1.495978707e11

    R = np.asarray(comet.coords, dtype=float)
    V = np.asarray(comet.Vastvec, dtype=float)

    # z-axis normal to orbital plane: z = V × R, then normalize
    zvec = np.cross(V, R)
    nz = np.linalg.norm(zvec)
    if not np.isfinite(nz) or nz < 1e-15:
        zvec = np.array([0.0, 0.0, 1.0], float)
    else:
        zvec = zvec / nz

    # x-axis ~ along heliocentric radius but with a small x-tweak to avoid bad geometry
    xvec = R.copy()
    xvec[0] *= 0.95  # same tweak as your Fortran: xvec(1) = xvec(1) * 0.95
    tmpvec = xvec.copy()
    nx = np.linalg.norm(tmpvec)
    if not np.isfinite(nx) or nx < 1e-15:
        xvec = np.array([1.0, 0.0, 0.0], float)
    else:
        xvec = xvec / nx

    # y-axis completes the right-handed triad
    yvec = np.cross(zvec, xvec)

    # Convert resolution from meters to AU for step vectors
    xstep = xvec * (resolution_m[0] / AU_M)
    ystep = yvec * (resolution_m[1] / AU_M)

    # Lower-left corner (Fortran: center - nt1*xstep/2 - nt2*ystep/2)
    origin = center - nt1 * xstep * 0.5 - nt2 * ystep * 0.5

    pts = np.empty((nt1, nt2), dtype=object)
    for j in range(nt2):          # ii = 1..nt2
        for i in range(nt1):      #  i = 1..nt1
            rvec = origin + (i + 1) * xstep + (j + 1) * ystep
            r, alpha, beta = _cart_to_spherical(rvec)
            pts[i, j] = Point(r=r, alpha=alpha, beta=beta, rvector=rvec.astype(float))
    return pts

from concurrent.futures import ProcessPoolExecutor, as_completed
import os

def _eval_row(j, NT1, origin, xstep, ystep,
              source, comet, muR, tnow, Rast_AU, PERICENTER, cloud_center):
    # Compute one row (fixed j) of densities
    row_d = np.empty(NT1, dtype=float)
    row_v = np.empty(NT1, dtype=float)
    row_s = np.empty(NT1, dtype=float)

    for i in range(NT1):
        rvec = origin + (i + 1) * xstep + (j + 1) * ystep
        r, alpha, beta = _cart_to_spherical(rvec)
        pt = Point(r=r, alpha=alpha, beta=beta, rvector=rvec.astype(float))

        row_d[i] = delta_ejection(pt, source, comet, muR=muR, dt=tnow, Rast_AU=Rast_AU)
        row_v[i] = v_integration(pt, source, comet, muR=muR, tnow=tnow, Rast_AU=Rast_AU, pericenter=PERICENTER)
        row_s[i] = simple_expansion(pt, source, cloudcentr=cloud_center, dt=tnow)

    return j, row_d, row_v, row_s




def main() -> int:
    # Resolve repository root from this file location
    here = Path(__file__).resolve()
    repo_root = here.parents[2]  # .../ (root)
    input_path = repo_root / "input_data_files" / "orbit_and_time_test.dat"
    results_dir = repo_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    # Read test parameters
    coords, vastvec, tnow, dtau_s = _read_test_input(input_path)
    vast = float(np.linalg.norm(vastvec))

    # Source and comet (mirrors Fortran block)
    r = float(np.linalg.norm(coords))
    alphaM = math.acos(float(coords[2]) / r)
    betaM = math.atan2(float(coords[1]), float(coords[0]))
    axis = coords / r

    ud = EjectionSpeedProperties(
        ud_shape=0,                         # uniform speed PDF
        umin=0.5 * MPS_to_AU_PER_DAY,       # 0.5 m/s
        umax=50.0 * MPS_to_AU_PER_DAY,      # 50 m/s
    )
    source = Source(
        r=r, alphaM=alphaM, betaM=betaM,
        rrM=coords,
        zeta=0.0, eta=0.0,
        symmetry_axis=axis,
        ejection_angle_distr=1,
        ud=ud,
        Nparticles=1.0e10,
        Tj=0.0,
        dtau=2.0 * (dtau_s / DAY_S),        # Fortran: 2*dtau / s_in_day
    )
    comet = Comet(coords=coords, Vastvec=vastvec, Vast=vast)
    print(comet)

    muR = GMSUN_AU3_PER_DAY2
    Rast_AU = 0.0

    # Plane through the middle of the cloud.
    # Fortran uses runge_kutta_point_position; here we approximate:
    cloud_center, _ = propagate_two_body(coords, vastvec, mu=muR, time=tnow, prefer_scipy=False)

    # --- build orbital-plane grid around the cloud center and compute densities ---
    # Resolution in *meters* (matches Fortran): umax [AU/day] * tnow [day] -> [AU], then * AU_M -> [m]
    resolution_m = (
        float(source.ud.umax * tnow * AU_M) / NT1,
        float(source.ud.umax * tnow * AU_M) / NT2,
    )
    points = orbital_plane_grid(NT1, NT2, resolution_m, comet, cloud_center)

        # --- Parallel evaluation over rows with ProcessPool ---
    # Recreate the same orbital-plane frame as in orbital_plane_grid:
    #   z = V x R; normalize
    zvec = np.cross(comet.Vastvec, comet.coords)
    nz = np.linalg.norm(zvec)
    zvec = zvec / nz if (np.isfinite(nz) and nz > 1e-15) else np.array([0.0, 0.0, 1.0])

    #   x = R; x[0]*=0.95; normalize (same degeneracy-avoidance as Fortran)
    xvec = comet.coords.copy()
    xvec[0] *= 0.95
    nx = np.linalg.norm(xvec)
    xvec = xvec / nx if (np.isfinite(nx) and nx > 1e-15) else np.array([1.0, 0.0, 0.0])

    #   y = z × x
    yvec = np.cross(zvec, xvec)

    # Step vectors from your existing resolution (meters) -> AU
    xstep = xvec * (resolution_m[0] / AU_M)
    ystep = yvec * (resolution_m[1] / AU_M)

    # Lower-left corner
    origin = cloud_center - NT1 * xstep * 0.5 - NT2 * ystep * 0.5

    dens_s = np.zeros((NT1, NT2), dtype=float)
    dens_d = np.zeros_like(dens_s)
    dens_v = np.zeros_like(dens_s)

    max_workers = min(os.cpu_count() or 1, NT2)
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = [
            ex.submit(
                _eval_row, j, NT1, origin, xstep, ystep,
                source, comet, muR, tnow, Rast_AU, PERICENTER, cloud_center
            )
            for j in range(NT2)
        ]
        for fut in as_completed(futures):
            j, row_d, row_v, row_s = fut.result()
            dens_d[:, j] = row_d
            dens_v[:, j] = row_v
            dens_s[:, j] = row_s


    # Exclude the center (consistent with Fortran)
    dx_AU = resolution_m[0] / AU_M
    cx, cy = NT1 // 2, NT2 // 2
    k = int(source.ud.umin * tnow / dx_AU) + 1
    x0 = max(0, cx - k); x1 = min(NT1, cx + k + 1)
    y0 = max(0, cy - k); y1 = min(NT2, cy + k + 1)
    dens_d[x0:x1, y0:y1] = 1.0
    dens_v[x0:x1, y0:y1] = 1.0
    dens_s[x0:x1, y0:y1] = 1.0

    # Save matrices
    def _save_matrix(path: Path, A: np.ndarray) -> None:
        np.savetxt(path, A, fmt="%.9e")

    _save_matrix(results_dir / "py_test_simple_exp_meth.dat", dens_s)
    _save_matrix(results_dir / "py_test_delta-eject_meth.dat", dens_d)
    _save_matrix(results_dir / "py_test_v-integr_meth.dat", dens_v)

    # Discrepancies
        # Discrepancies & stats (now also print detailed min/max in % as in Fortran)
    with np.errstate(divide="ignore", invalid="ignore"):
        test_d = dens_d / dens_v - 1.0   # delta-ejection vs v-integration
        test_s = dens_s / dens_d - 1.0   # simple expansion vs delta-ejection
        test_d[~np.isfinite(test_d)] = 0.0
        test_s[~np.isfinite(test_s)] = 0.0

    _save_matrix(results_dir / "py_test_delta-eject_vs_v-integr.dat", test_d)
    _save_matrix(results_dir / "py_test_simp_exp_vs_delta-eject.dat", test_s)

    # Print detailed statistics (percent)
    min_d = float(np.min(test_d)) * 100.0
    max_d = float(np.max(test_d)) * 100.0
    min_s = float(np.min(test_s)) * 100.0
    max_s = float(np.max(test_s)) * 100.0

    print(" difference between the delta-ejection solution and v-integration solution")
    print(f"min {min_d:8.1f}%")
    print(f"max {max_d:8.1f}%")
    print(" difference between the simple expansion solution and delta-ejection solution")
    print(f"min {min_s:8.1f}%")
    print(f"max {max_s:8.1f}%")

    # Recommendations (keep 5% threshold; wording mirrors Fortran)
    mean_abs_d = float(np.mean(np.abs(test_d)))
    mean_abs_s = float(np.mean(np.abs(test_s)))
    extreme_d = max(abs(min_d), abs(max_d))  # already in percent
    extreme_s = max(abs(min_s), abs(max_s))  # already in percent

    if (mean_abs_d < ACCURACY_PERCENT * 1e-2) and (extreme_d < ACCURACY_PERCENT):
        print(" delta-ejection method is applicable")
    else:
        print(" v-integration method is recommended")

    if (mean_abs_s < ACCURACY_PERCENT * 1e-2) and (extreme_s < ACCURACY_PERCENT):
        print(" simple expansion method is applicable too")
    else:
        print(" simple expansion method is NOT recommended")



if __name__ == "__main__":
    raise SystemExit(main())
