from __future__ import annotations

from typing import Iterable, Sequence
import numpy as np
import types

from python_interface.dudi_hc.models import Point, Source, Comet, EjectionSpeedProperties
from python_interface.dudi_hc import api


# ----------------------------------------------------------------------
# Local numeric helpers and constants
# (If you have a Python constants module mirroring const.f90,
# you can import from there instead of using these definitions.)
# ----------------------------------------------------------------------

PI = np.pi
HALFPI = 0.5 * PI
TWOPI = 2.0 * PI
DEG2RAD = PI / 180.0

# Physical units
AU = 1.495978707e11       # [m]
S_IN_DAY = 86400.0        # [s]
# Fortran AUdays2SI ≈ AU / seconds_per_day
AUdays2SI = AU / S_IN_DAY


def _linterpol(x: np.ndarray, y: np.ndarray, x_new: float) -> float:
    """
    1-D linear interpolation, Python equivalent of Fortran LiNTERPOL.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape")
    return float(np.interp(x_new, x, y))


# ----------------------------------------------------------------------
# get_maps_data
# ----------------------------------------------------------------------

def get_maps_data() -> tuple[np.ndarray, list[str]]:
    """
    Python translation of Fortran get_maps_data.

    Returns
    -------
    rhels : ndarray, shape (N,)
        Heliocentric distances [AU].
    fnames : list of str, length N
        Paths to the impact-ejecta map files.
    """
    N = 4
    rhels = np.array([0.96, 0.94, 0.93, 0.91], dtype=np.float64)
    fnames: list[str] = []

    for i in range(1, N + 1):
        # Fortran: write(ic,'(I1)') i+4  -> '5','6','7','8'
        ic = f"{i + 4:d}"
        fname = (
            "./input_data_files/impact_ejecta_maps/"
            f"COS3_Yield_ALL_map_{ic}.sav_yield_orb_pha.txt"
        )
        fnames.append(fname)

    return rhels, fnames


# ----------------------------------------------------------------------
# get_moving_sources
# ----------------------------------------------------------------------

def get_moving_sources(
    fname: str,
    Np: int,
    Nlin: int,
) -> tuple[list[Source], list[Comet]]:
    """
    Python translation of Fortran get_moving_sources.

    Parameters
    ----------
    fname : str
        Ephemeris file; each line contains:
        t, x, y, z, vx, vy, vz  (free format).
    Np : int
        Number of time steps / sources.
    Nlin : int
        Step between ephemeris lines. Intermediate entries are filled
        by linear interpolation.

    Returns
    -------
    sources : list[Source]
        One source per time step (length Np).
    comet : list[Comet]
        Comet state at each time step (length Np).
    """
    # ----- map file names and distances -----
    rhels, fnames = get_maps_data()
    if len(fnames) < 2:
        raise RuntimeError("Need at least two ratemap files.")

    # Fortran: mapind1=1, mapind2=mapind1+1 (1-based)
    mapind1 = 0
    mapind2 = 1

    rhel1 = api.read_first_ratemap(fnames[mapind1])
    rhel2 = api.read_ratemap(fnames[mapind2])

    def _copy_rmap2_to_rmap1() -> None:
        """Fortran 'rmap1 = rmap2' via api."""
        rmap2 = api.get_rmap2()
        api.set_rmap1(rmap2)

    # ----- read ephemeris and interpolate (Fortran logic) -----
    moment = np.zeros(Np, dtype=np.float64)
    coords = np.zeros((Np, 3), dtype=np.float64)
    Vastvec = np.zeros((Np, 3), dtype=np.float64)

    dNlin = float(Nlin)

    with open(fname, "r") as f:
        # First line: Fortran read(200,*) moment(1), coords, Vastvec
        first = f.readline()
        if not first:
            raise RuntimeError(f"Ephemeris file '{fname}' is empty.")
        vals = [float(tok) for tok in first.split()]
        if len(vals) < 7:
            raise RuntimeError(
                "Each ephemeris line must contain at least 7 numbers: "
                "t, x, y, z, vx, vy, vz"
            )

        moment[0] = vals[0]
        coords[0, :] = vals[1:4]
        Vastvec[0, :] = vals[4:7]

        # Fortran:
        #   do i = Nlin+1, Np, Nlin       (i is 1-based)
        #   read(200,*) moment(i), comet(i)%coords, comet(i)%Vastvec
        #   forall(ii = (i-Nlin+1):(i-1)) ...
        #
        # Python: i_fortran -> i_python = i_fortran - 1
        for i_fortran in range(Nlin + 1, Np + 1, Nlin):
            line = f.readline()
            if not line:
                raise RuntimeError(
                    "Ephemeris file ended prematurely while filling comet array."
                )
            vals = [float(tok) for tok in line.split()]
            if len(vals) < 7:
                raise RuntimeError(
                    "Each ephemeris line must contain at least 7 numbers: "
                    "t, x, y, z, vx, vy, vz"
                )

            i = i_fortran - 1  # 0-based
            moment[i] = vals[0]
            coords[i, :] = vals[1:4]
            Vastvec[i, :] = vals[4:7]

            if Nlin > 1:
                # Fortran: ii = (i-Nlin+1):(i-1) (1-based)
                for ii_fortran in range(i_fortran - Nlin + 1, i_fortran):
                    j = ii_fortran - 1  # 0-based index
                    # weight w = (ii - i + Nlin) / dNlin
                    w = float(ii_fortran - i_fortran + Nlin) / dNlin

                    moment[j] = (
                        moment[i - Nlin]
                        + (moment[i] - moment[i - Nlin]) * w
                    )
                    Vastvec[j, :] = (
                        Vastvec[i - Nlin, :]
                        + (Vastvec[i, :] - Vastvec[i - Nlin, :]) * w
                    )
                    coords[j, :] = (
                        coords[i - Nlin, :]
                        + (coords[i, :] - coords[i - Nlin, :]) * w
                    )

    # Shift timeline: moment(i) = moment(i) - moment(1)
    moment -= moment[0]

    # Phaethon radius [m], as in Fortran
    Rast = 2.9e3

    # Precompute dt between steps in days (used for Nparticles)
    dt_days = float(moment[1] - moment[0]) if Np > 1 else 0.0

    # ----- build Comet objects -----
    comets: list[Comet] = []
    for i in range(Np):
        Vast_i = float(np.linalg.norm(Vastvec[i, :]))
        comets.append(
            Comet(
                coords=coords[i, :].copy(),
                Vastvec=Vastvec[i, :].copy(),
                Vast=Vast_i,
            )
        )
        #print(coords[i, :], Vastvec[i, :])

    # ----- build Source objects (1 per time step) -----
    sources: list[Source] = []

    for i in range(Np):
        rrM = coords[i, :].astype(np.float64)
        r = float(np.linalg.norm(rrM))

        # update which maps we interpolate between when r drops below rhel2
        if r < rhel2 and mapind2 + 1 < len(fnames):
            rhel1 = rhel2
            mapind1 = mapind2
            mapind2 += 1
            _copy_rmap2_to_rmap1()
            rhel2 = api.read_ratemap(fnames[mapind2])

        # Interpolate ratemap for this r (Fortran ratematr_interpolate)
        api.ratematr_interpolate(rhel=r, rhel1=rhel1, rhel2=rhel2)

        # Angular coordinates
        if r > 0.0:
            alphaM = float(np.arccos(rrM[2] / r))
            betaM = float(np.arctan2(rrM[1], rrM[0]))
        else:
            alphaM = 0.0
            betaM = 0.0

        # Symmetry axis: Fortran broadcasts scalar rrM(1)/r into the array
        if r > 0.0:
            axis_val = rrM / r
        else:
            axis_val = 0.0
        symmetry_axis = np.array(
            axis_val,
            dtype=np.float64,
        )

        # Distributions / timing parameters (from Fortran)
        zeta = 0.0
        eta = 0.0
        ud = EjectionSpeedProperties(
            ud_shape=1,
            umin=2.0 / AUdays2SI,
            umax=2399.0 / AUdays2SI,
        )
        ejection_angle_distr = 3
        Tj = float(moment[i])
        dtau = 0.0

        # Integrate number density of impact ejecta over the map
        totrate = integrate_over_matrix()

        # Convert number density to flux (Szalay et al. 2016, Eq. 3)
        totrate = totrate / 0.31 / 7.2e-3 / 4.0 / PI * Rast**2

        # Convert flux to number of ejected particles in this time step
        Nparticles = float(totrate * dt_days * S_IN_DAY)

        sources.append(
            Source(
                r=r,
                alphaM=alphaM,
                betaM=betaM,
                rrM=rrM.copy(),
                zeta=zeta,
                eta=eta,
                symmetry_axis=symmetry_axis,
                ejection_angle_distr=ejection_angle_distr,
                ud=ud,
                Nparticles=Nparticles,
                Tj=Tj,
                dtau=dtau,
            )
        )

    return sources, comets




# ----------------------------------------------------------------------
# get_flyby_trajectory
# ----------------------------------------------------------------------

def get_flyby_trajectory(
    n1: int,
    resolution: float,
    CAdist: float,
    lastrM: np.ndarray,
) -> list[Point]:
    """
    Python translation of Fortran get_flyby_trajectory.

    Parameters
    ----------
    n1 : int
        Number of trajectory points.
    resolution : float
        Step along the trajectory [m].
    CAdist : float
        Closest-approach distance along the normal direction [AU].
    lastrM : array_like, shape (3,)
        Asteroid position vector in AU.

    Returns
    -------
    points : list[Point]
    """
    lastrM = np.asarray(lastrM, dtype=np.float64)
    if lastrM.shape != (3,):
        raise ValueError("lastrM must be a 3-vector in AU")

    # CS to compare with Szalay et al. 2019
    zvec = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    xvec = lastrM / np.linalg.norm(lastrM)
    zvec = zvec - zvec * float(np.dot(xvec, zvec))
    zvec = zvec / np.linalg.norm(zvec)
    yvec = np.cross(zvec, xvec)

    angle2xvec = 29.0 * DEG2RAD

    tmpvec = xvec * np.cos(angle2xvec) + yvec * np.sin(angle2xvec)
    tmpvec = tmpvec * (resolution / AU)      # m -> AU
    tmpnorm = -xvec * np.sin(angle2xvec) + yvec * np.cos(angle2xvec)

    CApoint = lastrM - tmpnorm * CAdist     # CAdist already in AU

    points: list[Point] = []
    for i in range(1, n1 + 1):
        rvector = CApoint + tmpvec * (i - n1 / 2.0)
        r = np.linalg.norm(rvector)
        alpha = float(np.arccos(rvector[2] / r))
        beta = float(np.arctan2(rvector[1], rvector[0]))
        p = Point(
            r=r,
            alpha=alpha,
            beta=beta,
            rvector=np.asarray(rvector, dtype=np.float64),
        )
        points.append(p)

    return points


# ----------------------------------------------------------------------
# integrate_over_matrix
# ----------------------------------------------------------------------

def integrate_over_matrix() -> float:
    """
    Vectorized Python translation of Fortran integrate_over_matrix.

    Integrates the current Fortran ratemap (already set via api.read_*
    and api.ratematr_interpolate) over the sphere.
    """
    ratemap = api.get_ratemap()   # shape (nlats, nlons)
    lats = api.get_lats()         # shape (nlats,)
    lons = api.get_lons()         # shape (nlons,)

    nlats, nlons = ratemap.shape
    if lats.shape[0] != nlats or lons.shape[0] != nlons:
        raise RuntimeError("Inconsistent ratemap / lats / lons dimensions")

    # ------------------------------------------------------------------
    # 1) Integrate over longitude for each latitude -> rint[lat]
    # ------------------------------------------------------------------
    # lons differences between adjacent longitudes
    dlons = np.diff(lons)                    # shape (nlons-1,)
    wrap = lons[0] - lons[-1] + TWOPI        # closing segment

    # (ratemap[:,1:] + ratemap[:,:-1]) has shape (nlats, nlons-1)
    # Broadcast dlons across lat dimension and sum over lon
    inner = (ratemap[:, 1:] + ratemap[:, :-1]) * dlons[np.newaxis, :]
    rint = inner.sum(axis=1) + (ratemap[:, 0] + ratemap[:, -1]) * wrap
    rint *= 0.5  # trapezoid in longitude

    # ------------------------------------------------------------------
    # 2) Integrate over latitude (bands between ii-1 and ii)
    # ------------------------------------------------------------------
    polangle = HALFPI - lats  # colatitude

    # main bands: Fortran ii = 2..nlats
    # vectorized:
    #   (rint[ii] + rint[ii-1]) * sin(polangle[ii]) * (polangle[ii-1] - polangle[ii]) / 2
    main = (
        (rint[1:] + rint[:-1])
        * np.sin(polangle[1:])
        * (polangle[:-1] - polangle[1:])
        * 0.5
    )
    integral = main.sum()

    # ------------------------------------------------------------------
    # 3) Top and bottom rings (same formula as Fortran, reuse rint[0] and rint[-1])
    # ------------------------------------------------------------------
    integral += (
        rint[0]
        * np.sin((-HALFPI - lats[0]) / 2.0)
        * (-HALFPI - lats[0])
    )

    integral += (
        rint[-1]
        * np.sin((HALFPI - lats[-1]) / 2.0)
        * (HALFPI - lats[-1])
    )

    return float(integral)



# ----------------------------------------------------------------------
# get_points (2-D grid)
# ----------------------------------------------------------------------

def get_points(
    n1: int,
    n2: int,
    resolution: tuple[float, float],
    lastrM: np.ndarray,
    cntrpx: float,
    cntrpy: float,
) -> np.ndarray:
    """
    Python translation of Fortran get_points.

    Returns a (n1, n2) array of Point objects.
    """
    resx, resy = resolution
    lastrM = np.asarray(lastrM, dtype=np.float64)
    if lastrM.shape != (3,):
        raise ValueError("lastrM must be a 3-vector in AU")

    # CS as in Szalay et al. 2019 comparison
    zvec = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    xvec = lastrM / np.linalg.norm(lastrM)
    zvec = zvec - zvec * float(np.dot(xvec, zvec))
    zvec = zvec / np.linalg.norm(zvec)
    yvec = np.cross(zvec, xvec)

    xvec = xvec * (resx / AU)
    yvec = yvec * (resy / AU)
    tmpvec = lastrM - n1 * xvec * cntrpx - n2 * yvec * cntrpy

    points = np.empty((n1, n2), dtype=object)

    for j in range(n2):      # Fortran ii = 1..n2
        for i in range(n1):  # Fortran i  = 1..n1
            rvector = tmpvec + (i + 1) * xvec + (j + 1) * yvec
            r = np.linalg.norm(rvector)
            alpha = float(np.arccos(rvector[2] / r))
            beta = float(np.arctan2(rvector[1], rvector[0]))
            points[i, j] = Point(
                r=r,
                alpha=alpha,
                beta=beta,
                rvector=np.asarray(rvector, dtype=np.float64),
            )

    return points


# ----------------------------------------------------------------------
# get_points_3d
# ----------------------------------------------------------------------

def get_points_3d(
    nx: int,
    ny: int,
    nz: int,
    resolution: tuple[float, float, float],
    lastrM: np.ndarray,
    cntrpx: float,
    cntrpy: float,
    cntrpz: float,
) -> np.ndarray:
    """
    Python translation of Fortran get_points_3d.

    Returns an (nx, ny, nz) array of Point objects.
    """
    resx, resy, resz = resolution
    lastrM = np.asarray(lastrM, dtype=np.float64)
    if lastrM.shape != (3,):
        raise ValueError("lastrM must be a 3-vector in AU")

    # Start with ecliptic Z
    zvec = np.array([0.0, 0.0, 1.0], dtype=np.float64)

    # x̂: projection of lastrM onto ecliptic plane
    xvec = lastrM.copy()
    xvec = xvec - zvec * float(np.dot(xvec, zvec))
    if np.linalg.norm(xvec) == 0.0:
        xvec = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    else:
        xvec = xvec / np.linalg.norm(xvec)

    # Make ẑ ⟂ x̂, then ŷ = ẑ × x̂
    zvec = zvec - xvec * float(np.dot(xvec, zvec))
    zvec = zvec / np.linalg.norm(zvec)
    yvec = np.cross(zvec, xvec)

    # Step vectors in AU
    xvec = xvec * (resx / AU)
    yvec = yvec * (resy / AU)
    zvec = zvec * (resz / AU)

    # Lower-front-left corner (according to cntrp*)
    tmpvec = (
        lastrM
        - nx * xvec * cntrpx
        - ny * yvec * cntrpy
        - nz * zvec * cntrpz
    )

    points = np.empty((nx, ny, nz), dtype=object)

    for k in range(nz):         # Fortran k = 1..nz
        for j in range(ny):     # Fortran j = 1..ny
            for i in range(nx):  # Fortran i = 1..nx
                rvector = tmpvec + (i + 1) * xvec + (j + 1) * yvec + (k + 1) * zvec
                r = np.linalg.norm(rvector)
                alpha = float(np.arccos(rvector[2] / r))
                beta = float(np.arctan2(rvector[1], rvector[0]))
                points[i, j, k] = Point(
                    r=r,
                    alpha=alpha,
                    beta=beta,
                    rvector=np.asarray(rvector, dtype=np.float64),
                )

    return points


# ----------------------------------------------------------------------
# beta_from_Rg
# ----------------------------------------------------------------------

def beta_from_Rg(Rg: float) -> float:
    """
    Python translation of Fortran beta_from_Rg.

    From the table of grain radii and beta values, interpolate beta
    corresponding to Rg.
    """
    N = 500
    Rgs = np.zeros(N, dtype=np.float64)
    betas = np.zeros(N, dtype=np.float64)

    with open("input_data_files/beta_vs_forsterite_Rg.dat", "r") as f:
        for i in range(N):
            line = f.readline()
            if not line:
                raise RuntimeError(
                    "beta_vs_forsterite_Rg.dat ended prematurely while reading data."
                )
            vals = [float(tok) for tok in line.split()]
            if len(vals) < 2:
                raise RuntimeError(
                    "Each line of beta_vs_forsterite_Rg.dat must contain at least "
                    "two numbers: Rg, beta."
                )
            Rgs[i] = vals[0]
            betas[i] = vals[1]

    return _linterpol(Rgs, betas, Rg)


import math
import os
from typing import Sequence
import numpy as np

# ----------------------------------------------------------------------
# Extra constants (mirror const.f90)
# ----------------------------------------------------------------------

# Gravitational parameter of the Sun in canonical units (AU^3/day^2)
GMsun = 0.0002959122082855908  # 1.327124400419393e20 / AU**3 * 86400**2

# Phaethon radius (used already in get_moving_sources)
RAST_METERS = 2.9e3
Rast_AU = RAST_METERS / AU  # AU was defined earlier in the script


# ----------------------------------------------------------------------
# Runge–Kutta propagator: runge_kutta_point_position
# ----------------------------------------------------------------------

def runge_kutta_point_position(
    r0: Sequence[float],
    v0: Sequence[float],
    mu: float,
    time: float,
) -> np.ndarray:
    """
    Optimized RK4 propagator for a 2-body orbit.

    Parameters
    ----------
    r0 : (3,) initial position (AU)
    v0 : (3,) initial velocity (AU/day)
    mu : scalar, gravitational parameter [AU^3/day^2]
    time : float, integration interval (days; can be negative)

    Returns
    -------
    r : (3,) final position after `time` (AU)
    """
    # Unpack into plain floats (much faster in a tight loop than tiny NumPy arrays)
    x, y, z = map(float, r0)
    vx, vy, vz = map(float, v0)

    # --- step size logic (same as Fortran, just written more clearly) ---
    Nstep = 200
    dt = time / float(Nstep)
    while dt > 3.0e-4:
        Nstep = int(Nstep * 1.2)
        dt = time / float(Nstep)

    dt2 = 0.5 * dt
    dt6 = dt / 6.0

    def accel(x: float, y: float, z: float) -> tuple[float, float, float]:
        """Gravitational acceleration -mu r / |r|^3."""
        r2 = x*x + y*y + z*z
        r = math.sqrt(r2)
        inv_r3 = mu / (r2 * r)  # mu / |r|^3
        ax = -inv_r3 * x
        ay = -inv_r3 * y
        az = -inv_r3 * z
        return ax, ay, az

    for _ in range(Nstep):
        # k1, l1
        k1x, k1y, k1z = accel(x, y, z)
        l1x, l1y, l1z = vx, vy, vz

        # k2, l2
        x2 = x + l1x * dt2
        y2 = y + l1y * dt2
        z2 = z + l1z * dt2
        k2x, k2y, k2z = accel(x2, y2, z2)
        l2x = vx + k1x * dt2
        l2y = vy + k1y * dt2
        l2z = vz + k1z * dt2

        # k3, l3
        x3 = x + l2x * dt2
        y3 = y + l2y * dt2
        z3 = z + l2z * dt2
        k3x, k3y, k3z = accel(x3, y3, z3)
        l3x = vx + k2x * dt2
        l3y = vy + k2y * dt2
        l3z = vz + k2z * dt2

        # k4, l4
        x4 = x + l3x * dt
        y4 = y + l3y * dt
        z4 = z + l3z * dt
        k4x, k4y, k4z = accel(x4, y4, z4)
        l4x = vx + k3x * dt
        l4y = vy + k3y * dt
        l4z = vz + k3z * dt

        # Update v, r
        vx += dt6 * (k1x + 2.0*k2x + 2.0*k3x + k4x)
        vy += dt6 * (k1y + 2.0*k2y + 2.0*k3y + k4y)
        vz += dt6 * (k1z + 2.0*k2z + 2.0*k3z + k4z)

        x  += dt6 * (l1x + 2.0*l2x + 2.0*l3x + l4x)
        y  += dt6 * (l1y + 2.0*l2y + 2.0*l3y + l4y)
        z  += dt6 * (l1z + 2.0*l2z + 2.0*l3z + l4z)

    return np.array([x, y, z], dtype=np.float64)


# ----------------------------------------------------------------------
# Matrix output: matrix_out (data_out.f90)
# ----------------------------------------------------------------------

def matrix_out(fname: str, image: np.ndarray) -> None:
    """
    Approximate Python version of data_out.matrix_out.

    Writes the 2D array `image` to `fname` as text, one row per line.
    """
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    arr = np.asarray(image, dtype=np.float64)
    # Fortran used ES12.4E2; here we use a standard exponential format.
    np.savetxt(fname, arr, fmt="%.4E")


# ----------------------------------------------------------------------
# Main driver: phaethon program
# ----------------------------------------------------------------------

def run_phaethon(
    eph_filename: str = "input_data_files/"
    "Phaethon_2025-02-22_last_int=10min_ECLIPJ2000.dat",
    Neph: int = 1000,
    Nlin: int = 10,
    n1: int = 200,
    n2: int = 200,
    centerpositionx: float = 0.5,
    centerpositiony: float = 0.5,
) -> None:
    Nrgs = 13
    Rgs = np.array(
        [
            0.55,
            0.2,
            0.3,
            0.42,
            0.1,
            0.67,
            0.85,
            1.0,
            1.2,
            2.5,
            4.0,
            6.0,
            10.0,
            99.0
        ],
        dtype=np.float64,
    )

    # Number of points along the asteroid trajectory:
    Nt = (Neph - 1) * Nlin + 1

    # ------------------------------------------------------------------
    # Input source parameters along the orbit
    # ------------------------------------------------------------------
    print("getting moving sources")
    sources, comet = get_moving_sources(eph_filename, Nt, Nlin)
    print("got moving sources")

    # Moment for which we compute the density
    tnow = float(sources[-1].Tj)

    # Resolution of the planar grid [m]
    resolution = np.array([5.0e3, 5.0e3], dtype=np.float64)

    # Build 2D grid of points (n1 × n2)
    points_grid = get_points(
        n1=n1,
        n2=n2,
        resolution=resolution,
        lastrM=np.asarray(comet[-1].coords, dtype=np.float64),
        cntrpx=centerpositionx,
        cntrpy=centerpositiony,
    )

    # Flatten points in the same order as Fortran loops: ii=1..n2, i=1..n1
    points_flat: list[Point] = [
        points_grid[i, j] for j in range(n2) for i in range(n1)
    ]

    # Density arrays
    density = np.zeros((n1, n2), dtype=np.float64)

    def _group_sources_by_ratemap(
        sources: Sequence[Source],
        rhels: np.ndarray,
        idt: int,
        Nt: int,
    ) -> list[tuple[int, int, int, int]]:
        """
        Compute contiguous blocks of i_t that share the same mapind2.

        Returns a list of tuples (block_start, block_end, mapind1, mapind2),
        where block_end is exclusive.
        """
        blocks: list[tuple[int, int, int, int]] = []

        # Initial indices as in your original code
        mapind1 = 0
        mapind2 = 1
        rhel1 = rhels[mapind1]
        rhel2 = rhels[mapind2]

        current_start = idt
        current_mapind1 = mapind1
        current_mapind2 = mapind2

        for i_t in range(idt, Nt - 1):
            r = sources[i_t].r

            # Advance mapind* until this r is between rhel1 and rhel2
            while mapind2 < len(rhels) - 1 and r < rhel2:
                mapind1 = mapind2
                mapind2 += 1
                rhel1 = rhels[mapind1]
                rhel2 = rhels[mapind2]

            # If the pair (mapind1,mapind2) changed, close previous block
            if (mapind1, mapind2) != (current_mapind1, current_mapind2):
                blocks.append((current_start, i_t, current_mapind1, current_mapind2))
                current_start = i_t
                current_mapind1 = mapind1
                current_mapind2 = mapind2

        # Close last block
        blocks.append((current_start, Nt - 1, current_mapind1, current_mapind2))
        return blocks

    # Load list of impact-ejecta maps (Szalay et al. 2019)
    rhels, fnames = get_maps_data()  # length Nmaps=4 in this setup
    

    # ------------------------------------------------------------------
    # Loop over particle radii (different beta and muR)
    # ------------------------------------------------------------------
    for i_R in range(0, Nrgs + 1):  # Fortran: i_R = 0, Nrgs
        density = np.zeros((n1, n2), dtype=float)

        # --- β from grain radius ---
        beta = beta_from_Rg(Rgs[i_R]) if i_R > 0 else 0.0      # ~0 for large grains
        muR = GMsun * (1.0 - beta)

        # --- time-limit dtlim2 ---
        dtlim2 = (
            resolution[0] * n1 * (1.0 - centerpositionx)
            / AU / float(sources[0].ud.umin)
        )

        # --- time-limit dtlim3 ---
        denom = GMsun - muR
        dtlim3 = (
            math.sqrt(
                2.0 * sources[0].r**2 * resolution[0] / AU
                * (1.0 - centerpositionx) * n1 / denom
            )
            if abs(denom) >= 1e-14
            else float("inf")
        )

        dt_limit = min(dtlim2, dtlim3)

        # --- earliest still-visible index idt ---
        # Equivalent to:
        # while idt < Nt-1 and tnow - Tj[idt] > dt_limit: idt += 1
        Tj = np.array([s.Tj for s in sources], dtype=float)
        age = tnow - Tj                                          # time since ejection

        # idt = first index where age <= dt_limit (or Nt-1 if none)
        mask = np.where(age <= dt_limit)[0]
        idt = int(mask[0]) if mask.size else Nt - 1


        # current state of ratemap indices / data
        mapind1 = 0
        mapind2 = 1
        rhel1 = rhels[mapind1]
        rhel2 = rhels[mapind2]

        # Ensure rmap1/rmap2 are in a known initial state
        rhel1 = api.read_first_ratemap(fnames[mapind1])
        rhel2 = api.read_ratemap(fnames[mapind2])

        blocks = _group_sources_by_ratemap(sources, rhels, idt, Nt)

        for block_start, block_end, b_mapind1, b_mapind2 in blocks:
            print(block_start, block_end, b_mapind1, b_mapind2)
            # If we need to move to a new map pair, do the same steps as before
            if b_mapind2 != mapind2:
                # Move mapind1 / mapind2 forward one by one, like original code
                while mapind2 < b_mapind2:
                    mapind1 = mapind2
                    mapind2 += 1

                    # rmap1 = rmap2
                    rmap2 = api.get_rmap2()
                    api.set_rmap1(rmap2)

                    rhel2 = api.read_ratemap(fnames[mapind2])

                rhel1 = rhels[mapind1]
                rhel2 = rhels[mapind2]

            # Now all sources in [block_start, block_end) use same mapind1/mapind2.
            # We can *freeze* the interpolation: one call only.
            api.ratematr_interpolate(
                rhel=rhel2,  # your proposed simplification
                rhel1=rhel1,
                rhel2=rhel2,
            )

            dens_flat = api.batch_over_points_sources(
                points=points_flat,
                sources_by_time=sources[block_start:block_end],
                comets_by_time=comet[block_start:block_end],
                muR=muR,
                tnow=tnow,
                Rast_AU=Rast_AU,
                pericenter=False,          # not used by delta_ejection
                method="simple_expansion",
                )

            # Replace the double loop with a reshape using Fortran order
            density[:, :] += np.asarray(dens_flat, float).reshape(
                (n1, n2), order="F"
            )


        # Output file name: "results/Rg= xx.xxmicron.dat"
        fnameout = f"results/Rg={Rgs[i_R]:5.2f}micron.dat"
        matrix_out(fnameout, density)
        print("result is in the file", fnameout)

        if i_R < Nrgs:
            print("calculations continue")


run_phaethon()