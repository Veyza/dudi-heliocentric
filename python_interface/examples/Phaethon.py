from __future__ import annotations

from typing import Iterable, Sequence
import numpy as np
import types

from python_interface.dudi_hc.models import Point, Source, Comet
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


def _norma3d(v: Iterable[float]) -> float:
    """Euclidean norm of a 3-vector (Fortran norma3d)."""
    arr = np.asarray(v, dtype=np.float64)
    if arr.shape != (3,):
        raise ValueError(f"Expected 3-vector, got shape {arr.shape}")
    return float(np.linalg.norm(arr))


def _vector_product(a: Iterable[float], b: Iterable[float]) -> np.ndarray:
    """3-D vector product (Fortran vector_product)."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != (3,) or b.shape != (3,):
        raise ValueError("vector_product expects two 3-vectors")
    return np.cross(a, b)


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
    sources: Sequence[Source],
    comet: Sequence[Comet],
) -> None:
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
    sources : sequence of Source
        Pre-allocated sequence of length Np, modified in place.
    comet : sequence of Comet
        Pre-allocated sequence of length Np, modified in place.
    """
    if len(sources) != Np or len(comet) != Np:
        raise ValueError("sources and comet must both have length Np")

    # ----- map file names and distances -----
    rhels, fnames = get_maps_data()
    if len(fnames) < 2:
        raise RuntimeError("Need at least two ratemap files.")

    # Fortran: mapind1=1, mapind2=mapind1+1 (1-based)
    mapind1 = 0
    mapind2 = 1

    # Use api wrappers to call Fortran:
    rhel1 = api.read_first_ratemap(fnames[mapind1])
    rhel2 = api.read_ratemap(fnames[mapind2])

    def _copy_rmap2_to_rmap1() -> None:
        """Fortran 'rmap1 = rmap2' via api."""
        rmap2 = api.get_rmap2()
        api.set_rmap1(rmap2)

    # ----- read ephemeris and interpolate -----
    moment = np.zeros(Np, dtype=np.float64)
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
        comet[0].coords = np.array(vals[1:4], dtype=np.float64)
        comet[0].Vastvec = np.array(vals[4:7], dtype=np.float64)

        # Fortran: do i = Nlin+1, Np, Nlin  (1-based)
        # Python index: i = i_fortran - 1
        for i in range(Nlin, Np, Nlin):
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

            moment[i] = vals[0]
            comet[i].coords = np.array(vals[1:4], dtype=np.float64)
            comet[i].Vastvec = np.array(vals[4:7], dtype=np.float64)

            if Nlin > 1:
                # Fortran: forall(ii = (i-Nlin+1):(i-1))
                for j in range(i - Nlin + 1, i):
                    # weight w = (ii - i + Nlin) / dNlin
                    w = float(j - i + Nlin) / dNlin

                    moment[j] = moment[i - Nlin] + (moment[i] - moment[i - Nlin]) * w
                    comet[j].Vastvec = (
                        comet[i - Nlin].Vastvec
                        + (comet[i].Vastvec - comet[i - Nlin].Vastvec) * w
                    )
                    comet[j].coords = (
                        comet[i - Nlin].coords
                        + (comet[i].coords - comet[i - Nlin].coords) * w
                    )

    # Shift timeline: moment(i) = moment(i) - moment(1)
    moment -= moment[0]

    # Phaethon radius, as in Fortran
    Rast = 2.9e3  # [m]

    # ----- fill source and comet properties for each time step -----
    for i in range(Np):
        rrM = np.asarray(comet[i].coords, dtype=np.float64)
        sources[i].rrM = rrM
        r = _norma3d(rrM)
        sources[i].r = r

        # Update which maps we interpolate between when r drops below rhel2
        if r < rhel2 and mapind2 + 1 < len(fnames):
            rhel1 = rhel2
            mapind1 = mapind2
            mapind2 += 1
            _copy_rmap2_to_rmap1()
            rhel2 = api.read_ratemap(fnames[mapind2])

        # Interpolate ratemap for this r (Fortran ratematr_interpolate)
        api.ratematr_interpolate(rhel=r, rhel1=rhel1, rhel2=rhel2)

        # Asteroid speed at position i
        comet[i].Vast = _norma3d(comet[i].Vastvec)

        # Angular coordinates
        sources[i].alphaM = float(np.arccos(rrM[2] / r))
        sources[i].betaM = float(np.arctan2(rrM[1], rrM[0]))

        # Symmetry axis: Fortran assigns a scalar to an array, broadcasting;
        # we mimic that behaviour exactly.
        axis_val = rrM[0] / r
        sources[i].symmetry_axis = np.array(
            [axis_val, axis_val, axis_val], dtype=np.float64
        )

        # Distributions / timing parameters
        sources[i].zeta = 0.0
        sources[i].eta = 0.0
        sources[i].ud.ud_shape = 1
        sources[i].ud.umin = 2.0 / AUdays2SI
        sources[i].ud.umax = 2399.0 / AUdays2SI
        sources[i].ejection_angle_distr = 3
        sources[i].Tj = float(moment[i])
        sources[i].dtau = 0.0

        # Integrate number density of impact ejecta over the map
        totrate = integrate_over_matrix()

        # Convert number density to flux (Szalay et al. 2016, Eq. 3)
        totrate = totrate / 0.31 / 7.2e-3 / 4.0 / PI * Rast**2

        # Convert flux to number of ejected particles
        dt_days = moment[1] - moment[0] if Np > 1 else 0.0
        sources[i].Nparticles = totrate * dt_days * S_IN_DAY


# ----------------------------------------------------------------------
# get_flyby_trajectory
# ----------------------------------------------------------------------

def get_flyby_trajectory(
    nt1: int,
    resolution: float,
    CAdist: float,
    lastrM: np.ndarray,
) -> list[Point]:
    """
    Python translation of Fortran get_flyby_trajectory.

    Parameters
    ----------
    nt1 : int
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
    xvec = lastrM / _norma3d(lastrM)
    zvec = zvec - zvec * float(np.dot(xvec, zvec))
    zvec = zvec / _norma3d(zvec)
    yvec = _vector_product(zvec, xvec)

    angle2xvec = 29.0 * DEG2RAD

    tmpvec = xvec * np.cos(angle2xvec) + yvec * np.sin(angle2xvec)
    tmpvec = tmpvec * (resolution / AU)      # m -> AU
    tmpnorm = -xvec * np.sin(angle2xvec) + yvec * np.cos(angle2xvec)

    CApoint = lastrM - tmpnorm * CAdist     # CAdist already in AU

    points: list[Point] = []
    for i in range(1, nt1 + 1):
        rvector = CApoint + tmpvec * (i - nt1 / 2.0)
        r = _norma3d(rvector)
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
    Python translation of Fortran integrate_over_matrix.

    Integrates the current Fortran ratemap (already set via api.read_* and
    api.ratematr_interpolate) over the sphere.
    """
    ratemap = api.get_ratemap()   # shape (nlats, nlons)
    lats = api.get_lats()         # shape (nlats,)
    lons = api.get_lons()         # shape (nlons,)

    nlats, nlons = ratemap.shape
    if lats.shape[0] != nlats or lons.shape[0] != nlons:
        raise RuntimeError("Inconsistent ratemap / lats / lons dimensions")

    integral = 0.0

    # Polar angle (colatitude) = π/2 - latitude
    polangle = HALFPI - lats

    # Main bands: Fortran ii = 2..nlats
    for ii in range(1, nlats):
        # ring at latitude index ii
        rint = 0.0
        for i in range(1, nlons):
            rint += (ratemap[ii, i] + ratemap[ii, i - 1]) * (lons[i] - lons[i - 1])
        # close ring
        rint += (ratemap[ii, 0] + ratemap[ii, nlons - 1]) * (
            lons[0] - lons[nlons - 1] + TWOPI
        )
        rint *= 0.5

        # ring at previous latitude (ii-1)
        rint1 = 0.0
        for i in range(1, nlons):
            rint1 += (ratemap[ii - 1, i] + ratemap[ii - 1, i - 1]) * (
                lons[i] - lons[i - 1]
            )
        rint1 += (ratemap[ii - 1, 0] + ratemap[ii - 1, nlons - 1]) * (
            lons[0] - lons[nlons - 1] + TWOPI
        )
        rint1 *= 0.5

        # trapezoid in colatitude with sin(polangle)
        integral += (rint + rint1) * np.sin(polangle[ii]) * (
            polangle[ii - 1] - polangle[ii]
        ) * 0.5

    # Uppermost ring (Fortran index 1)
    rint = 0.0
    for i in range(1, nlons):
        rint += (ratemap[0, i] + ratemap[0, i - 1]) * (lons[i] - lons[i - 1])
    rint += (ratemap[0, 0] + ratemap[0, nlons - 1]) * (
        lons[0] - lons[nlons - 1] + TWOPI
    )
    rint *= 0.5
    integral += (
        rint
        * np.sin((-HALFPI - lats[0]) / 2.0)
        * (-HALFPI - lats[0])
    )

    # Lowermost ring (Fortran index nlats)
    rint = 0.0
    for i in range(1, nlons):
        rint += (ratemap[nlats - 1, i] + ratemap[nlats - 1, i - 1]) * (
            lons[i] - lons[i - 1]
        )
    rint += (ratemap[nlats - 1, 0] + ratemap[nlats - 1, nlons - 1]) * (
        lons[0] - lons[nlons - 1] + TWOPI
    )
    integral += (
        rint
        * np.sin((HALFPI - lats[nlats - 1]) / 2.0)
        * (HALFPI - lats[nlats - 1])
    )

    return float(integral)


# ----------------------------------------------------------------------
# get_points (2-D grid)
# ----------------------------------------------------------------------

def get_points(
    nt1: int,
    nt2: int,
    resolution: tuple[float, float],
    lastrM: np.ndarray,
    cntrpx: float,
    cntrpy: float,
) -> np.ndarray:
    """
    Python translation of Fortran get_points.

    Returns a (nt1, nt2) array of Point objects.
    """
    resx, resy = resolution
    lastrM = np.asarray(lastrM, dtype=np.float64)
    if lastrM.shape != (3,):
        raise ValueError("lastrM must be a 3-vector in AU")

    # CS as in Szalay et al. 2019 comparison
    zvec = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    xvec = lastrM / _norma3d(lastrM)
    zvec = zvec - zvec * float(np.dot(xvec, zvec))
    zvec = zvec / _norma3d(zvec)
    yvec = _vector_product(zvec, xvec)

    xvec = xvec * (resx / AU)
    yvec = yvec * (resy / AU)
    tmpvec = lastrM - nt1 * xvec * cntrpx - nt2 * yvec * cntrpy

    points = np.empty((nt1, nt2), dtype=object)

    for j in range(nt2):      # Fortran ii = 1..nt2
        for i in range(nt1):  # Fortran i  = 1..nt1
            rvector = tmpvec + (i + 1) * xvec + (j + 1) * yvec
            r = _norma3d(rvector)
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
    if _norma3d(xvec) == 0.0:
        xvec = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    else:
        xvec = xvec / _norma3d(xvec)

    # Make ẑ ⟂ x̂, then ŷ = ẑ × x̂
    zvec = zvec - xvec * float(np.dot(xvec, zvec))
    zvec = zvec / _norma3d(zvec)
    yvec = _vector_product(zvec, xvec)

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
                r = _norma3d(rvector)
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
    Python translation of help.f90: runge_kutta_point_position.

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
    r = np.asarray(r0, dtype=np.float64).copy()
    v = np.asarray(v0, dtype=np.float64).copy()

    Nstep = 200
    dt = time / float(Nstep)

    # Same adaptive logic as in Fortran
    while dt > 3.0e-4:
        Nstep = int(Nstep * 1.2)
        dt = time / float(Nstep)

    for _ in range(Nstep):
        r2 = np.sum(r * r)
        k1 = -mu / math.sqrt(r2) ** 3 * r
        l1 = v

        r_tmp = r + l1 * dt / 2.0
        r2_tmp = np.sum(r_tmp * r_tmp)
        k2 = -mu / math.sqrt(r2_tmp) ** 3 * r_tmp
        l2 = v + k1 * dt / 2.0

        r_tmp = r + l2 * dt / 2.0
        r2_tmp = np.sum(r_tmp * r_tmp)
        k3 = -mu / math.sqrt(r2_tmp) ** 3 * r_tmp
        l3 = v + k2 * dt / 2.0

        r_tmp = r + l3 * dt
        r2_tmp = np.sum(r_tmp * r_tmp)
        k4 = -mu / math.sqrt(r2_tmp) ** 3 * r_tmp
        l4 = v + k3 * dt

        v = v + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        r = r + dt / 6.0 * (l1 + 2.0 * l2 + 2.0 * l3 + l4)

    return r


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
    Neph: int = 2000,
    Nlin: int = 10,
    nt1: int = 400,
    nt2: int = 400,
    centerpositionx: float = 0.5,
    centerpositiony: float = 0.5,
) -> None:
    Nrgs = 13
    Rgs = np.array(
        [
            0.1,
            0.2,
            0.3,
            0.42,
            0.55,
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
    Np = (Neph - 1) * Nlin + 1

    # ------------------------------------------------------------------
    # Allocate objects: Source and Comet
    # ------------------------------------------------------------------
    # Source dataclass requires:
    #   r, alphaM, betaM, rrM, zeta, eta, symmetry_axis,
    #   ejection_angle_distr, ud, Nparticles, Tj, dtau
    #
    # We create "empty" sources with sensible zero defaults; they will
    # be filled in by get_moving_sources.
    def _make_empty_source() -> Source:
        return Source(
            r=0.0,
            alphaM=0.0,
            betaM=0.0,
            rrM=np.zeros(3, dtype=np.float64),
            zeta=0.0,
            eta=0.0,
            symmetry_axis=np.array([1.0, 0.0, 0.0], dtype=np.float64),
            ejection_angle_distr=0,
            # ud is whatever type you used in models; a SimpleNamespace with
            # the right attributes works fine because api.py only accesses
            # ud.ud_shape, ud.umin, ud.umax.
            ud=types.SimpleNamespace(
                ud_shape=0,
                umin=0.0,
                umax=0.0,
            ),
            Nparticles=0.0,
            Tj=0.0,
            dtau=0.0,
        )

    # Comet dataclass (mirror of ephemeris type):
    #   coords(3), Vastvec(3), Vast
    def _make_empty_comet() -> Comet:
        return Comet(
            coords=np.zeros(3, dtype=np.float64),
            Vastvec=np.zeros(3, dtype=np.float64),
            Vast=0.0,
        )

    sources: list[Source] = [_make_empty_source() for _ in range(Np)]
    comet: list[Comet] = [_make_empty_comet() for _ in range(Np)]

    # ------------------------------------------------------------------
    # Input source parameters along the orbit
    # ------------------------------------------------------------------
    get_moving_sources(eph_filename, Np, Nlin, sources, comet)

    # Moment for which we compute the density
    tnow = float(sources[-1].Tj)

    # Resolution of the planar grid [m]
    resolution = np.array([5.0e3, 5.0e3], dtype=np.float64)

    # Build 2D grid of points (nt1 × nt2)
    points_grid = get_points(
        nt1=nt1,
        nt2=nt2,
        resolution=resolution,
        lastrM=np.asarray(comet[-1].coords, dtype=np.float64),
        cntrpx=centerpositionx,
        cntrpy=centerpositiony,
    )

    # Flatten points in the same order as Fortran loops: ii=1..nt2, i=1..nt1
    points_flat: list[Point] = [
        points_grid[i, j] for j in range(nt2) for i in range(nt1)
    ]

    # Load list of impact-ejecta maps (Szalay et al. 2019)
    rhels, fnames = get_maps_data()  # length Nmaps=4 in this setup

    # Density arrays
    density = np.zeros((nt1, nt2), dtype=np.float64)

    # ------------------------------------------------------------------
    # Loop over particle radii (different beta and muR)
    # ------------------------------------------------------------------
    for i_R in range(0, Nrgs + 1):  # Fortran: i_R = 0, Nrgs
        # Reset which ejecta maps we interpolate between
        mapind1 = 0  # Fortran index 1
        mapind2 = 1  # Fortran index 2

        rhel1 = api.read_first_ratemap(fnames[mapind1])  # rmap1
        rhel2 = api.read_ratemap(fnames[mapind2])        # rmap2

        density[:, :] = 0.0

        # beta from grain radius
        if i_R > 0:
            beta = beta_from_Rg(Rgs[i_R])
        else:
            beta = 0.0  # large grains ~100 μm

        muR = GMsun * (1.0 - beta)

        # Time limits (in days) for contributing ejecta
        dtlim2 = (
            resolution[0]
            * nt1
            * (1.0 - centerpositionx)
            / AU
            / float(sources[0].ud.umin)
        )

        denom = GMsun - muR
        if abs(denom) < 1e-14:
            dtlim3 = float("inf")
        else:
            dtlim3 = math.sqrt(
                2.0
                * sources[0].r**2
                * resolution[0]
                / AU
                * (1.0 - centerpositionx)
                * nt1
                / denom
            )

        dt_limit = min(dtlim2, dtlim3)

        # Find earliest index whose dust is still in field of view
        idt = 0  # Python 0-based; Fortran started from 1
        while (
            idt < Np - 1
            and tnow - float(sources[idt].Tj) > dt_limit
        ):
            idt += 1

        print("start index", idt + 1)  # report in Fortran-style 1-based

        # ------------------------------------------------------------------
        # Loop over active sources along the trajectory
        # ------------------------------------------------------------------
        for i_p in range(idt, Np - 1):  # Fortran: i_p = idt, Np-1
            # Choose which ejecta maps correspond to current heliocentric distance
            while (
                mapind2 < len(fnames) - 1
                and sources[i_p].r < rhel2
            ):
                rhel1 = rhel2
                mapind1 = mapind2
                mapind2 += 1

                # rmap1 = rmap2
                rmap2 = api.get_rmap2()
                api.set_rmap1(rmap2)

                rhel2 = api.read_ratemap(fnames[mapind2])

            # Interpolate ratemap for this heliocentric distance
            api.ratematr_interpolate(
                rhel=sources[i_p].r,
                rhel1=rhel1,
                rhel2=rhel2,
            )

            dt = tnow - float(sources[i_p].Tj)

            # Cloud centre (propagate along two-body orbit with muR)
            cloudcentr = runge_kutta_point_position(
                r0=comet[i_p].coords,
                v0=comet[i_p].Vastvec,
                mu=muR,
                time=dt,
            )

            # Set rMtmp in Fortran (used in distributions_fun)
            api.set_rMtmp(cloudcentr)

            # Compute densities at all points for this source
            dens_flat = api.batch_over_points(
                points=points_flat,
                source=sources[i_p],
                comet=comet[i_p],
                muR=muR,
                tnow=tnow,
                dt=dt,
                Rast_AU=Rast_AU,
                pericenter=False,
                cloudcentr=cloudcentr,
                method="simple_expansion",
            )

            # Accumulate into 2D array, matching Fortran (i,ii) ordering
            idx = 0
            for j in range(nt2):
                for i in range(nt1):
                    density[i, j] += float(dens_flat[idx])
                    idx += 1

        # Output file name: "results/Rg= xx.xxmicron.dat"
        fnameout = f"results/Rg={Rgs[i_R]:5.2f}micron.dat"
        matrix_out(fnameout, density)
        print("result is in the file", fnameout)

        if i_R < Nrgs:
            print("calculations continue")


run_phaethon()