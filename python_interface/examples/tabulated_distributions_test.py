#!/usr/bin/env python3
"""
Minimal sanity test for tabulated ejection speed (fu) and direction (fpsi)
distributed from Python into the Fortran module via the ctypes bridge.

What it does:
1) Builds a normalized power-law PDF for ejection speed between 1 and 200 m/s,
   converted to AU/day.
2) Builds a normalized directional PDF on the sphere:
      f(psi, lambdaM) = (1 + a*P2(cos psi)) / (4*pi)
   which integrates to 1 over the whole sphere.
3) Sends both tables to Fortran via:
      set_tabulated_fu(...)
      set_tabulated_fpsi(...)
4) Runs one scalar density call with ud_shape=10 and eject_distr=10 to ensure
   the code path is exercised (prints returned density).

Adjust import paths if your package name differs.
"""

import numpy as np

# ---- adjust these imports to your package layout ----
from python_interface.dudi_hc.api import (set_tabulated_fu, set_tabulated_fpsi,
    v_integration, delta_ejection, simple_expansion
)
from python_interface.dudi_hc.models import (
    Point, Source, Comet, EjectionSpeedProperties
)
# -----------------------------------------------------


AU_M = 149_597_870_700.0      # meters
DAY_S = 86_400.0              # seconds
MS_TO_AUDAY = DAY_S / AU_M    # (m/s) -> (AU/day)


def powerlaw_pdf(u: np.ndarray, umin: float, umax: float, k: float) -> np.ndarray:
    """
    Normalized power-law pdf on [umin, umax], proportional to u^(-k).
    """
    if not (umin > 0 and umax > umin):
        raise ValueError("Require 0 < umin < umax.")
    if np.any((u < umin) | (u > umax)):
        raise ValueError("u contains values outside [umin, umax].")

    if abs(k - 1.0) < 1e-14:
        # p(u) = C/u, C = 1/ln(umax/umin)
        C = 1.0 / np.log(umax / umin)
        return C / u

    # p(u) = C * u^{-k}, C = (1-k) / (umax^{1-k} - umin^{1-k})
    C = (1.0 - k) / (umax ** (1.0 - k) - umin ** (1.0 - k))
    return C * u ** (-k)


def main() -> None:
    # ============================================================
    # 1) Tabulated ejection speed PDF fu(u) on [1, 200] m/s
    # ============================================================
    umin = 1.0 * MS_TO_AUDAY
    umax = 200.0 * MS_TO_AUDAY
    k = 3.6  # power-law index

    u_tab = np.linspace(umin, umax, 50, dtype=np.float64)  # 50 points on u-interpolation grid
    # Convert grid to AU/day. Note: pdf transforms as p_x(x) = p_u(u) * du/dx.
    fu_tab = powerlaw_pdf(u_tab, umin, umax, k)

    # Quick normalization check in AU/day:
    norm_fu = np.trapezoid(fu_tab, u_tab)
    print(f"[check] fu normalization (AU/day grid): {norm_fu:.12f}")

    set_tabulated_fu(u_tab, fu_tab)
    print("[ok] set_tabulated_fu")

    # ============================================================
    # 2) Tabulated directional PDF fpsi(psi, lambdaM)
    #    f = (1 + a*P2(cos psi)) / (4*pi), integral over sphere = 1
    # ============================================================
    a = 0.5  # keep |a| < 1 for non-negativity

    Npsi = 181
    # Use a lambdaM grid that does NOT span [0, 2pi] to exercise the "outside" case logic.
    Nlam = 180
    psi_tab = np.linspace(0.0, np.pi, Npsi, dtype=np.float64)
    lambdaM_tab = np.linspace(0.1, 6.2, Nlam, dtype=np.float64)  # subset of [0, 2pi]

    # P2(x) = (1/2)*(3x^2 - 1)
    cospsi = np.cos(psi_tab)
    P2 = 0.5 * (3.0 * cospsi**2 - 1.0)

    # Axisymmetric in lambdaM (independent of lambdaM)
    fpsi_1d = (1.0 + a * P2) / (4.0 * np.pi)  # shape (Npsi,)
    fpsi_tab = np.repeat(fpsi_1d[:, None], Nlam, axis=1)  # shape (Npsi, Nlam)

    # Check normalization numerically over full sphere:
    # integral = ∫_0^{2pi} ∫_0^{pi} f(psi)*sin(psi) dpsi dlam
    sphere_norm = (2.0 * np.pi) * np.trapz(fpsi_1d * np.sin(psi_tab), psi_tab)
    print(f"[check] fpsi normalization over sphere: {sphere_norm:.12f}")

    set_tabulated_fpsi(psi_tab, lambdaM_tab, fpsi_tab)
    print("[ok] set_tabulated_fpsi")

    # ============================================================
    # 3) One scalar call that uses ud_shape=10 and eject_distr=10
    #    to ensure the tabulated branches are hit.
    # ============================================================
    # Simple geometry: point and source near 1 AU, same direction.
    pt = Point(r=1.0, alpha=0.0, beta=0.0, rvector=np.array([1.0, 0.0, 0.0], float))
    src = Source(
        r=1.0, alphaM=0.0, betaM=0.0,
        rrM=np.array([1.0, 0.0, 0.0], float),
        zeta=0.0, eta=0.0,
        symmetry_axis=np.array([0.0, 0.0, 1.0], float),
        ejection_angle_distr=10,     # 10 is the code for the tabulated distribution
        ud=EjectionSpeedProperties(ud_shape=10, umin=0.0, umax=1e-3),   # 10 is the code for the tabulated distribution
        Nparticles=1.0e10, Tj=0.0, dtau=1.0e-4,
    )
    cm = Comet(coords=np.zeros(3, float), Vastvec=np.zeros(3, float), Vast=0.0)

    muR = 1.0e-4      
    tnow = 2.0
    Rast_AU = 0.0
    pericenter = False

    dens = v_integration(pt, src, cm, muR=1e-4, tnow=1.0, Rast_AU=1e-6, pericenter=False)
    print(f"[result] density (v_integration, tabulated fu/fpsi): {dens:.6e}")


if __name__ == "__main__":
    main()
