import numpy as np

from python_interface.dudi_hc.models import (
    Point, Source, Comet,
    EjectionSpeedProperties,
    spherical_to_cartesian,
    normalize,  # optional, just to demo axis normalization
)

def main() -> None:
    # --- Build Point ---
    r, alpha, beta = 1.0, 1.0, 0.5  # AU, rad, rad
    rvector = spherical_to_cartesian(r, alpha, beta)  # AU
    p = Point(r=r, alpha=alpha, beta=beta, rvector=rvector)
    print("Point OK:", p)
    print("  rvector:", p.rvector, "shape:", p.rvector.shape)

    # --- Build Source ---
    rM, alphaM, betaM = 1.0, 1.0, 0.0
    rrM = spherical_to_cartesian(rM, alphaM, betaM)
    # You decide how to construct symmetry_axis; here we just normalize an example:
    sym_axis = normalize(np.array([0.1, 0.2, 0.97], dtype=np.float64))

    ud = EjectionSpeedProperties(ud_shape=1, umin=0.0, umax=0.01)  # AU/day
    s = Source(
        r=rM, alphaM=alphaM, betaM=betaM, rrM=rrM,
        zeta=0.3, eta=1.2,
        symmetry_axis=sym_axis,
        ejection_angle_distr=3,
        ud=ud,
    )
    print("Source OK:", s)
    print("  rrM:", s.rrM, "shape:", s.rrM.shape)
    print("  symmetry_axis norm:", np.linalg.norm(s.symmetry_axis))

    # --- Build Comet ---
    coords = np.array([1.0, 0.0, 0.0], dtype=np.float64)        # AU
    Vastvec = np.array([0.0, 1.0, 0.0], dtype=np.float64)       # AU/day
    Vast = float(np.linalg.norm(Vastvec))                       # AU/day
    c = Comet(coords=coords, Vastvec=Vastvec, Vast=Vast)
    print("Comet OK:", c)
    print("  coords:", c.coords, "shape:", c.coords.shape)
    print("  Vastvec:", c.Vastvec, "||V||:", np.linalg.norm(c.Vastvec), "Vast:", c.Vast)

    print("\nAll model constructions passed.")

if __name__ == "__main__":
    main()
