"""
Sample the model on a small angular grid and report min/max density.

Run:
  python python_interface/examples/grid_sample.py
"""
import numpy as np
from python_interface.dudi_hc.api import simple_expansion
from python_interface.dudi_hc.models import (
    Point, Source, EjectionSpeedProperties
)

def make_source():
    return Source(
        r=1.0, alphaM=0.0, betaM=0.0,
        rrM=np.array([1.0, 0.0, 0.0], float),
        zeta=0.0, eta=0.0,
        symmetry_axis=np.array([0.0, 0.0, 1.0], float),
        ejection_angle_distr=0,
        ud=EjectionSpeedProperties(ud_shape=0, umin=0.0, umax=1.0),
        Nparticles=1.0e10, Tj=0.0, dtau=1.0e-4,
    )

def main():
    src = make_source()
    cloudcentr = np.array([0.0, 0.0, 0.0], float)
    dt = 1e-3  # small positive to avoid any zero-division internals

    alphas = np.linspace(-0.2, 0.2, 9)
    betas  = np.linspace(-0.2, 0.2, 9)

    vals = []
    for a in alphas:
        for b in betas:
            pt = Point(r=1.0, alpha=a, beta=b, rvector=np.array([1.0, 0.0, 0.0], float))
            vals.append(simple_expansion(pt, src, cloudcentr=cloudcentr, dt=dt))

    vals = np.array(vals, float)
    print(f"Grid size: {vals.size}")
    print(f"min={vals.min():.6g}, max={vals.max():.6g}, mean={vals.mean():.6g}")

if __name__ == "__main__":
    main()
