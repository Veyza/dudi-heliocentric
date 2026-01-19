"""
Minimal, end-to-end example that calls all three routines and prints results.

Prereqs:
  1) Build the Fortran ctypes bridge:
       bash python_interface/fortran_bridge/build_ctypes_bridge.sh
  2) (optional) Install in editable mode:
       pip install -e .

Run:
  python python_interface/examples/minimal.py
"""
import numpy as np

from python_interface.dudi_hc.api import (
    v_integration, delta_ejection, simple_expansion
)
from python_interface.dudi_hc.models import (
    Point, Source, Comet, EjectionSpeedProperties
)

pt = Point(r=1.0, alpha=0.0, beta=0.0, rvector=np.array([1.0, 0.0, 0.0], float))
src = Source(
    r=1.0, alphaM=0.0, betaM=0.0,
    rrM=np.array([1.0, 0.0, 0.0], float),
    zeta=0.0, eta=0.0,
    symmetry_axis=np.array([0.0, 0.0, 1.0], float),
    ejection_angle_distr=0,
    ud=EjectionSpeedProperties(ud_shape=0, umin=0.0, umax=1e-3),
    Nparticles=1.0e10, Tj=0.0, dtau=1.0e-4,
)
cm = Comet(coords=np.zeros(3, float), Vastvec=np.zeros(3, float), Vast=0.0)

print("simple_expansion:", simple_expansion(pt, src, cloudcentr=[1.0, 0.0, 0.0], dt=1.0))
print("delta_ejection:",  delta_ejection(pt, src, cm, muR=1e-4, dt=1.0, Rast_AU=1e-6))
print("v_integration:",   v_integration(pt, src, cm, muR=1e-4, tnow=1.0, Rast_AU=1e-6, pericenter=False))
