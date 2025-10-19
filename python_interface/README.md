# Python Interface for DUDI-heliocentric
This optional interface provides a lightweight Python wrapper around the
Fortran core routines of DUDI-heliocentric, allowing users to run and
visualize results from Python without modifying the original Fortran code.

Python ↔ Fortran bridge (how it works & how to build)
Architecture at a glance
python_interface/dudi_hc/api.py     ← (Python-facing thin API; calls the bridge)
        │
        ├── imports
        ▼
python_interface/dudi_hc/_bridge_ctypes.py
        │  (loads one shared lib with ctypes and exposes 3 functions)
        │
        ├── ctypes.CDLL("libpy_dudihc_bridge.so")
        ▼
python_interface/fortran_bridge/py_bridge.f90
   (Fortran wrappers with `bind(C)`; rebuild derived types, call DUDIhc)
        │
        ▼
src/*.f90 (Fortran core: DUDIhc.f90, define_types.f90, …)

Key idea: Python never touches Fortran derived types directly. The C-bindable 
Fortran wrappers (py_bridge.f90) receive only C-friendly scalars/arrays, 
reconstruct the Fortran types (position_in_space, source_properties, ephemeris),
call the real routines in DUDIhc, and return a single real(8) density value to Python.
**Types & precision**
*The DUDIhc routines compute density as real(4). Our wrappers convert that 
to real(8) (double) just before returning so Python gets a normal float.
*Python vectors are passed as NumPy float64 contiguous arrays of shape (3,).
*Fortran LOGICAL inputs at the boundary are passed as C integers (0/1).

# Build (one command)
**Re-run this script any time you change .f90 sources or the wrapper.**
bash python_interface/fortran_bridge/build_ctypes_bridge.sh
What it does:
*Compiles the Fortran core in dependency order (src/*.f90) into objects.
*Compiles the C-bind wrapper py_bridge.f90.
*Links everything into one shared library:
"python_interface/dudi_hc/libpy_dudihc_bridge.so"
*Smoke tests that the ctypes loader sees the functions.

# Using it from Python
1. Low-level bridge:

from python_interface.dudi_hc import _bridge_ctypes as fb
density = fb.call_v_integration(
    point_r=..., point_alpha=..., point_beta=..., point_rvector=[...,...,...],
    src_r=..., src_alphaM=..., src_betaM=..., src_rrM=[...,...,...],
    src_zeta=..., src_eta=..., src_axis=[...,...,...],
    src_eject_distr=..., src_ud_shape=..., src_umin=..., src_umax=...,
    comet_coords=[...,...,...], comet_vastvec=[...,...,...], comet_vast=...,
    muR=..., tnow=..., Rast_AU=..., pericenter=True
)
print(density)  # Python float

(There are also call_delta_ejection(...) and call_simple_expansion(...) 
with the obvious signatures.)

2. High-level API:

from python_interface.dudi_hc.api import v_integration, Point, Source, Comet
#construct Point/Source/Comet, then:
density = v_integration(point, source, comet, muR=..., tnow=..., Rast_AU=..., pericenter=True)

# File locations you’ll care about

- src/ — Fortran core (unchanged scientific code)
- python_interface/fortran_bridge/py_bridge.f90 — C-bindable wrappers
- python_interface/fortran_bridge/build_ctypes_bridge.sh — rerunnable build script
- python_interface/dudi_hc/_bridge_ctypes.py — Python ctypes loader
- python_interface/dudi_hc/libpy_dudihc_bridge.so — built shared library (not tracked by git)

## Install / Build

### Prereqs (for maintainers)
- gfortran, OpenMP runtime
- Python 3.9+ and pip

### Build Fortran ctypes bridge

bash python_interface/fortran_bridge/build_ctypes_bridge.sh

### Install (editable dev mode)

pip install -e .

### Tests

python3 -m pytest -q

### Build a wheel + sdist

python3 -m pip install build
python3 -m build
The wheel bundles the compiled shared libraries, so end users can pip install 
without a Fortran compiler (on compatible platforms).

## Quick start (Python API)

1. Build the Fortran ctypes bridge (compiles Fortran and installs `.so` into `python_interface/dudi_hc/`):
bash python_interface/fortran_bridge/build_ctypes_bridge.sh
2. Install (editable dev mode is convenient while developing):
pip install -e .
3. Run examples:
python examples/minimal.py
python examples/grid_sample.py

## Troubleshooting

**ImportError: undefined symbol**
Rebuild the bridge so the .so matches the current Fortran sources:
bash python_interface/fortran_bridge/build_ctypes_bridge.sh
**No prints from Fortran DIAG mode**
Use the environment flag and ensure Fortran output is flushed (already wired in):
export HC_BRIDGE_DIAG=1
python examples/minimal.py
unset HC_BRIDGE_DIAG
**ABI issues / segfaults**
We fixed the scalar calling convention by using value in the bind(C) wrappers.
If you edit py_bridge.f90, rebuild and re-run.
