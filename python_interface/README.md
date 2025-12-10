# DUDI-HC Python Interface

## 1. What is this?

This directory provides the **Python interface** 
to the **DUDI-heliocentric (DUDI-HC)** model for dust number density calculations
around comets, asteroids, and atmosphereless bodies.  
The original DUDI-HC implementation is entirely in Fortran-95, and no 
changes were made to the physical or numerical core for this Python release.

To support calling DUDI-HC directly from Python, a small Fortran-2003 
layer is added using:

- `bind(C)` interfaces,
- C-compatible derived types,
- a compiled shared library exposed to Python through `ctypes`.

This Python interface provides:

- A documented and object-oriented high-level API 
  (`Point`, `Source`, `Comet`, `v_integration`, batching routines).
- Full access to all scientifically relevant features of DUDI-HC without writing
  in Fortran.
- 1-to-1 Python equivalents of the Fortran example programs 
  (select_method, example, Phaethon).


## 2. Architecture Overview

Python (your scripts, notebooks)
        |
        v      high-level API + data models + batching utilities
  python_interface/dudi_hc/api.py
  python_interface/dudi_hc/models.py
        |
        v      thin ctypes wrapper (NumPy <-> raw C arrays)
  python_interface/dudi_hc/ctypes_bridge.py
        |
        v      shared library exposing Fortran routines as C functions
  python_interface/dudi_hc/libpy_dudihc_bridge.so
        |
        v      Fortran-2003 wrapper module (bind(C))
  python_interface/fortran_bridge/py_bridge.f90 and helpers
        |
        v      Fortran-95 DUDI-HC scientific core (unchanged)


Python never touches Fortran derived types directly. The C-bindable Fortran 
wrappers (py_bridge.f90) receive only C-friendly scalars/arrays, reconstruct 
the Fortran types (position_in_space, source_properties, ephemeris), call the 
real routines in DUDIhc, and return a single real(8) density value to Python. 
**Types & precision**
The DUDIhc routines compute density as real(4). Our wrappers convert that 
to real(8) (double) just before returning so Python gets a normal float. 
Python vectors are passed as NumPy float64 contiguous arrays of shape (3,). 
Fortran LOGICAL inputs at the boundary are passed as C integers (0/1).

All computational kernels and all **OpenMP parallelism** remain entirely in 
**Fortran**, where they are most efficient. The Python layer is thin and 
designed only for data preparation, convenience, and analysis.


## 3. Installation and Build (Linux)

The model has been tested and validated on **Linux** systems. 
Source installation requires a working Fortran toolchain.

### 3.1. System prerequisites

Install:

- `gfortran` (Fortran 95/2003 compiler)
- OpenMP (`libgomp` – normally included with GCC)
- Python ≥ 3.9
- `pip`
- GitHub CLI or `git`

Example (Ubuntu-like):

sudo apt-get install gfortran libgomp1 python3 python3-pip git

### 3.2. Clone the repository

git clone https://github.com/Veyza/dudi-heliocentric.git
cd dudi-heliocentric

### 3.3. Build the Fortran ctypes bridge

bash python_interface/fortran_bridge/build_ctypes_bridge.sh

When finished, you should see:

python_interface/dudi_hc/libpy_dudihc_bridge.so

### 3.4. Install the Python package

Editable mode:

pip install -e .

### Tests

python3 -m pytest -q
