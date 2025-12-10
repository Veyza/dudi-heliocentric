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

## 4. Model input and solution method selection

For the physical meaning of all input quantities, users should refer to 
**Section 3 “Supplying Input Data” of the main package README**, which describes 
the required model parameters, units, and conventions used by DUDI-HC. The Python 
classes defined in python_interface/dudi_hc/models (Point, Source, Comet) 
are direct mirrors of the Fortran derived types defined in the DUDI-HC core and 
follow the same data structure and validation logic. Likewise, the high-level 
API functions exposed in dudi_hc.api use the same notation and argument 
structure as the corresponding Fortran routines in DUDIhc.f90, enabling users to
rely on the main README and the original Fortran documentation when preparing 
input for the Python interface.

For the explanation of the three numerical methods implemented in DUDI-HC 
(simple expansion, delta-ejection, and v-integration) users should refer to 
**Section 4 of the main README** and to the original A&A paper, which provide  
physical and mathematical context. The Python interface additionally provides 
a direct analogue of the Fortran routine **select_method**, enabling users to 
evaluate the applicability of each method from within Python.

## 5. Python examples

The Python examples can be found in the folder python_interface/examples and may
serve as templates for your own applications. Apart from `minimal.py`, which is 
a small Python-only demonstration, each script is a direct translation of the
corresponding Fortran example and reproduces the same physical setup and output
structure.

### `select_method.py`

`select_method.py` is the Python equivalent of the Fortran program
`examples/select_method.f90`. It reads orbit and timing parameters from
`input_data_files/orbit_and_time_test.dat`, constructs a 200×200 grid in the
orbital plane, and computes dust number density on that grid using all three
DUDI-HC methods: v-integration, delta-ejection, and simple expansion. The
script writes the resulting density matrices and discrepancy maps
(delta-ejection vs v-integration, simple expansion vs delta-ejection) to the
`results/` directory and prints method applicability recommendations based on
a 5% accuracy threshold, mirroring the logic of the original Fortran
`select_method.f90` routine.

### `minimal.py`

`minimal.py` is a compact, end-to-end sanity check and usage example. It builds
a single `Point`, `Source`, and `Comet`, then calls all three main routines
(`simple_expansion`, `delta_ejection`, and `v_integration`) once and prints
the resulting densities to stdout. It is intended as a minimal test that your
Fortran bridge, Python installation, and basic calling convention are working
correctly.

### `example.py`

`example.py` is the Python analogue of `examples/example.f90`. It reads 
ephemeris (`input_data/ephemeridae.dat`), constructs a set of sources on a
sphere around the comet at multiple time steps, and evaluates the dust number
density on a planar grid in the orbital plane. The script uses the
delta-ejection method (as in the Fortran example), performs batched
calculations over sources and points for efficiency, and writes the resulting
2D density matrix to `results/result.dat`, with the layout matching the
Fortran output for further plotting or post-processing.
The runtime of the `example.py` script is ~ 30 seconds with 4 OpenMP threads.
Plot the result with the command:
      `python3 scripts/show_image.py` 


### `Phaethon.py`

`Phaethon.py` is a full Python counterpart of the Fortran `Phaethon` program,
implementing the Phaethon dust environment case study described in the main
README and the original paper. It reads a time series of Phaethon ephemeris
data, builds moving sources and comet states along the trajectory, constructs
a 2D grid around Phaethon in a suitable Cartesian frame, and loops over a
set of grain radii to produce one density map per grain size. Each map is
written as a text matrix to `results/Rg=<value>micron.dat`.

In the Fortran implementation, the dust ejection–rate distribution over
Phaethon’s surface is stored in matrices defined in
`distributions_fun.f90` and used by the`ejection_direction_distribution` function. 
These same matrices are accessed and updated from the main Fortran program 
`phaethon.f90` to orchestrate time-variable dust ejection from Phaethon's surface.
These matrices are also exposed through the Python API and are set directly in 
the body of the `Phaethon.py` script.

For performance reasons, the Python implementation slightly modifies the
algorithm compared to the original Fortran version: instead of updating the
impact-ejecta ratemap at every time step, it groups time intervals into blocks
that share the same ratemap pair and updates the interpolation once per block
while performing batched density evaluations. This batching significantly
reduces runtime but leads to subtle but noticeable differences from the original
Fortran result shown in Fig. 10 of the paper. Also, the resolution of the dust 
number density map around Phaethon is reduced by a factor of two in the Python 
script compared to the Fortran implementation. It corresponds to a fourfold 
decrease in the number of FLOPs. As a result, the runtime of the Phaethon example
in Python is lower (4 minutes with 4 OpenMP threads).
Plot the result with the command:
     `python3 scripts/plot_Fig10.py`
     
## 6. Maintenance notes

### 6.1 Rebuilding the Fortran bridge (release and debug modes)

The shared library `libpy_dudihc_bridge.so` is generated by the build script:

`python_interface/fortran_bridge/build_ctypes_bridge.sh`

By default, the script compiles the Fortran core with optimization flags (-O2)
and OpenMP enabled.

To build the bridge in debug mode, which enables runtime checks, disables
optimizations, and may activate additional debugging output inside Fortran,
set the environment variable `DEBUG=1`:

`DEBUG=1 bash python_interface/fortran_bridge/build_ctypes_bridge.sh`

This produces a slower but more traceable version of the shared library, useful
when diagnosing issues such as NaNs, unexpected densities, mismatched array
shapes, or memory-related bugs.

After rebuilding the library, reinstall the Python package:
`pip install -e .`

### 6.2 Rebuilding or modifying Fortran sources

If **any** Fortran source file is changed, the shared library must be rebuilt.

### 6.3 Continuous Integration

The GitHub Actions workflows located in .github/workflows/ run tests on every
push or pull request.
A full CI run:

  1. Build the Fortran ctypes bridge,

  2. Install the Python package,

  3. Run pytest.

If CI is extended to include macOS or Windows, update the classifiers in
pyproject.toml accordingly and ensure the build script is adapted.
