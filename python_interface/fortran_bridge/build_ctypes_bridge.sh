#!/usr/bin/env bash
set -euo pipefail

# Simple, repeatable build that makes one shared lib for ctypes.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC_DIR="${ROOT_DIR}/src"
BRIDGE_DIR="${ROOT_DIR}/python_interface/fortran_bridge"
PKG_DIR="${ROOT_DIR}/python_interface/dudi_hc"
BUILD_DIR="${BRIDGE_DIR}/build_ctypes"

LIB_OUT="${PKG_DIR}/libpy_dudihc_bridge.so"

mkdir -p "${BUILD_DIR}" "${PKG_DIR}"

# Locale to read UTF-8 Fortran files
export LC_ALL=${LC_ALL:-C.UTF-8}
export LANG=${LANG:-C.UTF-8}

# Default: optimized build
FC=${FC:-gfortran}
debug=${DEBUG:-0}

if [[ "$debug" == "1" ]]; then
  echo "• DEBUG build: bounds/undefined checks enabled"
  FFLAGS="-O0 -g -fPIC -fopenmp -fcheck=all -finit-real=snan -finit-local-zero -fbacktrace"
  LDFLAGS="-shared -fopenmp -g"
else
  FFLAGS=${FFLAGS:-"-O2 -fPIC -fopenmp"}
  LDFLAGS=${LDFLAGS:-"-shared -fopenmp"}
fi


# Compile order matters for .mod files
FILES_IN_ORDER=(
  "const.f90"
  "nan_utils.f90"         
  "define_types.f90"
  "help.f90"              
  "distributions_fun.f90"
  "twobody_fun.f90"       
  "data_in.f90"
  "DUDIhc.f90"
  "batching.f90"
)

echo "• Building Fortran core + C-bind bridge → ${LIB_OUT}"
MODFLAGS="-J${BUILD_DIR} -I${BUILD_DIR}"

# 1) Compile core modules to objects
pushd "${SRC_DIR}" >/dev/null
for src in "${FILES_IN_ORDER[@]}"; do
  [[ -f "${src}" ]] || { echo "  (skip) ${src}"; continue; }
  echo "  ${FC} -c ${FFLAGS} ${MODFLAGS} -o '${BUILD_DIR}/${src%.f90}.o' '${src}'"
  ${FC} -c ${FFLAGS} ${MODFLAGS} -o "${BUILD_DIR}/${src%.f90}.o" "${src}"
done
popd >/dev/null

# 2) Compile the bridge (depends on define_types & DUDIhc modules)
echo "  ${FC} -c ${FFLAGS} ${MODFLAGS} -o '${BUILD_DIR}/py_bridge.o' '${BRIDGE_DIR}/py_bridge.f90'"
${FC} -c ${FFLAGS} ${MODFLAGS} -o "${BUILD_DIR}/py_bridge.o" "${BRIDGE_DIR}/py_bridge.f90"

# 3) Link into one shared library
echo "  ${FC} ${LDFLAGS} -o '${LIB_OUT}' ${BUILD_DIR}/*.o"
${FC} ${LDFLAGS} -o "${LIB_OUT}" ${BUILD_DIR}/*.o

echo "• Built: ${LIB_OUT}"

# 4) Quick symbol check (fail early if something went wrong)
if ! nm -D "${LIB_OUT}" | grep -qi " py_hc_v_integration$"; then
  echo "ERROR: py_hc_v_integration not found in ${LIB_OUT}"
  exit 1
fi

# 5) Import smoke test (via ctypes wrapper; you’ll add it next)
python3 - <<'PY'
from python_interface.dudi_hc import _bridge_ctypes as B
print("ctypes bridge loaded OK; functions:", [x for x in dir(B) if x.startswith("call_")])
PY

echo "✓ ctypes bridge build complete."

