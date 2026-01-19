import numpy as np
from python_interface.dudi_hc import api

# sizes
nx, ny = api.get_nlons(), api.get_nlats()

# lon bounds
api.set_lon_bounds(-180.0, 180.0)
print("bounds:", api.get_lon_bounds())

# lats/lons
api.set_lats(np.linspace(-89.0, 89.0, ny, dtype=np.float32))
api.set_lons(np.linspace(-179.0, 179.0, nx, dtype=np.float32))
print("lats shape:", api.get_lats().shape, api.get_lats().dtype)
print("lons shape:", api.get_lons().shape, api.get_lons().dtype)

# maps (float64, (nlons, nlats), Fortran-order recommended)
A = np.asfortranarray(np.zeros((nx, ny), dtype=np.float64))
api.set_rmap1(A)
api.set_rmap2(A + 1.0)
api.set_ratemap(A + 2.0)

# rMtmp
api.set_rmtmp([1.0, 0.0, 0.0])
print("rMtmp:", api.get_rmtmp())
