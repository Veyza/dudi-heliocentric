# python_interface/tests/test_distributions_fun_bridge.py
import numpy as np
import pytest

from python_interface.dudi_hc import api


def test_sizes_and_bounds_roundtrip():
    nx = api.get_nlons()
    ny = api.get_nlats()
    assert isinstance(nx, int) and nx > 0
    assert isinstance(ny, int) and ny > 0

    api.set_lon_bounds(-180.0, 180.0)
    lo, hi = api.get_lon_bounds()
    assert lo == pytest.approx(-180.0)
    assert hi == pytest.approx(180.0)


def test_lats_lons_roundtrip_and_dtype():
    nx = api.get_nlons()
    ny = api.get_nlats()

    lats = np.linspace(-89.0, 89.0, ny, dtype=np.float32)
    lons = np.linspace(-179.0, 179.0, nx, dtype=np.float32)

    # Set as float32
    api.set_lats(lats)
    api.set_lons(lons)

    # Get back; implementation returns float32 arrays
    lats_rt = api.get_lats()
    lons_rt = api.get_lons()
    assert lats_rt.dtype == np.float32
    assert lons_rt.dtype == np.float32
    assert lats_rt.shape == (ny,)
    assert lons_rt.shape == (nx,)
    np.testing.assert_allclose(lats_rt, lats, rtol=0, atol=0)

    # Also accept float64 on input (wrapper converts)
    api.set_lats(lats.astype(np.float64))
    api.set_lons(lons.astype(np.float64))
    np.testing.assert_allclose(api.get_lats(), lats, rtol=0, atol=0)
    np.testing.assert_allclose(api.get_lons(), lons, rtol=0, atol=0)


def test_rmtmp_roundtrip():
    v = np.array([1.0, -2.0, 3.5], dtype=np.float64)
    api.set_rmtmp(v)
    got = api.get_rmtmp()
    assert got.dtype == np.float64
    assert got.shape == (3,)
    np.testing.assert_allclose(got, v, rtol=0, atol=0)


def test_set_maps_shape_validation():
    nx = api.get_nlons()
    ny = api.get_nlats()

    # Correct shape works (should not raise)
    ok = np.asfortranarray(np.zeros((nx, ny), dtype=np.float64))
    api.set_rmap1(ok)
    api.set_rmap2(ok)
    api.set_ratemap(ok)

    # Wrong shapes should raise ValueError from the wrapper
    bad1 = np.zeros((ny, nx), dtype=np.float64)     # swapped
    bad2 = np.zeros((nx, ny + 1), dtype=np.float64) # off by one
    with pytest.raises(ValueError):
        api.set_rmap1(bad1)
    with pytest.raises(ValueError):
        api.set_rmap2(bad2)
