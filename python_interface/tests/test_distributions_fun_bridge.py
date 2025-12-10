# python_interface/tests/test_distributions_fun_bridge.py
import numpy as np
import pytest

from python_interface.dudi_hc import api


def test_sizes_and_bounds_roundtrip():
    api.set_lon_limits(-180.0, 180.0)
    lo, hi = api.get_lon_limits()
    assert lo == pytest.approx(-180.0)
    assert hi == pytest.approx(180.0)


def test_rmtmp_roundtrip():
    v = np.array([1.0, -2.0, 3.5], dtype=np.float64)
    api.set_rMtmp(v)
    got = api.get_rMtmp()
    assert got.dtype == np.float64
    assert got.shape == (3,)
    np.testing.assert_allclose(got, v, rtol=0, atol=0)
