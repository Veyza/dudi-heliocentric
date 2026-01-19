# python_interface/tests/test_api_bridge.py
import os
import math
import numpy as np
import pytest

from python_interface.dudi_hc.api import v_integration, delta_ejection, simple_expansion
from python_interface.dudi_hc.models import Point, Source, Comet, EjectionSpeedProperties

@pytest.fixture
def pt():
    return Point(
        r=1.0, alpha=0.0, beta=0.0,
        rvector=np.array([1.0, 0.0, 0.0], dtype=float)
    )

@pytest.fixture
def src():
    return Source(
        r=1.0, alphaM=0.0, betaM=0.0,
        rrM=np.array([1.0, 0.0, 0.0], dtype=float),
        zeta=0.0, eta=0.0,
        symmetry_axis=np.array([0.0, 0.0, 1.0], dtype=float),
        ejection_angle_distr=0,  # 0 is allowed
        ud=EjectionSpeedProperties(ud_shape=0, umin=0.0, umax=1.0),
        Nparticles=1.0e10, Tj=0.0, dtau=1.0e-4,
    )

@pytest.fixture
def cm():
    return Comet(
        coords=np.zeros(3, dtype=float),
        Vastvec=np.zeros(3, dtype=float),
        Vast=0.0,
    )

def _assert_finite(x: float):
    assert isinstance(x, float)
    assert math.isfinite(x)

def test_simple_expansion_finite(pt, src):
    y = simple_expansion(pt, src, cloudcentr=[0.0, 0.0, 0.0], dt=0.0)
    _assert_finite(y)

def test_delta_ejection_finite(pt, src, cm):
    y = delta_ejection(pt, src, cm, muR=1.0, dt=0.1, Rast_AU=1.0)
    _assert_finite(y)

def test_v_integration_finite(pt, src, cm):
    y = v_integration(pt, src, cm, muR=1.0, tnow=0.0, Rast_AU=1.0, pericenter=False)
    _assert_finite(y)

def test_diag_returns_zero_and_prints(monkeypatch, pt, src, cm, capfd):
    # Turn on the Fortran-side diagnostic to return 0.0 early
    monkeypatch.setenv("HC_BRIDGE_DIAG", "1")
    y1 = simple_expansion(pt, src, cloudcentr=[0.0, 0.0, 0.0], dt=0.0)
    y2 = delta_ejection(pt, src, cm, muR=1.0, dt=0.1, Rast_AU=1.0)
    y3 = v_integration(pt, src, cm, muR=1.0, tnow=0.0, Rast_AU=1.0, pericenter=False)

    # All three should be exactly 0.0 under DIAG mode
    assert y1 == 0.0
    assert y2 == 0.0
    assert y3 == 0.0

    # Capture low-level FD output produced by Fortran PRINT/printf
    out, err = capfd.readouterr()
    text = (out or "") + (err or "")
    # Sanity: our Fortran DIAG prints should be present
    assert "py_hc_simple_expansion" in text
    assert "py_hc_delta_ejection" in text
    assert "py_hc_v_integration" in text
