import numpy as np
import pytest

from python_interface.dudi_hc.models import (
    Point, Source, Comet, EjectionSpeedProperties, spherical_to_cartesian, normalize
)
from python_interface.dudi_hc import api

def build_ok_models():
    p = Point(1.0, 1.0, 0.5, spherical_to_cartesian(1.0, 1.0, 0.5))
    s = Source(
        r=1.0, alphaM=1.0, betaM=0.0,
        rrM=spherical_to_cartesian(1.0, 1.0, 0.0),
        zeta=0.3, eta=1.2,
        symmetry_axis=normalize(np.array([0.1, 0.2, 0.97], dtype=float)),
        ejection_angle_distr=3,
        ud=EjectionSpeedProperties(ud_shape=1, umin=0.0, umax=0.01),
    )
    Vastvec = np.array([0.0, 1.0, 0.0], dtype=float)
    c = Comet(coords=np.array([1.0, 0.0, 0.0], dtype=float), Vastvec=Vastvec, Vast=float(np.linalg.norm(Vastvec)))
    return p, s, c

def test_delta_ejection_stub_and_validation():
    p, s, c = build_ok_models()
    with pytest.raises(NotImplementedError):
        api.delta_ejection(p, s, c, muR=0.6, dt=0.1, Rast_AU=0.0)

def test_delta_ejection_rejects_negative_dt():
    p, s, c = build_ok_models()
    with pytest.raises(ValueError):
        api.delta_ejection(p, s, c, muR=0.6, dt=-0.1, Rast_AU=0.0)

def test_v_integration_rejects_non_bool_pericenter():
    p, s, c = build_ok_models()
    with pytest.raises(ValueError):
        api.v_integration(p, s, c, muR=0.6, tnow=0.0, Rast_AU=0.0, pericenter=0)  # int, not bool

def test_simple_expansion_requires_cloud_center_vec3():
    p, s, _ = build_ok_models()
    with pytest.raises(ValueError):
        api.simple_expansion(p, s, cloudcentr=np.array([1.0, 0.0], dtype=float), dt=0.1)
