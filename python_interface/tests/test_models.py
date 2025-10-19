import numpy as np
import pytest

from python_interface.dudi_hc.models import (
    Point, Source, Comet, EjectionSpeedProperties,
    spherical_to_cartesian, normalize, as_vec3
)

def test_point_ok():
    r, alpha, beta = 1.0, 1.0, 0.5
    p = Point(r=r, alpha=alpha, beta=beta, rvector=spherical_to_cartesian(r, alpha, beta))
    assert p.r == pytest.approx(1.0)
    assert p.rvector.shape == (3,)

def test_point_rejects_bad_rvector_shape():
    with pytest.raises(ValueError):
        Point(r=1.0, alpha=1.0, beta=0.5, rvector=np.array([1.0, 2.0], dtype=float))

def test_source_ok():
    rrM = spherical_to_cartesian(1.0, 1.0, 0.0)
    axis = normalize(np.array([0.1, 0.2, 0.97], dtype=float))
    ud = EjectionSpeedProperties(ud_shape=1, umin=0.0, umax=0.01)
    s = Source(
        r=1.0, alphaM=1.0, betaM=0.0, rrM=rrM,
        zeta=0.3, eta=1.2,
        symmetry_axis=axis,
        ejection_angle_distr=3,
        ud=ud
    )
    assert s.rrM.shape == (3,)
    assert np.isclose(np.linalg.norm(s.symmetry_axis), 1.0)

def test_source_rejects_non_unit_axis():
    rrM = spherical_to_cartesian(1.0, 1.0, 0.0)
    bad_axis = np.array([10.0, 0.0, 0.0], dtype=float)  # norm != 1
    ud = EjectionSpeedProperties(ud_shape=1, umin=0.0, umax=0.01)
    with pytest.raises(ValueError):
        Source(
            r=1.0, alphaM=1.0, betaM=0.0, rrM=rrM,
            zeta=0.3, eta=1.2,
            symmetry_axis=bad_axis,
            ejection_angle_distr=3,
            ud=ud
        )

def test_ejection_speed_properties_bounds():
    with pytest.raises(ValueError):
        EjectionSpeedProperties(ud_shape=1, umin=-1.0, umax=0.0)
    with pytest.raises(ValueError):
        EjectionSpeedProperties(ud_shape=1, umin=0.01, umax=0.0)

def test_comet_ok():
    coords = np.array([1.0, 0.0, 0.0], dtype=float)
    Vastvec = np.array([0.0, 1.0, 0.0], dtype=float)
    c = Comet(coords=coords, Vastvec=Vastvec, Vast=float(np.linalg.norm(Vastvec)))
    assert c.coords.shape == (3,)
    assert c.Vast >= 0.0

def test_comet_rejects_bad_vectors():
    with pytest.raises(ValueError):
        Comet(coords=np.array([1.0, 0.0], dtype=float), Vastvec=np.array([0.0, 1.0, 0.0], dtype=float), Vast=1.0)