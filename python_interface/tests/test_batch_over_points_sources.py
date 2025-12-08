# python_interface/tests/test_batch_over_points_sources.py
import numpy as np
import pytest

from python_interface.dudi_hc import api
from python_interface.dudi_hc.models import (
    Point, Source, Comet, EjectionSpeedProperties,
    spherical_to_cartesian, normalize,
)

def _build_one_source(Tj: float = 0.0) -> Source:
    return Source(
        r=1.0, alphaM=1.0, betaM=0.0,
        rrM=spherical_to_cartesian(1.0, 1.0, 0.0),
        zeta=0.3, eta=1.2,
        symmetry_axis=normalize(np.array([0.1, 0.2, 0.97], dtype=float)),
        ejection_angle_distr=3,
        ud=EjectionSpeedProperties(ud_shape=1, umin=0.0, umax=0.01),
        Nparticles=1e10,
        Tj=Tj,
        dtau=0.0004,
    )

def _build_comet() -> Comet:
    Vastvec = np.array([0.0, 1.0, 0.0], dtype=float)
    return Comet(
        coords=np.array([1.0, 0.0, 0.0], dtype=float),
        Vastvec=Vastvec,
        Vast=float(np.linalg.norm(Vastvec)),
    )

def test_batch_over_points_sources_empty_points_returns_empty():
    # Nt = 2, Ns = 1 (any small, valid configuration)
    sources_by_time = [[_build_one_source(Tj=0.0)], [_build_one_source(Tj=0.1)]]
    comets_by_time = [_build_comet(), _build_comet()]

    densities = api.batch_over_points_sources(
        points=[],  # <- important bit
        sources_by_time=sources_by_time,
        comets_by_time=comets_by_time,
        muR=0.6,
        tnow=0.2,
        Rast_AU=0.0,
        pericenter=False,
        method="delta_ejection",
    )

    assert isinstance(densities, np.ndarray)
    assert densities.shape == (0,)

def test_batch_over_points_sources_rejects_mismatched_comets_length():
    p = Point(
        r=1.0, alpha=1.0, beta=0.5,
        rvector=spherical_to_cartesian(1.0, 1.0, 0.5),
    )
    sources_by_time = [[_build_one_source(Tj=0.0)], [_build_one_source(Tj=0.1)]]  # Nt=2
    comets_by_time = [_build_comet()]  # Nt=1 -> mismatch

    with pytest.raises(ValueError):
        api.batch_over_points_sources(
            points=[p],
            sources_by_time=sources_by_time,
            comets_by_time=comets_by_time,
            muR=0.6,
            tnow=0.2,
            Rast_AU=0.0,
            pericenter=False,
            method="delta_ejection",
        )

 

## 4. `batch_over_points_sources` flat `sources_by_time` (1D) handling

# We consider a special case that wraps a flat list of `Source` into shape `(Nt, 1)` if `sources_by_time[0]` is a `Source`.

def test_batch_over_points_sources_accepts_flat_source_list():
    p = Point(
        r=1.0, alpha=1.0, beta=0.5,
        rvector=spherical_to_cartesian(1.0, 1.0, 0.5),
    )
    # Flat list: Nt=3, Ns=1 after auto-wrapping
    sources_flat = [
        _build_one_source(Tj=0.0),
        _build_one_source(Tj=0.1),
        _build_one_source(Tj=0.2),
    ]
    comets_by_time = [_build_comet(), _build_comet(), _build_comet()]

    densities = api.batch_over_points_sources(
        points=[p],
        sources_by_time=sources_flat,  # <- flat list
        comets_by_time=comets_by_time,
        muR=0.6,
        tnow=0.3,
        Rast_AU=0.0,
        pericenter=False,
        method="delta_ejection",
    )

    assert densities.shape == (1,)
    assert np.isfinite(densities[0])
