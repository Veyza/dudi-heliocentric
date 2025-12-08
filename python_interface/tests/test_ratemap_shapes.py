# python_interface/tests/test_ratemap_shapes.py
import numpy as np
import pytest

from python_interface.dudi_hc import api

def test_set_lats_rejects_wrong_length():
    nlats, _ = api.get_ratemap_dims()
    good = np.linspace(-89.0, 89.0, nlats, dtype=float)
    api.set_lats(good)  # should work

    bad = np.linspace(-89.0, 89.0, nlats + 1, dtype=float)
    with pytest.raises(ValueError):
        api.set_lats(bad)

def test_set_lons_rejects_wrong_length():
    _, nlons = api.get_ratemap_dims()
    good = np.linspace(-179.0, 179.0, nlons, dtype=float)
    api.set_lons(good)  # should work

    bad = np.linspace(-179.0, 179.0, nlons - 1, dtype=float)
    with pytest.raises(ValueError):
        api.set_lons(bad)

def test_set_ratemap_rejects_swapped_shape():
    nlats, nlons = api.get_ratemap_dims()

    # correct orientation (nlons, nlats) should work
    good = np.zeros((nlons, nlats), dtype=np.float64)
    api.set_ratemap(good)

    # swapped dims (nlats, nlons) must raise
    bad = np.zeros((nlats, nlons), dtype=np.float64)
    with pytest.raises(ValueError):
        api.set_ratemap(bad)
