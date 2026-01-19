# python_interface/tests/test_method_to_id.py
import pytest
from python_interface.dudi_hc import api
from python_interface.dudi_hc._bridge_ctypes import (
    METHOD_SIMPLE_EXPANSION,
    METHOD_DELTA_EJECTION,
    METHOD_V_INTEGRATION,
)

def test_method_to_id_accepts_ints_and_strings():
    # ints
    assert api._method_to_id(METHOD_SIMPLE_EXPANSION) == METHOD_SIMPLE_EXPANSION
    assert api._method_to_id(METHOD_DELTA_EJECTION) == METHOD_DELTA_EJECTION
    assert api._method_to_id(METHOD_V_INTEGRATION) == METHOD_V_INTEGRATION

    # canonical names + aliases
    assert api._method_to_id("simple_expansion") == METHOD_SIMPLE_EXPANSION
    assert api._method_to_id("simple") == METHOD_SIMPLE_EXPANSION
    assert api._method_to_id("delta_ejection") == METHOD_DELTA_EJECTION
    assert api._method_to_id("delta") == METHOD_DELTA_EJECTION
    assert api._method_to_id("v_integration") == METHOD_V_INTEGRATION
    assert api._method_to_id("vintegration") == METHOD_V_INTEGRATION
    assert api._method_to_id("v-int") == METHOD_V_INTEGRATION

def test_method_to_id_rejects_unknown():
    with pytest.raises(ValueError):
        api._method_to_id(999)

    with pytest.raises(ValueError):
        api._method_to_id("foobar")