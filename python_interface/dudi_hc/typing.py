"""
Typing aliases for the Python interface.

- FloatArray: np.ndarray of float64 (any shape)
- Vec3:       np.ndarray of float64 with expected shape (3,)
Shape is enforced at runtime where relevant (e.g., via as_vec3).
"""
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
Vec3 = NDArray[np.float64]
