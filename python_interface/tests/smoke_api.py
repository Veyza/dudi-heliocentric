import numpy as np
from python_interface.dudi_hc.models import (
    Point, Source, Comet, EjectionSpeedProperties, spherical_to_cartesian, normalize
)
from python_interface.dudi_hc import api

def build_models():
    # Point
    p = Point(r=1.0, alpha=1.0, beta=0.5, rvector=spherical_to_cartesian(1.0, 1.0, 0.5))
    # Source
    s = Source(
        r=1.0, alphaM=1.0, betaM=0.0,
        rrM=spherical_to_cartesian(1.0, 1.0, 0.0),
        zeta=0.3, eta=1.2,
        symmetry_axis=normalize(np.array([0.1, 0.2, 0.97], dtype=float)),
        ejection_angle_distr=3,
        ud=EjectionSpeedProperties(ud_shape=1, umin=0.0, umax=0.01),
    )
    # Comet
    Vastvec = np.array([0.0, 1.0, 0.0], dtype=float)
    c = Comet(coords=np.array([1.0, 0.0, 0.0], dtype=float), Vastvec=Vastvec, Vast=float(np.linalg.norm(Vastvec)))
    return p, s, c

def main():
    p, s, c = build_models()

    def expect_ni(callable_desc, fn, *args, **kwargs):
        try:
            fn(*args, **kwargs)
        except NotImplementedError as e:
            print(f"{callable_desc}: OK (stub) ->", e)
        else:
            raise AssertionError(f"{callable_desc}: expected NotImplementedError")

    expect_ni("v_integration", api.v_integration, p, s, c, muR=0.6, tnow=0.0, Rast_AU=0.0, pericenter=False)
    expect_ni("delta_ejection", api.delta_ejection, p, s, c, muR=0.6, dt=0.1, Rast_AU=0.0)
    expect_ni("simple_expansion", api.simple_expansion, p, s, cloudcentr=np.array([1.0, 0.0, 0.0], dtype=float), dt=0.05)

    print("All API stubs validated.")

if __name__ == "__main__":
    main()
