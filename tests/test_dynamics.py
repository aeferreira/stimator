import pytest
from stimator import read_model
from stimator import dynamics as dyn

demomodel = """
title a simple 2 step system
v1: A -> B, rate = V / (Km1 + A), V = 1, Km = 1

v2: B ->  , rate = V * c2 * B**3

V  = sqrt(4.0)
Km1 = 1
c2 = 0.2

find c2 in [0, 1.2]

init: B = 0.4, A = 1

-> vin = 2 * A * v1.Km
~ t1 = A + B + vin
~ t2 = v1.V * A * step(t, 1.0)
# ~ t3 = v1.V * A * max(t, 1.0)"""


@pytest.fixture
def m():
    return read_model(demomodel)

def test_genStoichiometryMatrix(m):
    N = dyn.genStoichiometryMatrix(m)
    nreactions = len(m.reactions) 
    nvars = len(m.varnames)
    assert N.shape == (nvars, nreactions)

