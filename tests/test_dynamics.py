import pytest
import numpy as np
import stimator as st
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
    return st.read_model(demomodel)


def test_genStoichiometryMatrix(m):
    N = dyn.genStoichiometryMatrix(m)
    nreactions = len(m.reactions) 
    nvars = len(m.varnames)
    # Stoichiometry matrix:
    #    v1  v2
    # A [-1.  0.]
    # B [ 1. -1.]
    assert isinstance(N, np.ndarray)
    assert N.shape == (nvars, nreactions)
    assert N[0, 0] == -1
    assert N[1, 1] == -1
    assert N[0, 1] == 0
    assert N[1, 0] == 1


def test_state2array(m):
    v = dyn.init2array(m)
    nvars = len(m.varnames)
    assert isinstance(v, np.ndarray)
    assert v.shape == (nvars, )


def test_identifiersInExpr():
    expr = 'v1.V / (Km1 + A) + V * c2 * B**3'
    allids = dyn.identifiersInExpr(expr)
    assert len(allids) == 7
    for name in ['v1', 'V', 'Km1', 'A', 'V', 'c2', 'B']:
        assert name in allids


def test_rate_strings(m):
    rs = dyn.rates_strings(m, fully_qualified=False)
    assert rs['v1'] == 'V / (Km1 + A)'
    assert rs['v2'] == 'V * c2 * B**3'
    rs = dyn.rates_strings(m, fully_qualified=True)
    assert rs['v1'] == 'v1.V / (Km1 + A)'
    assert rs['v2'] == 'V * c2 * B**3'


def test_dXdt_strings(m):
    dxdt_strs = dyn.dXdt_strings(m)
    assert dxdt_strs['A'] == '-v1.V/(A + Km1)'
    assert dxdt_strs['B'] == '-B**3*V*c2 + v1.V/(A + Km1)'


def test_gen_canonical_symbmap(m):
    symbmap = dyn._gen_canonical_symbmap(m)['s_table']
    assert len(symbmap) == 7
    for name in 'A B V Km1 c2 v1.Km v1.V'.split():
        assert name in symbmap
        assert symbmap[name].startswith('_symbol_Id')

def test_string_differentiation(m):
    symbols = dyn._gen_canonical_symbmap(m)
    dif = dyn._differentiate_expr
    dxdt_strs = dyn.dXdt_strings(m)
    expA = dxdt_strs['A']
    assert expA == '-v1.V/(A + Km1)'
    # variable A 
    # expression = -v1.V/(A + Km1)
    # d / d A = v1.V/(A + Km1)**2
    # d / d B = 0.0
    # ---
    # d / d V = 0.0
    # d / d Km1 = v1.V/(A + Km1)**2
    # d / d c2 = 0.0
    # d / d v1.Km = 0.0
    # d / d v1.V = -1/(A + Km1)
    assert dif(expA, 'A', symbols) == "v1.V/(A + Km1)**2"
    assert dif(expA, 'B', symbols) == "0.0"
    assert dif(expA, 'v1.V', symbols) == "-1/(A + Km1)"
    assert dif(expA, 'c2', symbols) == "0.0"

def test_Jacobian_strings(m):
    nvars = len(m.varnames)
    j_strings = dyn.Jacobian_strings(m)
    # assert j_strings.shape == (len(vnames), len(vnames))
    assert len(j_strings) == nvars
    assert len(j_strings[0]) == nvars
    # (d dA/dt / d A) = v1.V/(A + Km1)**2
    # (d dA/dt / d B) = 0.0
    # (d dB/dt / d A) = -v1.V/(A + Km1)**2
    # (d dB/dt / d B) = -3*B**2*V*c2
    assert j_strings[0][0] == "v1.V/(A + Km1)**2"
    assert j_strings[0][1] == "0.0"
    assert j_strings[1][0] == "-v1.V/(A + Km1)**2"
    assert j_strings[1][1] == "-3*B**2*V*c2"