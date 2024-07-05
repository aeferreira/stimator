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


def test_genStoichiometryMatrix(m: st.Model):
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


def test_state2array(m: st.Model):
    v = dyn.init2array(m)
    nvars = len(m.varnames)
    assert isinstance(v, np.ndarray)
    assert v.shape == (nvars,)


def test_identifiersInExpr():
    expr = "v1.V / (Km1 + A) + V * c2 * B**3"
    allids = dyn.identifiersInExpr(expr)
    assert len(allids) == 7
    for name in ["v1", "V", "Km1", "A", "V", "c2", "B"]:
        assert name in allids


def test_rate_strings(m: st.Model):
    rs = dyn.rates_strings(m, fully_qualified=False)
    assert rs["v1"] == "V / (Km1 + A)"
    assert rs["v2"] == "V * c2 * B**3"
    rs = dyn.rates_strings(m, fully_qualified=True)
    assert rs["v1"] == "v1.V / (Km1 + A)"
    assert rs["v2"] == "V * c2 * B**3"


def test_dXdt_strings(m: st.Model):
    dxdt_strs = dyn.dXdt_strings(m)
    assert dxdt_strs["A"] == "-v1.V/(A + Km1)"
    assert dxdt_strs["B"] == "-B**3*V*c2 + v1.V/(A + Km1)"


def test_gen_canonical_symbmap(m: st.Model):
    symbmap = dyn._gen_canonical_symbmap(m)["s_table"]
    assert len(symbmap) == 7
    for name in "A B V Km1 c2 v1.Km v1.V".split():
        assert name in symbmap
        assert symbmap[name].startswith("_symbol_Id")


def test_string_differentiation(m: st.Model):
    symbols = dyn._gen_canonical_symbmap(m)
    dif = dyn._differentiate_expr
    dxdt_strs = dyn.dXdt_strings(m)
    expA = dxdt_strs["A"]
    assert expA == "-v1.V/(A + Km1)"
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
    assert dif(expA, "A", symbols) == "v1.V/(A + Km1)**2"
    assert dif(expA, "B", symbols) == "0.0"
    assert dif(expA, "v1.V", symbols) == "-1/(A + Km1)"
    assert dif(expA, "c2", symbols) == "0.0"


def test_Jacobian_strings(m: st.Model):
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


def test_dfdp_strings(m: st.Model):
    parnames = "c2 v1.V".split()
    dfdp_strs = dyn.dfdp_strings(m, parnames)
    nvars = len(m.varnames)
    assert len(dfdp_strs) == nvars
    assert len(dfdp_strs[0]) == 2
    # (d dA/dt / d c2) = 0.0
    # (d dA/dt / d v1.V) = -1/(A + Km1)
    # (d dB/dt / d c2) = -B**3*V
    # (d dB/dt / d v1.V) = 1/(A + Km1)
    assert dfdp_strs[0][0] == "0.0"
    assert dfdp_strs[0][1] == "-1/(A + Km1)"
    assert dfdp_strs[1][0] == "-B**3*V"
    assert dfdp_strs[1][1] == "1/(A + Km1)"

def test_dfdp_strings_with_unknown(m: st.Model):
    parnames = "c3 v1.V".split()
    dfdp_strs = dyn.dfdp_strings(m, parnames)
    nvars = len(m.varnames)
    assert len(dfdp_strs) == nvars
    assert len(dfdp_strs[0]) == 2
    # (d dA/dt / d c3) = 0.0
    # (d dA/dt / d v1.V) = -1/(A + Km1)
    # (d dB/dt / d c3) = 0.0
    # (d dB/dt / d v1.V) = 1/(A + Km1)
    assert dfdp_strs[0][0] == "0.0"
    assert dfdp_strs[0][1] == "-1/(A + Km1)"
    assert dfdp_strs[1][0] == "0.0"
    assert dfdp_strs[1][1] == "1/(A + Km1)"


def test_gen_calc_symbmap(m: st.Model):
    symbmap = dyn._gen_calc_symbmap(m)
    assert len(symbmap) == 8
    # A        --> variables[0]
    # B        --> variables[1]
    # vin      --> input_variables[0]
    # V        --> 2
    # Km1      --> 1
    # c2       --> 0.2
    # v1.Km    --> 1
    # v1.V     --> 1
    for name in "A B".split():
        assert name in symbmap
        assert symbmap[name].startswith("variables[")
    for name in "vin".split():
        assert name in symbmap
        assert symbmap[name].startswith("input_variables[")
    for name in "V Km1 c2 v1.Km v1.V".split():
        assert name in symbmap
        assert isinstance(float(symbmap[name]), float)

def test_gen_calc_symbmap_with_uncertain(m: st.Model):
    symbmap = dyn._gen_calc_symbmap(m, with_uncertain=True)
    assert len(symbmap) == 8
    # A        --> variables[0]
    # B        --> variables[1]
    # vin      --> input_variables[0]
    # c2       --> m_Parameters[0]
    # V        --> 2
    # Km1      --> 1
    # v1.Km    --> 1
    # v1.V     --> 1
    for name in "A B".split():
        assert symbmap[name].startswith("variables[")
    for name in "vin".split():
        assert symbmap[name].startswith("input_variables[")
    for name in "c2".split():
        assert symbmap[name].startswith("m_Parameters[")
    for name in "V Km1 v1.Km v1.V".split():
        assert isinstance(float(symbmap[name]), float)

# ********** Testing calc_string **************************
# calcstring for v1 = v1.V / (Km1 + A)
#     1 / (1 + variables[0])
# calcstring for v2 = V * c2 * B**3
#     2 * 0.2 * variables[1]**3
# calcstring for t1 = A + B + vin
#     variables[0] + variables[1] + input_variables[0]
# calcstring for t2 = v1.V * A * step(t, 1.0)
#     1 * variables[0] * step(t, 1.0)
# calcstring for vin = 2 * A * v1.Km
#     2 * variables[0] * 1
# calcstring for v2 with uncertain parameters:
#          2 * m_Parameters[0] * variables[1]**3

def test_calc_string(m: st.Model):
    symbmap = dyn._gen_calc_symbmap(m, with_uncertain=False)
    calc_strs = {}
    for v in (
        m.reactions.v1,
        m.reactions.v2,
        m.transformations.t1,
        m.transformations.t2,
        m.input_variables.vin,
    ):
        vstr = v(fully_qualified=True)
        calc_strs[v.name] = dyn.calc_string(vstr, symbmap)
