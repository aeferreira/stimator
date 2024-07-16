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


class dummy:
    pass


v1 = dummy()
v1.V = 1.01

TVs_FULLY = {
    s: 1.01
    for s in ("V", "Km1", "A", "V", "c2", "B", "d_A_d_init_B", "d_B_d_init_B")
}
TVs_FULLY["v1"] = v1


def same_expr(expr1, expr2):
    v1 = eval(expr1, TVs_FULLY)
    v2 = eval(expr2, TVs_FULLY)
    return v1 == pytest.approx(v2)


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
    assert same_expr(rs["v1"], "V / (Km1 + A)")
    assert same_expr(rs["v2"], "V * c2 * B**3")
    rs = dyn.rates_strings(m, fully_qualified=True)
    assert same_expr(rs["v1"], "v1.V / (Km1 + A)")
    assert same_expr(rs["v2"], "V * c2 * B**3")


def test_dXdt_strings(m: st.Model):
    dxdt_strs = dyn.dXdt_strings(m)
    assert same_expr(dxdt_strs["A"], "-v1.V/(A + Km1)")
    assert same_expr(dxdt_strs["B"], "-B**3*V*c2 + v1.V/(A + Km1)")


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
    assert same_expr(expA, "-v1.V/(A + Km1)")
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
    assert same_expr(dif(expA, "A", symbols), "v1.V/(A + Km1)**2")
    assert same_expr(dif(expA, "B", symbols), "0.0")
    assert same_expr(dif(expA, "v1.V", symbols), "-1/(A + Km1)")
    assert same_expr(dif(expA, "c2", symbols), "0.0")


def test_Jacobian_strings(m: st.Model):
    nvars = len(m.varnames)
    j_strings = dyn.Jacobian_strings(m)
    assert len(j_strings) == nvars
    assert len(j_strings[0]) == nvars
    # (d dA/dt / d A) = v1.V/(A + Km1)**2
    # (d dA/dt / d B) = 0.0
    # (d dB/dt / d A) = -v1.V/(A + Km1)**2
    # (d dB/dt / d B) = -3*B**2*V*c2
    assert same_expr(j_strings[0][0], "v1.V/(A + Km1)**2")
    assert same_expr(j_strings[0][1], "0.0")
    assert same_expr(j_strings[1][0], "-v1.V/(A + Km1)**2")
    assert same_expr(j_strings[1][1], "-3*B**2*V*c2")


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
    assert same_expr(dfdp_strs[0][0], "0.0")
    assert same_expr(dfdp_strs[0][1], "-1/(A + Km1)")
    assert same_expr(dfdp_strs[1][0], "-B**3*V")
    assert same_expr(dfdp_strs[1][1], "1/(A + Km1)")


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
    assert same_expr(dfdp_strs[0][0], "0.0")
    assert same_expr(dfdp_strs[0][1], "-1/(A + Km1)")
    assert same_expr(dfdp_strs[1][0], "0.0")
    assert same_expr(dfdp_strs[1][1], "1/(A + Km1)")


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
    # v1 = v1.V / (Km1 + A)
    # v2 = V * c2 * B**3
    # calcstring for t1 = A + B + vin
    # calcstring for t2 = v1.V * A * step(t, 1.0)
    # calcstring for vin = 2 * A * v1.Km
    assert calc_strs["v1"] == "1 / (1 + variables[0])"
    assert calc_strs["v2"] == "2 * 0.2 * variables[1]**3"
    assert (
        calc_strs["t1"] == "variables[0] + variables[1] + input_variables[0]"
    )
    assert calc_strs["t2"] == "1 * variables[0] * step(t, 1.0)"
    assert calc_strs["vin"] == "2 * variables[0] * 1"


def test_calc_string_uncertain(m: st.Model):
    symbmap = dyn._gen_calc_symbmap(m, with_uncertain=True)
    unc_v2 = dyn.calc_string(m.reactions.v2(fully_qualified=True), symbmap)
    assert unc_v2 == "2 * m_Parameters[0] * variables[1]**3"


def test_all_rates_func(m: st.Model):
    func = dyn.all_rates_func(m)
    # Operating point
    varvalues = 1.0, 0.4
    # at t == 0
    t = 0.0
    ivs, vs, ts = func(varvalues, t)
    assert vs == pytest.approx((0.5, 0.0256))
    assert ts == pytest.approx((3.4, 0.0))
    assert ivs == pytest.approx((2.0,))
    # at t == 2.0
    t = 2.0
    ivs, vs, ts = func(varvalues, t)
    assert vs == pytest.approx((0.5, 0.0256))
    assert ts == pytest.approx((3.4, 1.0))
    assert ivs == pytest.approx((2.0,))


def test_add_dSdt_to_model(m: st.Model):
    vnames = m.varnames
    dxdtstrs = dyn.dXdt_strings(m)
    # before adding sensitivities
    assert len(vnames) == 2
    assert m.get_init("A") == 1.0
    assert m.get_init("B") == 0.4
    assert same_expr(dxdtstrs["A"], "-v1.V/(A + Km1)")
    assert same_expr(dxdtstrs["B"], "-B**3*V*c2 + v1.V/(A + Km1)")
    # after adding sensitivities
    pars = "Km2 v1.V init.B".split()
    Snames = dyn.add_dSdt_to_model(m, pars)
    assert len(Snames) == 2 * 3
    assert ("A", "v1.V", "d_A_d_v1_V") in Snames
    vnames = m.varnames
    dxdtstrs = dyn.dXdt_strings(m)
    assert len(vnames) == 2 + 2 * 3
    assert m.get_init("A") == 1.0
    assert m.get_init("B") == 0.4
    assert m.get_init("d_A_d_init_B") == 0.0
    assert m.get_init("d_B_d_init_B") == 1.0
    assert same_expr(dxdtstrs["A"], "-v1.V/(A + Km1)")
    assert same_expr(dxdtstrs["B"], "-B**3*V*c2 + v1.V/(A + Km1)")
    assert same_expr(
        dxdtstrs["d_A_d_init_B"], "v1.V*d_A_d_init_B/(A + Km1)**2"
    )
    assert same_expr(
        dxdtstrs["d_B_d_init_B"],
        "-3*B**2*c2*d_B_d_init_B*V - v1.V*d_A_d_init_B/(A + Km1)**2",
    )
