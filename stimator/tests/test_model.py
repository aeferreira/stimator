from stimator import *
from stimator.model import isPairOfNums
from nose.tools import *

def test_M1():
    """test Model __init__ empty ()"""
    m = Model()
    assert isinstance(m, Model)

def test_M2():
    """test Model __init__ with title"""
    m = Model("My first model")
    assert isinstance(m, Model)
    assert m['title'] == "My first model"

def test_react1():
    """test react(string, int or float)"""
    m = Model("My first model")
    m.v1 = react("A->B", 4)
    m.v2 = react("B->C", 2.0)
    assert isinstance(m.v1, model.Reaction)
    assert isinstance(m.v2, model.Reaction)
    assert get_name(m.v1) == 'v1'
    assert get_name(m.v2) == 'v2'
    assert m.v1()== str(float(4))+ "*A"
    assert m.v2()== str(float(2.0))+"*B"
    check, msg = m.checkRates()
    assert check 

def test_react2():
    """test react(string, string)"""
    m = Model("My first model")
    m.v1 = react("A->B", " 4*A/(p1+A)-B ")
    m.p1 = 2
    assert get_name(m.v1) == 'v1'
    assert isinstance(m.v1, model.Reaction)
    assert m.v1()== "4*A/(p1+A)-B"
    check, msg = m.checkRates()
    assert check 

def test_react2b():
    """test react(string, string) with math functions"""
    m = Model("My first model")
    m.v1 = react("A->B", " 4*sqrt(A)/(p1+sin(A))-B ")
    m.p1 = 2
    assert get_name(m.v1) == 'v1'
    assert isinstance(m.v1, model.Reaction)
    assert m.v1()== "4*sqrt(A)/(p1+sin(A))-B"
    check, msg = m.checkRates()
    assert check 

def test_react2c():
    """test react(string, string) with kinetics functions"""
    m = Model("My first model")
    m.v1 = react("A->B", " 4*A*step(t,1.0)")
    m.p1 = 2
    assert get_name(m.v1) == 'v1'
    assert isinstance(m.v1, model.Reaction)
    assert m.v1()== "4*A*step(t,1.0)"
    check, msg = m.checkRates()
    assert check 

@raises(model.BadStoichError)
def test_react3():
    """test Bad stoichiometry"""
    m = Model("My first model")
    m.v1 = react("A->##B", " 4*A/(p1+A)-B ")

def test_react4():
    """test Bad rate law (unknown ID)"""
    m = Model("My first model")
    m.v1 = react("A->B", " 4*A/(p2+A)-B ")
    m.p1 = 2
    assert get_name(m.v1) == 'v1'
    assert isinstance(m.v1, model.Reaction)
    assert m.v1()== "4*A/(p2+A)-B"
    check, msg = m.checkRates()
    assert not check 

def test_react5():
    """test Bad rate law (malformed expression)"""
    m = Model("My first model")
    m.v1 = react("A->B", " 4*A/(p1+A-B ")
    m.p1 = 2
    assert get_name(m.v1) == 'v1'
    assert isinstance(m.v1, model.Reaction)
    assert m.v1()== "4*A/(p1+A-B"
    check, msg = m.checkRates()
    assert not check 

def test_react6():
    """test Bad rate law (fp overflow)"""
    m = Model("My first model")
    m.v1 = react("A->B", " 1e100**10000 * 4*A/(p1+A)-B ")
    m.p1 = 2
    assert get_name(m.v1) == 'v1'
    assert isinstance(m.v1, model.Reaction)
    assert m.v1()== "1e100**10000 * 4*A/(p1+A)-B"
    check, msg = m.checkRates()
    assert not check 

def test_par1():
    """test assignment of parameters"""
    m = Model("My first model")
    m.p1 = 4
    m.p2 = 3.0
    assert isinstance(m.p1, model.ConstValue)
    assert get_name(m.p1) == "p1"
    assert isinstance(m.p2, model.ConstValue)
    assert get_name(m.p2) == "p2"
    assert m.p1 == 4.0
    assert m.p2 == 3.0

def test_par_in_rates1():
    """test assignment of parameters 'local' to reactions"""
    m = Model("My first model")
    m.v1 = react("A->B", " p2*A/(p1+A)-B ", pars={'p1':4})
    m.p2 = 3.0
    assert get_name(m.v1) == 'v1'
    assert isinstance(m.v1, model.Reaction)
    assert m.v1()== "p2*A/(p1+A)-B"
    check, msg = m.checkRates()
    assert check 
    assert isinstance(m.v1.p1, model.ConstValue)
    assert get_name(m.v1.p1) == "p1"
    assert isinstance(m.p2, model.ConstValue)
    assert get_name(m.p2) == "p2"
    assert m.v1.p1 == 4.0
    assert m.p2 == 3.0

def test_par_from_rates1():
    """test rate expressions with parameters 'local' to reactions"""
    m = Model("My first model")
    m.v1 = react("A->B", " p2*A/(p1+A)-B ", pars={'p1':4})
    m.v2 = react("B->C", "2*v1.p1*B")
    m.p2 = 3.0
    assert get_name(m.v2) == 'v2'
    assert isinstance(m.v1, model.Reaction)
    assert m.v2()== "2*v1.p1*B"
    check, msg = m.checkRates()
    assert check 
    assert isinstance(m.v1.p1, model.ConstValue)
    assert get_name(m.v1.p1) == "p1"
    assert isinstance(m.p2, model.ConstValue)
    assert get_name(m.p2) == "p2"
    assert m.v1.p1 == 4.0
    assert m.p2 == 3.0

def test_isPairOfNums():
    """test isPairOfNums function"""
    p1 = 1,10.0
    p2 = (1, 10.0)
    p3 = [1,10.0]
    p4 = [2.3, 4.4]
    p5 = 3
    p6 = (4)
    p7 = ("2.3", 4.5)
    p8 = [8.8]
    p9 = (1,2.3,4)
    assert isPairOfNums(p1)
    assert isPairOfNums(p2)
    assert isPairOfNums(p3)
    assert isPairOfNums(p4)
    assert not isPairOfNums(p5)
    assert not isPairOfNums(p6)
    assert not isPairOfNums(p7)
    assert not isPairOfNums(p8)
    assert not isPairOfNums(p9)
    
def test_par2():
    """test assignment of parameters with bounds"""
    m = Model("My first model")
    m.p1 = 4
    m.p2 = 3.0
    m.p1 = 1,10 #tuple or list
    m.p2 = [1, 9.5]
    m.p3 = 5
    m.p4 = 6
    m.p4.uncertainty(1, 8.5) # or uncertainty function
    m.p5 = 0,10 # bounds create parameter with midpoint
    assert m.p1 == 4.0
    assert m.p2 == 3.0
    assert m.p3 == 5.0
    assert m.p4 == 6.0
    assert m.p3.bounds is None
    assert isinstance(m.p1.bounds, model.Bounds)
    assert m.p1.bounds.min == 1.0
    assert m.p1.bounds.max == 10.0
    assert isinstance(m.p2.bounds, model.Bounds)
    assert m.p2.bounds.min == 1.0
    assert m.p2.bounds.max == 9.5
    assert isinstance(m.p4.bounds, model.Bounds)
    assert m.p4.bounds.min == 1.0
    assert m.p4.bounds.max == 8.5
    assert isinstance(m.p5, model.ConstValue)
    assert get_name(m.p5) == "p5"
    assert m.p5 == 5.0
    assert m.p5.bounds.min == 0.0
    assert m.p5.bounds.max == 10.0
    m.p4.uncertainty()
    assert m.p4.bounds is None

def test_par_in_rates2():
    """test assignment of 'local' parameters with bounds"""
    m = Model("My first model")
    m.v1 = react("A->B", " p2*A/(p1+A)-B ", pars={'p1':4})
    m.p2 = 3.0
    m.v1.p1 = 1,10 #tuple or list
    m.p2 = [1, 9.5]
    m.p3 = 5
    assert m.v1.p1 == 4.0
    assert m.p2 == 3.0
    assert m.p3 == 5.0
    assert m.p3.bounds is None
    assert isinstance(m.v1.p1.bounds, model.Bounds)
    assert m.v1.p1.bounds.min == 1.0
    assert m.v1.p1.bounds.max == 10.0
    assert isinstance(m.p2.bounds, model.Bounds)
    assert m.p2.bounds.min == 1.0
    assert m.p2.bounds.max == 9.5

def test_transf1():
    """test transf(int or float)"""
    m = Model("My first model")
    m.t1 = transf(4)
    m.t2 = transf(2.0)
    assert isinstance(m.t1, model.Transformation)
    assert isinstance(m.t2, model.Transformation)
    assert get_name(m.t1) == 't1'
    assert get_name(m.t2)== 't2'
    assert m.t1()== str(float(4))
    assert m.t2()== str(float(2.0))
    check, msg = m.checkRates()
    assert check 

def test_transf2():
    """test transf(string)"""
    m = Model("My first model")
    m.v1 = react("A+B -> C"  , 3)
    m.t1 = transf(" p2*A/(p1+A)-B ", dict(p2=3))
    m.p1 = 2
    assert isinstance(m.t1, model.Transformation)
    assert get_name(m.t1) == 't1'
    assert m.t1()== "p2*A/(p1+A)-B"
    check, msg = m.checkRates()
    print msg
    assert check 
    assert isinstance(m.t1.p2, model.ConstValue)
    assert get_name(m.t1.p2) == "p2"
    assert m.t1.p2 == 3.0

def test_printmodel():
    """test print(model)"""
    import math
    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.v4 = react("B   ->  "  , "2*B")
    m.t1 = transf("A*4 + C")
    m.t2 = transf("sqrt(2*A)")
    m.D  = variable("-2 * D")
    m.B  = 2.2
    m.myconstant = 2 * m.B / 1.1 # should be 4.0
    m.V3 = 0.5
    m.V3 = [0.1, 1.0]
    m.Km3 = 4
    m.init = state(A = 1.0, C = 1, D = 1)
    m.afterwards = state(A = 1.0, C = 2, D = 1)
    m.afterwards.C.uncertainty(1,3)
    #print should not raise an Exception
    print m

def test_clonemodel():
    """test model.clone()"""
    import math
    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.v4 = react("B   ->  "  , "2*B")
    m.t1 = transf("A*4 + C")
    m.t2 = transf("sqrt(2*A)")
    m.D  = variable("-2 * D")
    m.B  = 2.2
    m.myconstant = 2 * m.B / 1.1 # should be 4.0
    m.V3 = 0.5
    m.V3 = [0.1, 1.0]
    m.Km3 = 4
    m.init = state(A = 1.0, C = 1, D = 1)
    m.afterwards = state(A = 1.0, C = 2, D = 1)
    m.afterwards.C.uncertainty(1,3)
    m2 = m.clone()
    assert m == m2

def test_init1():
    """test assignment of states"""
    m = Model("My first model")
    m.p1 = 4
    m.p2 = 3.0
    m.init = state(x = 1, y = 2.0)
    m.end = state(x = 2, y = 4.0)
    assert isinstance(m.init, model.StateArray)
    assert isinstance(m.end, model.StateArray)
    assert m.init.x == 1.0
    assert m.init.y == 2.0

def test_iter_reactions():
    """test iteration of reactions using reactions()"""
    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = 1.0)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.v4 = react("B   ->  "  , "2*B")
    m.D  = variable("-2 * D")
    rr = reactions(m)
    assert isinstance(rr, list)
    assert len(rr) == 5
    names = [get_name(v) for v in reactions(m)]
    rates = [v() for v in reactions(m)]
    reags = [v._reagents for v in reactions(m)]
    assert names[0] == 'v1'
    assert names[1] == 'v2'
    assert names[2] == 'v3'
    assert names[3] == 'v4'
    assert names[4] == 'd_D_dt'
    assert rates[0] == '3.0*A*B'
    assert rates[3] == '2*B'
    assert reags[0][0][0] == 'A'
    assert reags[0][0][1] == 1.0
    assert reags[0][1][0] == 'B'
    assert reags[0][1][1] == 1.0
    assert len(reags[1]) == 0
    assert len(reags[2]) == 1
    assert len(reags[3]) == 1

def test_iter_transf():
    """test iteration of transformations using transformations()"""
    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = 1.0)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.D  = variable("-2 * D")
    m.t1 = transf("A*4 + C")
    m.t2 = transf("sqrt(2*A)")
    rr = transformations(m)
    assert isinstance(rr, list)
    assert len(rr) == 2
    names = [get_name(v) for v in transformations(m)]
    rates = [v() for v in transformations(m)]
    assert names[0] == 't1'
    assert names[1] == 't2'
    assert rates[0] == 'A*4 + C'
    assert rates[1] == 'sqrt(2*A)'

def test_iter_variables():
    """test iteration of variables using variables() and varnames()"""
    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = 1.0)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.D  = variable("-2 * D")
    xx = varnames(m)
    assert isinstance(xx, list)
    assert len(xx) == 4
    assert xx == ['A', 'B', 'C', 'D']

def test_iter_extvariables():
    """test iteration of external variables using extvariables()"""
    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = 1.0)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.B  = 0.5
    xx = extvariables(m)
    assert isinstance(xx, list)
    assert len(xx) == 1
    assert xx == ['B']

def test_iter_parameters():
    """test iteration of parameters using parameters()"""
    import math
    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.v4 = react("B   ->  "  , "2*B")
    m.D  = variable("-2 * D")
    m.B  = 2.2
    m.myconstant = 2 * m.B / 1.1 # should be 4.0
    m.V3 = 0.5
    m.V3 = [0.1, 1.0]
    m.Km3 = 4
    m.init = state(A = 1.0, C = 1, D = 1)
    m.afterwards = state(A = 1.0, C = 2, D = 1)
    m.afterwards.C.uncertainty(1,3)
    pp = parameters(m)
    assert isinstance(pp, list)
    assert len(pp) == 4
    names = [get_name(x) for x in parameters(m)]
    assert names.sort() == ['B', 'myconstant','Km3', 'V3'].sort()
    values = [x for x in parameters(m)]
    values.sort()
    should_values = [2.2, 4.0, 0.5, 4.0]
    should_values.sort()
    for v1,v2 in zip(values, should_values):
        assert_almost_equal(v1, v2)

def test_iter_uncertain():
    """test iteration of uncertain parameters using uncertain()"""
    import math
    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.v4 = react("B   ->  "  , "2*B")
    m.D  = variable("-2 * D")
    m.B  = 2.2
    m.myconstant = 2 * m.B / 1.1 # should be 4.0
    m.V3 = 0.5
    m.V3 = [0.1, 1.0]
    m.Km3 = 4
    m.Km3.uncertainty(0,5)
    m.init = state(A = 1.0, C = 1, D = 1)
    m.afterwards = state(A = 1.0, C = 2, D = 1)
    m.afterwards.C.uncertainty(1,3)
    uu = uncertain(m)
    assert isinstance(uu, list)
    assert len(uu) == 3
    names = [get_name(x) for x in uncertain(m)]
    for n in ['V3', 'Km3', 'afterwards.C']:
        assert n in names
    should_values = {'V3':(0.1, 1.0), 'Km3':(0.0,5.0), 'afterwards.C':(1.0,3.0)}
    for i in range(len(uu)):
        assert_almost_equal(uu[i].min, should_values[names[i]][0])
        assert_almost_equal(uu[i].max, should_values[names[i]][1])

def test_iter_state():
    """test iteration of states"""
    import math
    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.v4 = react("B   ->  "  , "2*B")
    m.D  = variable("-2 * D")
    m.B  = 2.2
    m.myconstant = 2 * m.B / 1.1 # should be 4.0
    m.V3 = 0.5
    m.V3 = [0.1, 1.0]
    m.Km3 = 4
    m.Km3.uncertainty(0,5)
    m.init = state(A = 1.0, C = 1, D = 1)
    m.afterwards = state(C = 2, D = 1)
    m.afterwards.C.uncertainty(1,3)
    initnames = [i[0] for i in m.init]
    assert len(initnames) == len(varnames(m))
    afternames = [i[0] for i in m.afterwards]
    assert len(afternames) == 2
    initlist = [(i,j) for i,j in m.init]
    assert ('A', 1.0) in m.init
    assert ('D', 1.0) in m.init

def test_reassignment1():
    """test reassignment of parameters with bounds"""
    m = Model("My first model")
    m.p1 = 4
    m.p2 = 3.0
    m.p1 = 1,10 #tuple or list
    m.p2 = [1, 9.5]
    m.p3 = 5
    m.p4 = 6
    m.p4.uncertainty(1, 8.5) # or uncertainty function
    m.p5 = 0,10 # bounds create parameter with midpoint
    # start reassignments
    m.p1 = 5
    assert m.p1 == 5.0
    assert isinstance(m.p1.bounds, model.Bounds)
    assert m.p1.bounds.min == 1.0
    assert m.p1.bounds.max == 10.0

    m.p3 = 0,3
    assert m.p3 == 5.0
    assert isinstance(m.p3.bounds, model.Bounds)
    assert m.p3.bounds.min == 0.0
    assert m.p3.bounds.max == 3.0
    
    m.p4 = 1,10
    assert m.p4 == 6.0
    assert isinstance(m.p4.bounds, model.Bounds)
    assert m.p4.bounds.min == 1.0
    assert m.p4.bounds.max == 10.0

def test_reassignment2():
    """test reassignment of reactions"""
    m = Model("My first model")
    m.v1 = react("A->B", 4)
    m.v2 = react("B->C", 2.0)
    assert isinstance(m.v1, model.Reaction)
    assert isinstance(m.v2, model.Reaction)
    assert get_name(m.v1) == 'v1'
    assert get_name(m.v2) == 'v2'
    assert m.v1()== str(float(4))+ "*A"
    assert m.v2()== str(float(2.0))+"*B"
    check, msg = m.checkRates()
    assert check 
    m.v2 = react("D->C", 2.0)
    assert isinstance(m.v1, model.Reaction)
    assert isinstance(m.v2, model.Reaction)
    assert get_name(m.v1) == 'v1'
    assert get_name(m.v2) == 'v2'
    assert m.v1()== str(float(4))+ "*A"
    assert m.v2()== str(float(2.0))+"*D"
    check, msg = m.checkRates()
    assert check 

def test_reassignment3():
    """test change of variables by reassignment of reactions"""
    m = Model("My first model")
    m.v1 = react("A->B", 4)
    m.v2 = react("B->C", 2.0)
    xx = varnames(m)
    assert len(xx) == 3
    assert xx == ['A', 'B', 'C']
    check, msg = m.checkRates()
    assert check 
    m.v2 = react("B->D", 2.0)
    xx = varnames(m)
    assert len(xx) == 3
    assert xx == ['A', 'B', 'D']
    check, msg = m.checkRates()
    assert check 

@raises(model.BadTypeComponent)
def test_reassignment4():
    """test illegal type reassignment (reactions)"""
    m = Model("My first model")
    m.v1 = react("A->B", 4)
    m.v2 = react("B->C", 2.0)
    xx = varnames(m)
    assert len(xx) == 3
    assert xx == ['A', 'B', 'C']
    check, msg = m.checkRates()
    assert check 
    m.v2 = 3.14
    xx = varnames(m)
    assert len(xx) == 3
    assert xx == ['A', 'B', 'D']
    check, msg = m.checkRates()
    assert check 

@raises(model.BadTypeComponent)
def test_reassignment5():
    """test illegal type reassignment (parameters)"""
    m = Model("My first model")
    m.v1 = react("A->B", 4)
    m.v2 = react("B->C", 2.0)
    m.Km = 4
    check, msg = m.checkRates()
    assert check 
    m.Km = react("B->D", 2.0)
    xx = variables(m)
    assert len(xx) == 3
    names = [x for x in varnames(m)]
    assert names == ['A', 'B', 'D']
    check, msg = m.checkRates()
    assert check 

@raises(model.BadTypeComponent)
def test_par2():
    """test illegal type assignment"""
    m = Model("My first model")
    m.v1 = react("A->B", 4)
    m.v2 = react("B->C", 2.0)
    m.Km = [9,10,13,45]
    check, msg = m.checkRates()
    assert check 

def test_meta1():
    """test Model metadata"""
    m = Model("My first model")
    m.v1 = react("A->B", 4)
    m.v2 = react("B->C", 2.0)
    check, msg = m.checkRates()
    assert check 
    m['where'] = 'in model'
    m['for what'] = 'testing'
    assert m['where'] == 'in model'
    assert m['for what'] == 'testing'
    assert m['title'] == 'My first model'
    assert m['nonexistent'] is None
    del m['where']
    assert m['where'] is None
