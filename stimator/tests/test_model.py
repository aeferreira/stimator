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
    assert m.title == "My first model"

def test_react1():
    """test react(string, int or float)"""
    m = Model("My first model")
    m.v1 = react("A->B", 4)
    m.v2 = react("B->C", 2.0)
    assert isinstance(m.v1, model.Reaction)
    assert isinstance(m.v2, model.Reaction)
    assert m.v1.name == 'v1'
    assert m.v2.name == 'v2'
    assert m.v1.rate== str(float(4))+ "*A"
    assert m.v2.rate== str(float(2.0))+"*B"
    check, msg = m.checkRates()
    assert check 

def test_react2():
    """test react(string, string)"""
    m = Model("My first model")
    m.v1 = react("A->B", " 4*A/(p1+A)-B ")
    m.p1 = 2
    assert m.v1.name == 'v1'
    assert isinstance(m.v1, model.Reaction)
    assert m.v1.rate== "4*A/(p1+A)-B"
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
    assert m.v1.name == 'v1'
    assert isinstance(m.v1, model.Reaction)
    assert m.v1.rate== "4*A/(p2+A)-B"
    check, msg = m.checkRates()
    assert not check 

def test_react5():
    """test Bad rate law (malformed expression)"""
    m = Model("My first model")
    m.v1 = react("A->B", " 4*A/(p1+A-B ")
    m.p1 = 2
    assert m.v1.name == 'v1'
    assert isinstance(m.v1, model.Reaction)
    assert m.v1.rate== "4*A/(p1+A-B"
    check, msg = m.checkRates()
    assert not check 

def test_react6():
    """test Bad rate law (fp overflow)"""
    m = Model("My first model")
    m.v1 = react("A->B", " 1e100**10000 * 4*A/(p1+A)-B ")
    m.p1 = 2
    assert m.v1.name == 'v1'
    assert isinstance(m.v1, model.Reaction)
    assert m.v1.rate== "1e100**10000 * 4*A/(p1+A)-B"
    check, msg = m.checkRates()
    assert not check 

def test_par1():
    """test assignment of parameters"""
    m = Model("My first model")
    m.p1 = 4
    m.p2 = 3.0
    assert isinstance(m.p1, model.ConstValue)
    assert m.p1.name == "p1"
    assert isinstance(m.p2, model.ConstValue)
    assert m.p2.name == "p2"
    assert m.p1 == 4.0
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
    assert m.p5.name == "p5"
    assert m.p5 == 5.0
    assert m.p5.bounds.min == 0.0
    assert m.p5.bounds.max == 10.0
    m.p4.uncertainty()
    assert m.p4.bounds is None

def test_transf1():
    """test transf(int or float)"""
    m = Model("My first model")
    m.t1 = transf(4)
    m.t2 = transf(2.0)
    assert isinstance(m.t1, model.Transformation)
    assert isinstance(m.t2, model.Transformation)
    assert m.t1.name == 't1'
    assert m.t2.name== 't2'
    assert m.t1.rate== str(float(4))
    assert m.t2.rate== str(float(2.0))
    check, msg = m.checkRates()
    assert check 

def test_transf2():
    """test transf(string)"""
    m = Model("My first model")
    m.v1 = react("A+B -> C"  , 3)
    m.t1 = transf(" 4*A/(p1+A)-B ")
    m.p1 = 2
    assert isinstance(m.t1, model.Transformation)
    assert m.t1.name == 't1'
    assert m.t1.rate== "4*A/(p1+A)-B"
    check, msg = m.checkRates()
    print msg
    assert check 

def test_printmodel():
    """test print(model)"""
    import math
    def force(A, t):
        return t*A
    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.v4 = react("B   ->  "  , "2*input1")
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
    m.input1 = forcing(force)
    #print should not raise an Exception
    print m

def test_clonemodel():
    """test model.clone()"""
    import math
    def force(A, t):
        return t*A
    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.v4 = react("B   ->  "  , "2*input1")
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
    m.input1 = forcing(force)
    mstr = str(m)
    m2 = m.clone()
    m2str = str(m2)
    assert mstr == m2str

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

