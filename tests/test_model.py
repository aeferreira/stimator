import pytest
from stimator import *

def assert_almost_equal(x, y):
    if abs(x-y) < 0.0001:
        return True
    return False

def conservation(total, A):
    return total - A

def test_Model_init_1():
    """test Model __init__()"""
    m = Model()
    assert isinstance(m, Model)
    m = Model("My first model")
    assert isinstance(m, Model)
    assert m.metadata['title'] == "My first model"
    assert m.name == "My first model"

def test_Model_register_kin_func():
    """test Model register_kin_func()"""
    m = Model("My first model")
    assert isinstance(m, Model)
    before = len(m._usable_functions)
    m.register_kin_func(conservation)
    assert len(m._usable_functions) == before + 1
    assert callable(m._usable_functions['conservation'])

def test_set_reaction1():
    """test Model.set_reaction(string, string, int or float)"""
    m = Model("My first model")
    m.set_reaction('v1',"A->B", 4)
    m.set_reaction('v2',"B->C", 2.0)
    m.set_reaction('v3',"2C->D", 2.0)
    m.set_reaction('v4',"0C+2D->E", 2.0)
    m.set_reaction('v5',"->E", 2.0)
    assert isinstance(m.reactions.v1, model.Reaction)
    assert isinstance(m.reactions.v2, model.Reaction)
    assert isinstance(m.reactions.v3, model.Reaction)
    assert isinstance(m.reactions.v4, model.Reaction)
    assert isinstance(m.reactions.v5, model.Reaction)
    assert m.reactions.v1.name == 'v1'
    assert m.reactions.v2.name == 'v2'
    assert m.reactions.v3.name == 'v3'
    assert m.reactions.v4.name == 'v4'
    assert m.reactions.v5.name == 'v5'
    cd = {'A': 3.14, 'B': 2*3.14, 'C': 3.0*3.14, 'D': 4*3.14}
    assert_almost_equal(eval('4.0*A', cd), eval(m.reactions.v1(), cd))
    assert_almost_equal(eval('2.0*B', cd), eval(m.reactions.v2(), cd))
    assert_almost_equal(eval('2.0*C**2.0', cd), eval(m.reactions.v3(), cd))
    assert_almost_equal(eval('2.0*D**2', cd), eval(m.reactions.v4(), cd))
    assert_almost_equal(eval('2.0', cd), eval(m.reactions.v5(), cd))
    check, msg = m.checkRates()
    assert check
    assert 'v1' in m.reactions
    assert 'x1' not in m.reactions

def test_set_reaction1b():
    """test Model.set_reaction() and Reaction properties"""
    m = Model("My first model")
    m.set_reaction('v1',"2 A + 3.2 B -> C + 4 D + 0.5 E", 4)
    m.set_reaction('v2',"B->C", 2.0)
    m.set_reaction('v3',"2C<=>D", 2.0)
    m.set_reaction('v4',"0C+2D->E", 2.0)
    name = m.reactions.v1.name
    reagents = m.reactions.v1.reagents
    products = m.reactions.v1.products
    stoichiometry = m.reactions.v1.stoichiometry
    stoichiometry_string = m.reactions.v1.stoichiometry_string
    stoichiometry = m.reactions.v3.stoichiometry
    stoichiometry_string = m.reactions.v3.stoichiometry_string
    stoichiometry = m.reactions.v4.stoichiometry
    stoichiometry_string = m.reactions.v4.stoichiometry_string
    check, msg = m.checkRates()
    assert check 

def test_set_reaction2():
    """test Model.set_reaction(string, string, string)"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", " p2*A/(p1+A)-B ", {'p2': 4})
    m.setp('p1', 2)
    assert (m.reactions.v1.name) == 'v1'
    assert isinstance(m.reactions.v1, model.Reaction)
    assert m.reactions.v1()== "p2*A/(p1+A)-B"
    assert m.reactions.v1(fully_qualified=True) == "v1.p2*A/(p1+A)-B"
    check, msg = m.checkRates()
    assert check 

def test_set_reaction_absent():
    """test bad attribute of Model.reactions"""
    with pytest.raises(AttributeError):
        m = Model("My first model")
        m.set_reaction('v1', "A->B", " 4*A/(p1+A)-B ")
        m.setp('p1', 2)
        assert (m.reactions.v2.name) == 'v1'
        assert isinstance(m.reactions.v2, model.Reaction)
    
def test_set_reaction2b():
    """test Model.react(string, string) with math functions"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", " 4*sqrt(A)/(p1+sin(A))-B ")
    m.setp('p1', 2)
    assert (m.reactions.v1.name) == 'v1'
    assert isinstance(m.reactions.v1, model.Reaction)
    assert m.reactions.v1() == "4*sqrt(A)/(p1+sin(A))-B"
    check, msg = m.checkRates()
    assert check 

def test_set_reaction2c():
    """test Model.react(string, string) with kinetics functions"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", " 4*A*step(t,1.0)")
    m.setp('p1', 2)
    assert (m.reactions.v1.name) == 'v1'
    assert isinstance(m.reactions.v1, model.Reaction)
    assert m.reactions.v1() == "4*A*step(t,1.0)"
    check, msg = m.checkRates()
    assert check 


def test_set_reaction3():
    """test Bad stoichiometry"""
    with pytest.raises(model.BadStoichError):
        m = Model("My first model")
        m.set_reaction('v1', "A ## B", " 4*A/(p1+A)-B ")

def test_set_reaction3b():
    """test Bad stoichiometry again"""
    with pytest.raises(model.BadStoichError):
        m = Model("My first model")
        m.set_reaction('v1', "A->##B", " 4*A/(p1+A)-B ")

def test_set_reaction4():
    """test Bad rate law (unknown ID)"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", " 4*A/(p2+A)-B ")
    m.setp('p1',2)
    assert (m.reactions.v1.name) == 'v1'
    assert isinstance(m.reactions.v1, model.Reaction)
    assert m.reactions.v1() == "4*A/(p2+A)-B"
    check, msg = m.checkRates()
    assert not check 

def test_set_reaction5():
    """test Bad rate law (malformed expression)"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", " 4*A/(p1+A-B ")
    m.setp('p1',2)
    assert (m.reactions.v1.name) == 'v1'
    assert isinstance(m.reactions.v1, model.Reaction)
    assert m.reactions.v1() == "4*A/(p1+A-B"
    check, msg = m.checkRates()
    assert not check 

def test_set_reaction6():
    """test Bad rate law (fp overflow)"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", " 1e100**10000 * 4*A/(p1+A)-B ")
    m.setp('p1',2)
    assert (m.reactions.v1.name) == 'v1'
    assert isinstance(m.reactions.v1, model.Reaction)
    assert m.reactions.v1() == "1e100**10000 * 4*A/(p1+A)-B"
    check, msg = m.checkRates()
    assert not check 

def test_par1():
    """test assignment of parameters"""
    m = Model("My first model")
    m.setp('p1', 4)
    m.parameters.p2 = 3.0
    m.parameters.p3 = 3.0
    m.setp({'p10': 0, 'p11':1})
    m.setp([('p12', 2), ('p13', 3)])
    assert isinstance(m.parameters.p1, model.ConstValue)
    assert (m.parameters.p1.name) == "p1"
    assert isinstance(m.parameters.p2, model.ConstValue)
    assert (m.parameters.p2.name) == "p2"
    assert m.parameters.p1 == 4.0
    assert m.parameters.p2 == 3.0
    assert m.parameters.p10 == 0.0
    assert m.parameters.p11 == 1.0
    assert m.parameters.p12 == 2.0
    assert m.parameters.p13 == 3.0
    assert m.getp('p1') == 4.0
    assert 'p1' in m.parameters
    assert 'p8' not in m.parameters

def test_par1b():
    """test assignment and deletion of parameters"""
    m = Model("My first model")
    m.setp('p1', 4)
    m.parameters.p2 = 3.0
    m.parameters.p3 = 3.0
    with pytest.raises(AttributeError):
        assert isinstance(m.parameters.p1, model.ConstValue)
        assert (m.parameters.p1.name) == "p1"
        assert isinstance(m.parameters.p2, model.ConstValue)
        assert (m.parameters.p2.name) == "p2"
        assert m.parameters.p1 == 4.0
        assert m.parameters.p2 == 3.0
        m.setp('p1', None)
        assert m.getp('p1') is None

def test_par1c():
    """test assignment of bad parameter values"""
    m = Model("My first model")
    m.setp('p1', 4)
    with pytest.raises(model.BadTypeComponent):
        assert isinstance(m.parameters.p1, model.ConstValue)
        assert (m.parameters.p1.name) == "p1"
        assert m.parameters.p1 == 4.0
        m.setp('p2', 'bb')

def test_par1d():
    """test assignment of new parameters to None"""
    m = Model("My first model")
    m.setp('p1', 4)
    with pytest.raises(model.BadTypeComponent):
        assert isinstance(m.parameters.p1, model.ConstValue)
        assert (m.parameters.p1.name) == "p1"
        assert m.parameters.p1 == 4.0
        m.setp('p2', None)

def test_par_in_rates1():
    """test assignment of parameters 'local' to reactions"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", " p2*A/(p1+A)-B ", pars={'p1':4})
    m.setp('p2', 3.0)
    m.setp('v1.p1', 5)
    assert (m.reactions.v1.name) == 'v1'
    assert isinstance(m.reactions.v1, model.Reaction)
    assert m.reactions.v1()== "p2*A/(p1+A)-B"
    check, msg = m.checkRates()
    assert check 
    assert isinstance(m.parameters.v1.p1, model.ConstValue)
    assert (m.parameters.v1.p1.name) == "p1"
    assert isinstance(m.parameters.p2, model.ConstValue)
    assert (m.parameters.p2.name) == "p2"
    assert m.parameters.v1.p1 == 5.0
    assert m.parameters.p2 == 3.0
    assert m.getp('v1.p1') == 5.0

def test_par_in_rates1b():
    """testing absent parameters 'local' to reactions"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", " p2*A/(p1+A)-B ", pars={'p1':4})
    m.setp('p2', 3.0)
    m.setp('v1.p1', 5)
    check, msg = m.checkRates()
    assert check 
    with pytest.raises(AttributeError):
        assert (m.parameters.v1.p1.name) == "p1"
        assert (m.parameters.p2.name) == "p2"
        assert m.parameters.v1.p1 == 5.0
        assert m.parameters.v1.p2 == 5.0


def test_par_in_rates1c():
    """testing absent parameters 'local' to reactions with getp()"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", " p2*A/(p1+A)-B ", pars={'p1':4})
    m.setp('p2', 3.0)
    m.setp('v1.p1', 5)
    check, msg = m.checkRates()
    assert check 
    with pytest.raises(AttributeError):
        assert (m.parameters.v1.p1.name) == "p1"
        assert (m.parameters.p2.name) == "p2"
        assert m.parameters.v1.p1 == 5.0
        assert m.getp('v1.p2') == 5.0


def test_par_in_rates1d():
    """testing bad parameters 'local' to reactions with getp()"""
    m = Model("My first model")
    with pytest.raises(TypeError):
        m.set_reaction('v1', "A->B", " p2*A/(p1+A)-B ", pars={'p1':'bb'})
        m.setp('p2', 3.0)
        m.setp('v1.p1', 5)
        check, msg = m.checkRates()
        assert check 
        assert (m.parameters.v1.p1.name) == "p1"
        assert (m.parameters.p2.name) == "p2"
        assert m.parameters.v1.p1 == 5.0
        assert m.getp('v1.p2') == 5.0

def test_par_from_rates1():
    """test rate expressions with parameters 'local' to reactions"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", " p2*A/(p1+A)-B ", pars=(('p1', 4.0), ('p2',2)))
    m.set_reaction('v2', "B->C", "2*v1.p1*B")
    m.setp('p2', 3.0)
    assert (m.reactions.v2.name) == 'v2'
    assert isinstance(m.reactions.v1, model.Reaction)
    assert m.reactions.v2()== "2*v1.p1*B"
    check, msg = m.checkRates()
    assert check
    assert len(m.parameters.v1) == 2
    assert 'p1' in m.parameters.v1
    assert 'q2' not in m.parameters.v1
    pars_v1 = m.reactions.v1.parameters
    assert len(pars_v1) == 2
    pars_v1_iter = [par for par in m.reactions.v1]
    assert len(pars_v1_iter) == 2
    assert isinstance(m.parameters.v1.p1, model.ConstValue)
    assert (m.parameters.v1.p1.name) == "p1"
    assert isinstance(m.parameters.p2, model.ConstValue)
    assert (m.parameters.p2.name) == "p2"
    assert m.parameters.v1.p1 == 4.0
    assert m.parameters.p2 == 3.0

def test_bounds():
    """test assignment of parameters with bounds"""
    m = Model("My first model")
    m.parameters.p1 = 4
    m.parameters.p2 = 3.0
    m.parameters.p1.set_bounds((1,10)) #tuple or list
    m.set_bounds('p2', [1, 9.5])
    m.parameters.p3 = 5
    m.parameters.p4 = 6
    m.set_bounds('p4',(1, 8.5)) # or uncertainty function
    m.setp('p5',5)
    m.parameters.p5.bounds = model.Bounds('?',0,10)
    m.setp('p6',6)
    assert m.parameters.p1 == 4.0
    assert m.parameters.p2 == 3.0
    assert m.parameters.p3 == 5.0
    assert m.parameters.p4 == 6.0
    assert m.parameters.p3.bounds is None
    assert isinstance(m.parameters.p1.bounds, model.Bounds)
    assert m.parameters.p1.bounds.lower == 1.0
    assert m.parameters.p1.bounds.upper == 10.0
    bb = m.parameters.p2.get_bounds()
    assert bb == (1.0, 9.5)
    bb = m.parameters.p6.get_bounds()
    assert bb is None
    assert isinstance(m.parameters.p4.bounds, model.Bounds)
    assert m.parameters.p4.bounds.lower == 1.0
    assert m.parameters.p4.bounds.upper == 8.5
    assert isinstance(m.parameters.p5, model.ConstValue)
    assert (m.parameters.p5.name) == "p5"
    assert m.parameters.p5 == 5.0
    assert m.parameters.p5.bounds.lower == 0.0
    assert m.parameters.p5.bounds.upper == 10.0
    m.parameters.p4.reset_bounds()
    assert m.parameters.p4.bounds is None
    m.parameters.p5.set_bounds(None)
    assert m.parameters.p5.bounds is None
    bb = m.get_bounds('p5')
    assert bb is None
    m.reset_bounds('p2')
    bb = m.get_bounds('p2')
    assert bb is None
    # this creates a new parameter
    m.set_bounds('p7', (0,10))
    assert m.parameters.p7 == 5.0
    assert m.parameters.p7.bounds.lower == 0.0
    assert m.parameters.p7.bounds.upper == 10.0

def test_bounds2():
    """test assignment of parameters wrong bounds"""
    m = Model("My first model")
    m.parameters.p1 = 4
    m.parameters.p2 = 3.0
    with pytest.raises(model.BadTypeComponent):    
        m.parameters.p1.set_bounds((1,10,5))
        assert m.parameters.p1.bounds is None

def test_par_in_rates2():
    """test assignment of 'local' parameters with bounds"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", " p1*A/(q1+A)-B ", pars={'p1':4, 'q1':5})
    m.setp('p2', 3.0)
    m.parameters.v1.p1.set_bounds((1,10)) #tuple or list
    m.set_bounds('v1.q1', [2, 5]) #tuple or list
    m.parameters.p2.set_bounds([1, 9.5])
    m.parameters.p3 = 5
    assert m.parameters.v1.p1 == 4.0
    assert m.parameters.p2 == 3.0
    assert m.parameters.p3 == 5.0
    assert m.parameters.p3.bounds is None
    assert isinstance(m.parameters.v1.p1.bounds, model.Bounds)
    assert isinstance(m.parameters.v1.q1.bounds, model.Bounds)
    assert m.parameters.v1.p1.bounds.lower == 1.0
    assert m.parameters.v1.p1.bounds.upper == 10.0
    assert m.parameters.v1.q1.bounds.lower == 2.0
    assert m.parameters.v1.q1.bounds.upper == 5.0
    bb = m.get_bounds('v1.q1')
    assert bb == (2,5)
    assert isinstance(m.parameters.p2.bounds, model.Bounds)
    assert m.parameters.p2.bounds.lower == 1.0
    assert m.parameters.p2.bounds.upper == 9.5
    m.reset_bounds('v1.q1')
    assert m.parameters.v1.q1.bounds is None

def test_transf1():
    """test transf(int or float)"""
    m = Model("My first model")
    m.set_transformation('t1',4)
    m.set_transformation('t2',2.0)
    assert isinstance(m.transformations.t1, model.Transformation)
    assert isinstance(m.transformations.t2, model.Transformation)
    assert (m.transformations.t1.name) == 't1'
    assert (m.transformations.t2.name)== 't2'
    assert m.transformations.t1() == str(float(4))
    assert m.transformations.t2() == str(float(2.0))
    check, msg = m.checkRates()
    assert check 

def test_transf2():
    """test transf(string)"""
    m = Model("My first model")
    m.set_reaction('v1', "A+B -> C", 3)
    m.set_transformation('t1', " p2*A/(p1+A)-B ", dict(p2=3))
    m.parameters.p1 = 2
    assert isinstance(m.transformations.t1, model.Transformation)
    assert (m.transformations.t1.name) == 't1'
    assert m.transformations.t1() == "p2*A/(p1+A)-B"
    check, msg = m.checkRates()
    assert check 
    assert isinstance(m.parameters.t1.p2, model.ConstValue)
    assert (m.parameters.t1.p2.name) == "p2"
    assert m.parameters.t1.p2 == 3.0
    assert m.getp('t1.p2') == 3.0
    m.setp('t1.p2', 5.0)
    assert m.getp('t1.p2') == 5.0

def test_set_input_var():
    """test Model.set_input_var()"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", " 4*A/(p1+A)-B ")
    m.set_input_var('p1', '3*A')
    assert (m.input_variables.p1.name) == 'p1'
    assert isinstance(m.input_variables.p1, model.Input_Variable)
    assert m.reactions.v1()== "4*A/(p1+A)-B"
    assert m.input_variables.p1()== "3*A"
    m.set_input_var('p2', 2.0)
    assert (m.input_variables.p2.name) == 'p2'
    assert isinstance(m.input_variables.p2, model.Input_Variable)
    assert m.input_variables.p2()== "2.0"
    check, msg = m.checkRates()
    assert check

def test_printmodel():
    """test print(model)"""
    import math
    m = Model('My first model')
    m.set_reaction('v1', "A+B -> C"  , 3)
    m.set_reaction('v2', "    -> A"  , rate = math.sqrt(4.0)/2)
    m.set_reaction('v3', "C   ->  "  , "V3 * C / (Km3 + C)", {'Km3':4})
    m.set_reaction('v4', "B   ->  "  , "2*B")
    m.set_transformation('t1', "k1*A*4 + C", {'k1':3.14})
    m.set_transformation('t2', "sqrt(2*A)")
    m.set_variable_dXdt('D',"-2 * D")
    m.parameters.B  = 2.2
    m.parameters.myconstant = 2 * m.parameters.B / 1.1 # should be 4.0
    m.parameters.V3 = 0.5
    m.parameters.V3.set_bounds([0.1, 1.0])
    m.set_init(A = 1.0, C = 1, D = 1)
    m.init.C.set_bounds((1,3))
    # print() should not raise an Exception
    print (m)

def test_clonemodel():
    """test model.copy()"""
    import math
    m = Model('My first model')
    m.set_reaction('v1', "A+B -> C"  , 3)
    m.set_reaction('v2', "    -> A"  , rate = math.sqrt(4.0)/2)
    m.set_reaction('v3', "C   ->  "  , "V3 * C / (Km3 + C)")
    m.set_reaction('v4', "B   ->  "  , "2*B")
    m.set_transformation('t1', "A*4 + C")
    m.set_transformation('t2', "sqrt(2*A)")
    m.set_variable_dXdt('D',"-2 * D")
    m.set_input_var('p1', '3*A')
    m.parameters.B  = 2.2
    m.parameters.myconstant = 2 * m.parameters.B / 1.1 # should be 4.0
    m.parameters.V3 = 0.5
    m.parameters.V3.set_bounds([0.1, 1.0])
    m.parameters.Km3 = 4
    m.set_init(A = 1.0, C = 1, D = 1)
    m.init.C.set_bounds((1,3))
    m2 = m.copy()
    assert m2 == m

def test_model__eq__():
    """test model.__eq__"""
    import math
    m = Model('My first model')
    m.set_reaction('v1', "A+B -> C"  , 3)
    m.set_reaction('v2', "    -> A"  , rate = math.sqrt(4.0)/2)
    m.set_reaction('v3', "C   ->  "  , "V3 * C / (Km3 + C)", {'Km3': 4})
    m.set_reaction('v4', "B   ->  "  , "2*B")
    m.set_transformation('t1', "A*4 + C")
    m.set_transformation('t2', "sqrt(2*A)")
    m.set_variable_dXdt('D',"-2 * D")
    m.parameters.B  = 2.2
    m.parameters.myconstant = 2 * m.parameters.B / 1.1 # should be 4.0
    m.parameters.V3 = 0.5
    m.parameters.V3.set_bounds([0.1, 1.0])
    m.set_init(A = 1.0, C = 1, D = 1)
    m.init.C.set_bounds((1,3))
    m1 = m.copy()
    assert m1 == m
    m2 = m.copy(new_title='A different title')
    assert m2 != m
    m3 = m.copy()
    m3.set_transformation('t3', "sqrt(2*A**2)")
    assert m3 != m
    m4 = m.copy()
    m4.parameters.V3 = 0.6
    assert m4 != m
    m5 = m.copy()
    m5.parameters.V3.reset_bounds()
    assert m5 != m
    m6 = m.copy()
    m6.parameters.v3.Km3 = 5
    assert m6 != m
    m7 = m.copy()
    m7.reactions.v1.setp('newKm3', 5)
    assert m7 != m
    m8 = m.copy()
    m8.set_reaction('v1', "A+B -> 2C"  , 3)
    assert m8 != m
    m9 = m.copy()
    m9.set_reaction('v1', "D -> C"  , 3)
    assert m9 != m
    m10 = m.copy()
    m10.set_reaction('v1', "A+B -> C"  , 3.3)
    assert m10 != m
    m11 = m.copy()
    m11.metadata['!!'] = 'B C'
    assert m10 != m

def test_set_init1():
    """test set_init()"""
    m = Model("My first model")
    m.parameters.p1 = 4
    m.parameters.p2 = 3.0
    m.set_init(x = 1, y = 2.0)
    m.set_init('z', 3.0)
    m.set_init({'a': 0, 'b':1})
    m.set_init([('c', 2), ('d', 3)])
    m.init.aaa = 1
    assert m.init.x == 1.0 
    assert m.init.y == 2.0
    assert m.init.z == 3.0
    assert m.init.a == 0.0 
    assert m.init.aaa == 1.0 
    assert m.init.b == 1.0
    assert m.init.c == 2.0
    assert m.init.d == 3.0
    assert len(m.init) == 8
    assert 'x' in m.init
    assert 'e' not in m.init
    m.reset_init()
    assert m.init.x == 0.0
    assert m.init.y == 0.0

def test_getinit1():
    """test get_init() with None or sequence"""
    m = Model("My first model")
    m.parameters.p1 = 4
    m.parameters.p2 = 3.0
    m.set_init(x = 1, y = 2.0)
    result = m.get_init()
    assert result['x'] == 1.0
    assert result['y'] == 2.0
    result = m.get_init(list('xyz'))
    assert result['x'] == 1.0
    assert result['y'] == 2.0
    assert result['z'] == 0.0

def test_getinit2():
    """test get_init() with names"""
    m = Model("My first model")
    m.parameters.p1 = 4
    m.parameters.p2 = 3.0
    m.set_init(x = 1, y = 2.0)
    x = m.get_init('x')
    y = m.get_init('y')
    assert x == 1.0
    assert y == 2.0
    z = m.get_init('z')
    assert z == 0.0

def test_getinit3():
    """test get_init(), retrieving attributes"""
    m = Model("My first model")
    m.parameters.p1 = 4
    m.parameters.p2 = 3.0
    m.set_init(x = 1, y = 2.0)
    m.set_bounds('init.x', (0.8, 1.0))
    name_x = m.get_init('x').name
    bounds_x = m.get_init('x').bounds
    assert name_x == 'x'
    assert bounds_x.lower == 0.8
    assert bounds_x.upper == 1.0

def test_iter_reactions():
    """test iteration of reactions using reactions()"""
    import math
    m = Model('My first model')
    m.set_reaction('v1', "A+B -> C"  , 3)
    m.set_reaction('v2', "    -> A"  , rate = math.sqrt(4.0)/2)
    m.set_reaction('v3', "C   ->  "  , "V3 * C / (Km3 + C)")
    m.set_reaction('v4', "B   ->  "  , "2*B")
    m.set_transformation('t1', "A*4 + C")
    m.set_transformation('t2', "sqrt(2*A)")
    m.set_variable_dXdt('D',"-2 * D")
    assert len(m.reactions) == 5
    names = [v.name for v in m.reactions]
    rates = [v() for v in m.reactions]
    reags = [v._reagents for v in m.reactions]
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
    import math
    m = Model('My first model')
    m.set_reaction('v1', "A+B -> C"  , 3)
    m.set_reaction('v2', "    -> A"  , rate = math.sqrt(4.0)/2)
    m.set_reaction('v3', "C   ->  "  , "V3 * C / (Km3 + C)")
    m.set_reaction('v4', "B   ->  "  , "2*B")
    m.set_transformation('t1', "A*4 + C")
    m.set_transformation('t2', "sqrt(2*A)")
    m.set_variable_dXdt('D',"-2 * D")
    assert len(m.transformations) == 2
    names = [v.name for v in m.transformations]
    rates = [v() for v in m.transformations]
    assert names[0] == 't1'
    assert names[1] == 't2'
    assert rates[0] == 'A*4 + C'
    assert rates[1] == 'sqrt(2*A)'

def test_iter_variables():
    """test iteration of variables using and varnames"""
    import math
    m = Model('My first model')
    m.set_reaction('v1', "A+B -> C"  , 3)
    m.set_reaction('v2', "    -> A"  , rate = math.sqrt(4.0)/2)
    m.set_reaction('v3', "C   ->  "  , "V3 * C / (Km3 + C)")
    m.set_reaction('v4', "B   ->  "  , "2*B")
    m.set_transformation('t1', "A*4 + C")
    m.set_transformation('t2', "sqrt(2*A)")
    m.set_variable_dXdt('D',"-2 * D")
    xx = m.varnames
    assert isinstance(xx, list)
    assert len(xx) == 4
    assert xx == ['A', 'B', 'C', 'D']

def test_iter_extvariables():
    """test iteration of external variables using extvariables"""
    import math
    m = Model('My first model')
    m.set_reaction('v1', "A+B -> C"  , 3)
    m.set_reaction('v2', "    -> A"  , rate = math.sqrt(4.0)/2)
    m.set_reaction('v3', "C   ->  "  , "V3 * C / (Km3 + C)")
    m.set_reaction('v4', "B   ->  "  , "2*B")
    m.set_transformation('t1', "A*4 + C")
    m.set_transformation('t2', "sqrt(2*A)")
    m.set_variable_dXdt('D',"-2 * D")
    m.parameters.B  = 2.2
    m.parameters.myconstant = 2 * m.parameters.B / 1.1 # should be 4.0
    m.parameters.V3 = 0.5
    m.parameters.V3.set_bounds([0.1, 1.0])
    m.parameters.Km3 = 4
    m.set_init(A = 1.0, C = 1, D = 1)
    m.init.C.set_bounds((1,3))
    xx = m.extvariables
    assert isinstance(xx, list)
    assert len(xx) == 1
    assert xx == ['B']

def test_iter_parameters():
    """test iteration of parameters using parameters()"""
    import math
    m = Model('My first model')
    m.set_reaction('v1', "A+B -> C"  , 3)
    m.set_reaction('v2', "    -> A"  , rate = math.sqrt(4.0)/2)
    localp = {'V3': 4, 'Km3': 0.5}
    m.set_reaction('v3', "C   ->  "  , "V3 * C / (Km3 + C)", localp)
    m.set_reaction('v4', "B   ->  "  , "2*B")
    m.set_transformation('t1', "A*4 + C")
    m.set_transformation('t2', "sqrt(2*A)")
    m.set_variable_dXdt('D',"-2 * D")
    m.parameters.B  = 2.2
    m.parameters.myconstant = 2 * m.parameters.B / 1.1 # should be 4.0
    m.parameters.v3.V3 = 0.5
    m.parameters.v3.V3.set_bounds([0.1, 1.0])
    m.parameters.v3.Km3 = 4
    m.set_init(A = 1.0, C = 1, D = 1)
    pp = m.parameters
    assert len(pp) == 4
    names = [x.name for x in m.parameters]
    assert names.sort() == ['B', 'myconstant','v3.Km3', 'v3.V3'].sort()
    values = [x for x in m.parameters]
    values.sort()
    should_values = [2.2, 4.0, 0.5, 4.0]
    should_values.sort()
    for v1,v2 in zip(values, should_values):
        assert_almost_equal(v1, v2)

def test_iter_uncertain():
    """test iteration of uncertain parameters using uncertain()"""
    import math
    m = Model('My first model')
    m.set_reaction('v1', "A+B -> C"  , 3)
    m.set_reaction('v2', "    -> A"  , rate = math.sqrt(4.0)/2)
    m.set_reaction('v3', "C   ->  "  , "V3 * C / (Km3 + C)")
    m.set_reaction('v4', "B   ->  "  , "2*B")
    m.set_transformation('t1', "A*4 + C")
    m.set_transformation('t2', "sqrt(2*A)")
    m.set_variable_dXdt('D',"-2 * D")
    m.setp('B', 2.2)
    m.parameters.myconstant = 2 * m.getp('B') / 1.1 # should be 4.0
    m.parameters.V3 = 0.5
    m.parameters.V3.set_bounds([0.1, 1.0])
    m.parameters.Km3 = 4
    m.parameters.Km3.set_bounds((0,5))
    m.set_init(A = 1.0, C = 1, D = 1)
    m.init.A.set_bounds((1,3))
    uu = m.with_bounds
    assert len(uu) == 3
    names = [x.name for x in m.with_bounds]
    for n in ['V3', 'Km3', 'init.A']:
        assert n in names
    should_values = {'V3':(0.1, 1.0), 'Km3':(0.0,5.0), 'init.A':(1.0,3.0)}
    for b in uu:
        assert_almost_equal(b.bounds.lower, should_values[b.name][0])
        assert_almost_equal(b.bounds.upper, should_values[b.name][1])

def test_reassignment2():
    """test reassignment of reactions"""
    m = Model("My first model")
    m.set_reaction('v1', "A -> B"  , 4)
    m.set_reaction('v2', "B -> C"  , 2.0)
    assert isinstance(m.reactions.v1, model.Reaction)
    assert isinstance(m.reactions.v2, model.Reaction)
    assert (m.reactions.v1.name) == 'v1'
    assert (m.reactions.v2.name) == 'v2'
    assert m.reactions.v1()== str(float(4))+ "*A"
    assert m.reactions.v2()== str(float(2.0))+"*B"
    check, msg = m.checkRates()
    assert check 
    m.set_reaction('v2',"D->C", 2.0)
    assert isinstance(m.reactions.v1, model.Reaction)
    assert isinstance(m.reactions.v2, model.Reaction)
    assert (m.reactions.v1.name) == 'v1'
    assert (m.reactions.v2.name) == 'v2'
    assert m.reactions.v1()== str(float(4))+ "*A"
    assert m.reactions.v2()== str(float(2.0))+"*D"
    check, msg = m.checkRates()
    assert check 

def test_reassignment3():
    """test change of variables by reassignment of reactions"""
    m = Model("My first model")
    m.set_reaction('v1', "A -> B"  , 4)
    m.set_reaction('v2', "B -> C"  , 2.0)
    assert len(m.varnames) == 3
    assert m.varnames == ['A', 'B', 'C']
    check, msg = m.checkRates()
    assert check 
    m.set_reaction('v2',"B->D", 2.0)
    xx = m.varnames
    assert len(xx) == 3
    assert xx == ['A', 'B', 'D']
    check, msg = m.checkRates()
    assert check 

def test_illegal_type1():
    """test illegal type assignment"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", 4)
    m.set_reaction('v2', "B->C", 2.0)
    with pytest.raises(model.BadTypeComponent):
        m.parameters.Km = [9,10,13,45]

def test_meta1():
    """test Model metadata"""
    m = Model("My first model")
    m.set_reaction('v1', "A->B", 4)
    m.set_reaction('v2', "B->C", 2.0)
    check, msg = m.checkRates()
    assert check 
    m.metadata['where'] = 'in model'
    m.metadata['for what'] = 'testing'
    assert m.metadata['where'] == 'in model'
    assert m.metadata['for what'] == 'testing'
    assert m.metadata['title'] == 'My first model'
    assert m.metadata.get('nonexistent', None) is None
    del m.metadata['where']
    assert m.metadata.get('where', None) is None

if __name__ == '__main__':
    pytest.main()