import sys
import os.path

#append parent directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir))

from stimator import *
from stimator import model


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
    assert isinstance(m, Model)
    assert isinstance(m.v1, model.Reaction)
    assert isinstance(m.v2, model.Reaction)
    assert m.v1.rate== str(float(4))+ "*A"
    assert m.v2.rate== str(float(2.0))+"*B"

def test_react2():
    """test react(string, string)"""
    m = Model("My first model")
    m.v1 = react("A->B", "4*A/(2+A)")
    assert isinstance(m, Model)
    assert isinstance(m.v1, model.Reaction)
    assert m.v1.rate== "4*A/(2+A)"
