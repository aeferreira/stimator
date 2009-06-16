import sys
import os.path

#append parent directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir))

from stimator import *

m = Model("My first model")
m.v1 = react("A->B", 4)
m.v2 = react("B->C", 2.0)
print m.v1.rate
print m.v2.rate
