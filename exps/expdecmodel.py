""" second experiment with function introspection and function decoration """

import sys
import os.path

#append parent directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

import inspect
import tokenize
import StringIO
from decmodel import model
from stimator.analysis import *


#print '---------------- EXAMPLE: CICR model ------------------'

def step (t, t_stimulus, B):
    if t < t_stimulus:
        return 0.0
    else:
        return B

@model
def cicr(title = "Calcium Spikes"):
    v0         = ' -> Ca', 1
    v1         = " -> Ca", \
                 Bstep*k1
    
    Bstep      = 0.4
    #Bstep      = forcing(step)
    k1         = 7.3
    B          = 0.4
    t_stimulus = 1.0
    
    export     = "Ca ->       ", 10
    leak       = "CaComp -> Ca", 1
    
    v2         = "Ca -> CaComp", \
                  65 * Ca**2 / (1+Ca**2)
    
    v3         = "CaComp -> Ca", \
                  500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)
    
    init       = state(Ca = 0.1, CaComp = 0.63655)

cicr.title= "Calcium Spikes"
print cicr

svars = solve(cicr, tf = 6.0, npoints = 1000, title = 'X')
sdxdt = solve(cicr, tf = 6.0, npoints = 5000, title = 'dX/dt') #TODO: solutions should be 'clonable'.

transformation = cicr.getdXdt()
sdxdt.apply_transf(transformation)

# print '------------ LAST POINT: -------------'
# print sdxdt.last

plot([svars,sdxdt])
