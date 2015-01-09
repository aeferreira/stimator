# -*- coding:utf8 -*-
class ExampleModel(object):
    pass

lorentz = ExampleModel()
lorentz.desc = "Lorentz model. Deterministic chaos"
lorentz.text = """title Lorentz model: sensitivity to initial conditions
x' = 10*(y-x)
y' = x*(28-z)-y
z' = x*y - (8/3)*z
init: x = 1, y = 1, z = 1
"""

rossler = ExampleModel()
rossler.desc = "Rossler model. Deterministic chaos"
rossler.text = """title Rossler

X1' = X2 - X3
X2' = 0.36 * X2 - X1
X3' = X1 *X3 - 22.5 * X3 - 49.6 * X1 + 1117.8

init : (X1 = 19.0, X2 = 47, X3 = 50)
tf:200
~ x3 = X3 -50.0
~ x1 = X1 -18.0
~ x2 = X2 -50.0
"""

ca = ExampleModel()
ca.desc = "CICR model. Calcium spikes"
ca.text = """title CICR model. Calcium spikes
v0         = -> Ca, 1
v1         = -> Ca, k1*B*step(t, t_stimulus)

k1         = 7.3
B          = 0.4
t_stimulus = 1.0
    
export     = Ca ->  , 10 ..
leak       = CaComp -> Ca, 1 ..
    
v2         = Ca -> CaComp, 65 * Ca**2 / (1+Ca**2)
v3         = CaComp -> Ca, 500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)
init       : Ca = 0.1, CaComp = 0.63655
"""

branched = ExampleModel()
branched.desc = "Branched pathway"
branched.text = """title Branched pathway

v1 = A -> B, k1*A,      k1 = 10
v2 = B -> C, k2*B**0.5, k2 = 5
v3 = C -> D, k3*C**0.5, k3 = 2
v4 = C -> E, k4*C**0.5, k4 = 8
v5 = D ->  , k5*D**0.5, k5 = 1.25
v6 = E ->  , k6*E**0.5, k6 = 5
A  = 0.5
init : (B = 2, C = 0.25, D = 0.64, E = 0.64)
!! B D
tf: 10
"""

glos_hta = ExampleModel()
glos_hta.desc = "Glyoxalase system in L. Infantum"
glos_hta.text = """title Glyoxalase system in L. Infantum
variables SDLTSH HTA  # variables (the order matches the timecourse files)

#reactions (with stoichiometry and rate)
glx1 : HTA -> SDLTSH, rate = V1*HTA/(Km1 + HTA), V1 = 0.00005
glx2 : SDLTSH ->, V2*SDLTSH/(Km2 + SDLTSH)

find glx1.V1  in [0.00001, 0.0001]  # parameters to find and bounding intervals
find Km1 in [0.01, 1]
find V2  in [0.00001, 0.0001]
find Km2 in [0.01, 1]
pi = 3.14159
init : (SDLTSH = 7.69231E-05, HTA = 0.1357)

timecourse TSH2a.txt  # timecourses to load
timecourse TSH2b.txt
tf:4500
generations = 200   # maximum generations for GA
genomesize = 80     # population size in GA
"""

glyoxalases = ExampleModel()
glyoxalases.desc = "Glyoxalase system in L. Infantum"
glyoxalases.text = """title Glyoxalase system
glo1 = HTA -> SDLTSH  , rate = V1*HTA/(Km1 + HTA)
glo2 = SDLTSH -> DLac , rate = V2*SDLTSH/(Km2 + SDLTSH)
V1  = 2.57594e-05
Km1 = 0.252531
V2  = 2.23416e-05
Km2 = 0.0980973
init :(SDLTSH = 7.69231E-05, HTA = 0.1357)
tf: 4030
"""
