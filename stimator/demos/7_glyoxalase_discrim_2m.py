from __future__ import print_function
from stimator import read_model
from stimator.GDE3solver import GDE3Solver

model1_desc = """title model 1
rf: mgo + gsh -> hta, 0.34..
rr: hta -> mgo + gsh, 1.01..
r1: hta -> sdlt, kcat1 * e1 * hta / (km1 + hta)
r2: sdlt -> gsh, kcat2 * e2 * sdlt / (km2 + sdlt)
fake1: e1 ->, 0
fake2: e2 ->, 0

kcat1 = 8586
km1   = 0.223

kcat2 = 315
km2   = 2.86

init: mgo = 2.86, hta = 0, sdlt = 0, gsh = 4, e1 = 2e-3, e2 = 4e-4
"""

model2_desc = """title model 2
rf: mgo + gsh -> hta, 0.34..
rr: hta -> mgo + gsh, 1.01..
r1: mgo + gsh -> sdlt, kcat1 *e1 * mgo * gsh / ((km11 + gsh)*(km12 + mgo))
r2: sdlt -> gsh, kcat2 * e2 * sdlt / (km2 + sdlt)
fake1: e1 ->, 0
fake2: e2 ->, 0

kcat1 = 17046

km11  = 0.875
km12  = 1.178

kcat2 = 315
km2   = 2.86

init: mgo = 2.86, hta = 0, sdlt = 0, gsh = 4, e1 = 2e-3, e2 = 4e-4
"""

def compute(obj):

    npoints = 240
    t0 = 0.0
    tf = 120
    models = [read_model(model1_desc), read_model(model2_desc)]
    
    initial_opt = (('gsh', 0.1, 1.0), ('mgo', 0.1, 1.0))
    observed = ['sdlt']
    
    ## oOpt = {"mgo":[0.1, 1], "gsh":[0.1, 1]}, 
    ##         "e1":[1.9e-3, 2.0e-3], "e2":[3.9e-4, 4.0e-4]}
    
    populationSize = 200
    max_generations = 100
    DEStrategy = 'Rand1Bin'
    diffScale = 0.5
    crossoverProb = 0.7
    cutoffEnergy = 0 #Not used in multiobjective optimization
    useClassRandomNumberMethods = True
    dump_generations = None # do not generate generation log file
    # generate log for all generations
    #dump_generations = list(range(maxGenerations))

    print('==========================================================')
    print('Design of discriminatory experiment (initial conditions)\n')
    print ('Metrics used: {}\n'.format(obj))
    print ('observed variables: {}\n'.format(observed))
    print ('initial values to optimize:')
    for name, min_v, max_v in initial_opt:
        print (name, 'in [%g, %g]' % (min_v, max_v))
    
    solver = GDE3Solver(models, 
                       initial_opt, 
                       obj, 
                       observed, 
                       npoints, t0, tf, 
                       populationSize,
                       DEStrategy, 
                       diffScale, 
                       crossoverProb, 
                       cutoffEnergy, 
                       useClassRandomNumberMethods,
                       dump_generations=dump_generations,
                       max_generations=max_generations)#, dif = '-')
    solver.run()
    
    print ('-----------------------------------------------\nFinal front:')
    with open('glos_2m_{}_final_pop.txt'.format(obj), 'w') as f:
        sstr = ' '.join(solver.toOptKeys)
        ostr = ' '.join([str(d) for d in solver.model_indexes])
        print ('{} ----> {}'.format(sstr, ostr))
        for s,o in zip(solver.pop, solver.population_energies):
            sstr = ' '.join(["%6.4g"%i for i in s])
            ostr = ' '.join(["%6.4g"%i for i in o])
            print('{} ----> {}'.format(sstr, ostr))
            print('{} {}'.format(sstr, ostr), file=f)

if __name__ == "__main__":
    obj_funcs = ['extKL']
    #obj_funcs = ['extKL', 'KL', 'L2', 'L2_midpoint_weights']
    for obj in obj_funcs:
        compute(obj)
