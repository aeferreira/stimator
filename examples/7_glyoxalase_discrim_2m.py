from __future__ import print_function
from stimator import Model, read_model
from stimator.GDE3solver import GDE3Solver
from time import time

def compute(obj_func):

    npoints = 240
    t0 = 0.0
    tf = 120
    m1 = read_model("""title model 1
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
init  = state(mgo  = 2.86, hta = 0, sdlt = 0, gsh  = 4, e1 = 2e-3, e2   = 4e-4)
""")

    m2 = m1.copy(new_title = 'model 2')
    m2.set_reaction('r1', "mgo + gsh -> sdlt"  , "kcat1 *e1 * mgo * gsh / ((km11 + gsh)*(km12 + mgo))")
    m2.parameters.kcat1 = 17046
    m2.parameters.km11  = 0.875
    m2.parameters.km12  = 1.178
    
    models   = [m1, m2]
    
    initial_opt = (('gsh', 0.1, 1.0), ('mgo', 0.1, 1.0))
    observed    = ['sdlt']
    
    ## oOpt = {"mgo":[0.1, 1], "gsh":[0.1, 1]}, 
    ##         "e1":[1.9e-3, 2.0e-3], "e2":[3.9e-4, 4.0e-4]}
    
    objectiveFunction = obj_func
    populationSize = 200
    maxGenerations = 100
    DEStrategy = 'Rand1Bin'
    diffScale = 0.5
    crossoverProb = 0.7
    cutoffEnergy = 0 #Not used in multiobjective optimization
    useClassRandomNumberMethods = True
    dump_generations = None # do not generate generation log file
    #dump_generations = list(range(maxGenerations)) #generate log for all generations

    print ('==========================================================')
    print ('Design of discriminatory experiment (initial conditions)\n')
    print ('Metrics used: %s\n' % objectiveFunction)
    print ('observed variables: %s\n' % observed)
    print ('initial values to optimize:')
    for name, min_v, max_v in initial_opt:
        print (name, 'in [%g, %g]' %(min_v, max_v))
    
    solver = GDE3Solver(models, 
                       initial_opt, 
                       objectiveFunction, 
                       observed, 
                       npoints, t0, tf, 
                       populationSize, maxGenerations, 
                       DEStrategy, 
                       diffScale, 
                       crossoverProb, 
                       cutoffEnergy, 
                       useClassRandomNumberMethods,
                       dump_generations = dump_generations)#, dif = '-')
    solver.run()
    
    print ('-----------------------------------------------')
    print ('Final front:')
    f = open ('glyoxalase_discrim_2m_%s_final_pop.txt'%objectiveFunction, 'w')
    
    sstr = ' '.join(solver.toOptKeys)
    ostr = ' '.join([str(d) for d in solver.model_indexes])
    print ('%s ----> %s' % (sstr, ostr))
    for s,o in zip(solver.population, solver.population_energies):
        sstr = ' '.join(["%6.4g"%i for i in s])
        ostr = ' '.join(["%6.4g"%i for i in o])
        print('%s ----> %s'%(sstr, ostr))
        print('%s %s'%(sstr, ostr), file=f)
    f.close()

if __name__ == "__main__":
    obj_funcs = ['extKL']
    #obj_funcs = ['extKL', 'KL', 'L2', 'L2_midpoint_weights']
    for obj in obj_funcs:
        compute(obj)
