from stimator import Model, react, state, read_model
from stimator.GDE3solver import GDE3Solver
from time import time
from stimator import utils


def compute():

    npoints = 240
    t0 = 0.0
    tf = 120
    
    ## m1 = read_model("""
## title model 1
## rf: mgo + gsh -> hta, 0.34..
## rr: hta -> mgo + gsh, 1.01..
## r1: hta -> sdlt,      kcat1 * e1 * hta  / (km1 + hta)
## r2: sdlt -> gsh,      kcat2 * e2 * sdlt / (km2 + sdlt)
## fake1: e1->,          0 ..
## fake2: e2 ->,         0 ..
## kcat1 = 8586
## km1 = 0.223
## kcat2 = 315
## km2 = 2.86
## init = state(mgo = 2.86, hta = 0, sdlt = 0, gsh = 4, e1 = 2e-3, e2 = 4e-4)
## """)
    m1 = Model('model 1')
    m1.rf    = react("mgo + gsh -> hta", 0.34)
    m1.rr    = react("hta -> mgo + gsh", 1.01)
    m1.r1    = react("hta -> sdlt",      "kcat1 * e1 * hta / (km1 + hta)")
    m1.r2    = react("sdlt -> gsh",      "kcat2 * e2 * sdlt / (km2 + sdlt)")
    m1.fake1 = react("e1 ->", "0")
    m1.fake2 = react("e2 ->", "0")
    m1.kcat1 = 8586
    m1.km1   = 0.223
    m1.kcat2 = 315
    m1.km2   = 2.86
    m1.init  = state(mgo  = 2.86, 
                     hta  = 0, 
                     sdlt = 0, 
                     gsh  = 4, 
                     e1   = 2e-3, 
                     e2   = 4e-4)

    m2 = m1.clone()
    m2['title'] = 'model 2'
    m2.r1 = react("mgo + gsh -> sdlt"  , "kcat1 *e1 * mgo * gsh / ((km11 + gsh)*(km12 + mgo))")
    m2.kcat1 = 17046
    m2.km11  = 0.875
    m2.km12  = 1.178
    
    models   = [m1, m2]
    
    toOpt    = {"mgo":[0.1, 1], "gsh":[0.1, 1]}#, "e1":[1.9e-3, 2.0e-3], "e2":[3.9e-4, 4.0e-4]}
    
    observed = 'sdlt'
    
    objectiveFunction = 'KL'
    populationSize = 200
    maxGenerations = 400
    DEStrategy = 'Rand1Bin'
    diffScale = 0.5
    crossoverProb = 0.7
    cutoffEnergy = 0 #For now. Adjust for each objective function.
    useClassRandomNumberMethods = True
    simulatedError = 3
    absoluteMeasurementError = 0.00175

    print 'initial values to optimize:'
    for k in toOpt.keys():
        print k, 'in', toOpt[k]
    
    allTime1 = time()
    solver = GDE3Solver(models, 
                       toOpt, 
                       objectiveFunction, 
                       observed, 
                       npoints, t0, tf, 
                       populationSize, maxGenerations, 
                       DEStrategy, 
                       diffScale, 
                       crossoverProb, 
                       cutoffEnergy, 
                       useClassRandomNumberMethods)#, dif = '-')
    solver.Solve()
    finalSolutions = (solver.population, solver.population_energies)
    
    print '============================================='
    print "Finished!"
    ttime = time() - allTime1
    print "Total time: %g s (%s)"% (ttime, utils.s2HMS(ttime))
    
    print '%d generations'%(solver.generation-1)
    print
    print 'Final front:'
    
    f = open ('glyoxalase_discrim_2m_KLext_final_front.txt', 'w')
    
    sstr = ' '.join(solver.toOptKeys)
    ostr = ' '.join([str(d) for d in solver.model_indexes])
    print '[%s] ----> [%s]' % (sstr, ostr)
    for s,o in zip(*finalSolutions):
        print s, '---->',o
        sstr = ' '.join([str(i) for i in s])
        ostr = ' '.join([str(i) for i in o])
        print >> f, '%s %s'%(sstr, ostr)
    
    f.close()

    utils.write2file('glyoxalase_discrim_2m_KLext_times.txt', '\n'.join([str(t) for t in solver.gen_times]))
    
if __name__ == "__main__":
    compute()
