from stimator import Model, react, state
from stimator.GDE3solver import GDE3Solver
from time import time
from stimator.utils import write2file

def compute():

    npoints = 240
    t0 = 0.0
    tf = 120
    
    m1 = Model('model 1')
    m1.rf = react("mgo + gsh -> hta", 0.34)
    m1.rr = react("hta -> mgo + gsh", 1.01)
    m1.r1 = react("hta -> sdlt"  , "kcat1 * e1 * hta / (km1 + hta)")
    m1.r2 = react("sdlt -> gsh", "kcat2 * e2 * sdlt / (km2 + sdlt)")
    m1.fake1 = react("e1 ->", "0")
    m1.fake2 = react("e2 ->", "0")
    m1.kcat1 = 8586
    m1.km1 = 0.223
    m1.kcat2 = 315
    m1.km2 = 2.86
    m1.init = state(mgo = 2.86, hta = 0, sdlt = 0, gsh = 4, e1 = 0, e2 = 0)

    m2 = m1.clone()
    m2[title] = 'model 2'
    m2.r1 = react("mgo + gsh -> sdlt"  , "kcat1 *e1 * mgo * gsh / ((km11 + gsh)*(km12 + mgo))")
    m2.kcat1 = 17046
    m2.km11 = 0.875
    m2.km12 = 1.178
    
    models = [m1, m2]
    
    toOpt = {"mgo":[0.1, 1], "gsh":[0.1, 1], "e1":[0, 2*10**-3], "e2":[0, 4*10**-4]}#, "hta":[0, 1.2]}
    
    observed = 'sdlt'
    
    biasedCurveNumber = 3
    
    biasStandardDeviation = 0.03

    objectiveFunction = 'kremling'
    populationSize = 100
    maxGenerations = 200
    DEStrategy = 'Rand1Bin'
    diffScale = 0.5
    crossoverProb = 0.7
    cutoffEnergy = 0 #For now. Adjust for each objective function.
    useClassRandomNumberMethods = True
    simulatedError = 3
    absoluteMeasurementError = 0.00175
    allTime1 = time()
    
    solver = GDE3Solver(models, 
                           absoluteMeasurementError, 
                           toOpt, 
                           objectiveFunction, 
                           observed, npoints, t0, tf, 
                           populationSize, maxGenerations, 
                           biasedCurveNumber, 
                           biasStandardDeviation, 
                           DEStrategy, 
                           diffScale, 
                           crossoverProb, 
                           cutoffEnergy, 
                           useClassRandomNumberMethods)#, dif = '-')
    solver.Solve()
##     allSolutions = (solver.completeListOfSolutions,
##                     solver.completeListOfObjectives)
    finalSolutions = (solver.fronts, solver.frontObj)
    
    print 
    print "Finished!"
    print "Total time", time() - allTime1
    
##     print '\n\nfinalSolutions\n', finalSolutions
##     print "%d fronts"%(len(finalSolutions[0]))
    print '%d generations'%(solver.generation-1)
##     print 'front lengths:'
##     print [len(i) for i in solver.completeListOfSolutions]
##     print
    print 'Final front:'
    
    f = open ('glyoxalase_discrim_2m_L2w_final_front.txt', 'w')
    
    for s,o in zip(finalSolutions[0][-1], finalSolutions[1][-1]):
        print s, '---->',o
        sstr = ' '.join([str(i) for i in s])
        ostr = ' '.join([str(i) for i in o])
        f.write('%s %s\n'%(sstr, ostr))
    
    f.close()

    write2file('glyoxalase_discrim_2m_L2w_times.txt', '\n'.join([str(t) for t in solver.ftimes]))

if __name__ == "__main__":
    compute()
