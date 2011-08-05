from model import *
from analysis import *
from GDE3solver import *
from util import write2file

def compute():

    npoints = 240
    t0 = 0.0
    tf = 120
    
    #TODO: set estimated parameter values
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
    m1.init = state(mgo = 2.86, hta = 0, sdlt = 0, gsh = 4, e1 = 0, e2 = 0) #init state must however be defined here so that there are values to replace by the function generateRSS
    #m1.lenEstimatedParameters = len(("v1", "km1"))

    m2 = m1.clone()
    m2.title = 'model 2'
    m2.r1 = react("mgo + gsh -> sdlt"  , "kcat1 *e1 * mgo * gsh / ((km11 + gsh)*(km12 + mgo))")
    m2.kcat1 = 17046
    m2.km11 = 0.875
    m2.km12 = 1.178
    #m2.v2 = 0.658
    #m2.km2 = 0.0347
    #m2.init = state(mgo = 2.86, hta = 0, sdlt = 0, gsh = 4)
    #m2.lenEstimatedParameters = len(("v1", "km11", "km12"))
    
    #~ m3 = m1.clone()
    #~ m3.title = 'model 3'
    #~ m3.r1 = react("hta -> sdlt"  , "kcat1 * hta / (km1 * (1 + gsh/Ki) + hta)")
    #~ #m3.r12= react("mgo + gsh -> sdlt"  , "v1b * mgo * gsh / ((km11 + gsh)*(km12 + mgo))")
    #~ #m3.v1a = 0.164
    #~ m3.km1 = 0.208
    #~ m3.Ki = 4.106
    #~ #m3.v1b = 0.025
    #~ #m3.km11 = 0.481
    #~ #m3.km12 = 0.301
    #~ #m3.init = state(mgo = 2.86, hta = 0, sdlt = 0, gsh =4)
    #~ #m3.lenEstimatedParameters = len(("v1a", "v1b", "km1", "km11", "km12"))
    
    models = [m1, m2]#, m3]
    
    toOpt = {"mgo":[0.1, 1], "gsh":[0.1, 1], "e1":[0, 2*10**-3], "e2":[0, 4*10**-4]}#, "hta":[0, 1.2]}
    
    observed = 'sdlt'
    
    biasedCurveNumber = 3
    
    biasStandardDeviation = 0.03

    #TODO: how to set the energy functions to be used in the optimization?
    objectiveFunction = 'L2'#, 'AICm3') #Other examples: ('AICm1', 'AICm2')... and then the adequate function is selected among the ones defined in the class.
    populationSize = 200
    maxGenerations = 5000
    DEStrategy = 'Rand1Bin'
    diffScale = 0.5
    crossoverProb = 0.7
    cutoffEnergy = 0 #For now. Adjust for each objective function.
    useClassRandomNumberMethods = True
    simulatedError = 3
    absoluteMeasurementError = 0.00175
    
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
    allSolutions = (solver.completeListOfSolutions,
                    solver.completeListOfObjectives)
    finalSolutions = (solver.fronts, solver.frontObj)
    
    print 
    print "Finished!"
    print "Total time", time() - allTime1
    
##     print '\n\nfinalSolutions\n', finalSolutions
##     print "%d fronts"%(len(finalSolutions[0]))
    print '%d generations'%(solver.generation-1)
    print 'front lengths:'
    print [len(i) for i in solver.completeListOfSolutions]
    print
    print 'Final front:'
    
    f = open ('results/final_front_L2.txt', 'w')
    
    for s,o in zip(finalSolutions[0][-1], finalSolutions[1][-1]):
        print s, '---->',o
        sstr = ' '.join([str(i) for i in s])
        ostr = ' '.join([str(i) for i in o])
        f.write('%s %s\n'%(sstr, ostr))
    
    f.close()
    
    write2file('results/times.txt', '\n'.join([str(t) for t in solver.ftimes]))
    
##     aString = 'Asolutions = {'
##     for i in allSolutions[0]:
##         aString += '{'
##         for j in i:
##             aString += str(j) + ','
##         aString = aString[:-1] + '},'
##     aString = aString[:-1] + '}'
##     write2file(r'results/candSolsKrem2M.txt', aString)

##     bString = '{'
##     for i in allSolutions[1]:
##         bString += '{'
##         for j in i:
##             bString += str(j) + ','
##         bString = bString[:-1] + '},'
##     bString = bString[:-1] + '}'
##     write2file(r'results/candObjsKrem2M.txt', bString)
##     
##     fString = 'solutions = {'
##     for i in finalSolutions[0]:
##         fString += '{'
##         for j in i:
##             fString += str(j) + ','
##         fString = fString[:-1] + '},'
##     fString = fString[:-1] + '}'
##     write2file(r'results/finalSolsKrem2M.txt', fString)

##     gString = '{'
##     for i in finalSolutions[1]:
##         gString += '{'
##         for j in i:
##             gString += str(j) + ','
##         gString = gString[:-1] + '},'
##     gString = gString[:-1] + '}'
##     write2file(r'results/finalObjsKrem2M.txt', gString)
    
if __name__ == "__main__":
    compute()
