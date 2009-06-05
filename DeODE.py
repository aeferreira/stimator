#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2009 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from time import time
import DE
from numpy import *
from scipy import integrate
from model import *
import timecourse

DUMP_PARS_2FILES = False

#----------------------------------------------------------------------------
#         Class to perform DE optimization for ODE systems
#----------------------------------------------------------------------------

class DeODESolver(DE.DESolver):
    """Overides energy function and report functions.
    
    The energy function solves ODEs and computes a least-squares score.
    Ticker functions are called on generation completion and when optimization finishes.
    """
    
    def __init__(self, model, optSettings, timecoursecollection, aGenerationTicker, anEndComputationTicker):
        self.model = model
        self.tc    = timecoursecollection
        self.timecoursedata   = self.tc.data
        self.generationTicker = aGenerationTicker
        self.endTicker        = anEndComputationTicker

        pars = model.uncertain
        mins = array([u.min for u in pars])
        maxs = array([u.max for u in pars])
        
        DE.DESolver.__init__(self, len(pars), # number of parameters
                             int(optSettings['genomesize']),  # genome size
                             int(optSettings['generations']), # max number of generations
                             mins, maxs,              # min and max parameter values
                             "Best2Exp",              # DE strategy
                             0.7, 0.6, 0.0,           # DiffScale, Crossover Prob, Cut off Energy
                             True)                    # use class random number methods

        # cutoffEnergy is 0.1% of deviation from data
        self.cutoffEnergy =  1.0e-6*sum([nansum(abs(tc[:,1:])) for tc in self.timecoursedata])
        
        # scale times to maximum time in data
        scale = float(max([ (tc[-1,0]-tc[0,0]) for tc in self.timecoursedata]))
        
        self.calcDerivs = model.scaled_dXdt(scale, with_uncertain=True)
        
        # store initial values and (scaled) time points
        self.X0 = []
        self.times = []
        for data in self.timecoursedata:
            y0 = copy(data[0, 1:]) # variables are in columns 1 to end
            self.X0.append(y0)
            
            t  = data[:, 0]        # times are in columns 0
            t0 = t[0]
            times = (t-t0)/scale+t0  # this scales time points
            self.times.append(times)

        self.timecourse_scores = empty(len(self.timecoursedata))
        
        # open files to write parameter progression
        if DUMP_PARS_2FILES:
            self.parfilehandes = [open(par[0]+".par", 'w') for par in model.parameters]

    def externalEnergyFunction(self,trial):
        for par in range(self.parameterCount):
            if trial[par] > self.maxInitialValue[par] or trial[par] < self.minInitialValue[par]:
                return 1.0E300
        
        self.model.set_uncertain(trial)
        salg=integrate._odepack.odeint

        for i in range(len(self.timecoursedata)):
            y0 = copy(self.X0[i])
            t  = copy(self.times[i])

#           Y, infodict = integrate.odeint(self.calcDerivs, y0, t, full_output=True, printmessg=False)
            output = salg(self.calcDerivs, y0, t, (), None, 0, -1, -1, 0, None, 
                            None, None, 0.0, 0.0, 0.0, 0, 0, 0, 12, 5)
            if output[-1] < 0: return (1.0E300)
            #~ if infodict['message'] != 'Integration successful.':
                #~ return (1.0E300)
            Y = output[0]
            S = (Y- self.timecoursedata[i][:, 1:])**2
            score = nansum(S)
            self.timecourse_scores[i]=score
        
        gscore = self.timecourse_scores.sum()
        return gscore

    def reportGeneration (self):
        if not self.generationTicker:
            print self.reportGenerationString()
        else:
            self.generationTicker(self.generation, float(self.bestEnergy))
        if DUMP_PARS_2FILES:
            for par in range(self.parameterCount):
                parvector = [str(self.population[k][par]) for k in range(self.populationSize)]
                print >>self.parfilehandes[par], " ".join(parvector)
            
    
    def reportFinal (self):
        if self.exitCode==0: outCode = -1 
        else: 
            outCode = self.exitCode
            self.optimum = self.generateOptimumData()
            #reportstring = self.reportFinalString()
        if not self.endTicker:
            print self.reportFinalString()
        else:
            self.endTicker(outCode)
        if DUMP_PARS_2FILES:
            for par in self.parfilehandes:
                par.close()

    def generateOptimumData (self):
        timecoursedata = self.timecoursedata
        bestData = {'parameters'       : {'name':"parameters"}, 
                    'optimization'     : {'name':"D.E. optimization"}, 
                    'timecourses'      : {'name':"timecourses"},
                    'best timecourses' : {'name':"best timecourses"}}
        bestData['parameters']['data'] = [(self.model.uncertain[i].name, "%g"%value) for (i,value) in enumerate(self.bestSolution)]
        bestData['parameters']['format'] = "%s\t%s"
        bestData['parameters']['header'] = None
        
        bestData['optimization']['data'] = [('Final Score', "%g"% self.bestEnergy),
                                ('Generations', "%d"% self.generation),
                                ('Exit by    ', "%s"% self.exitCodeStrings[self.exitCode])]
        bestData['optimization']['format'] = "%s\t%s"
        bestData['optimization']['header'] = None

        #TODO: Store initial solver parameters?

        #generate best time-courses
        bestData['timecourses']['data'] = []
        bestData['best timecourses']['data'] = []
        besttimecoursedata = bestData['best timecourses']['data']
        for (i,data) in enumerate(timecoursedata):
            y0 = copy(self.X0[i])
            t  = self.times[i]
            Y, infodict = integrate.odeint(self.calcDerivs, y0, t, full_output=True, printmessg=False)
            besttimecoursedata.append(Y)
            if infodict['message'] != 'Integration successful.':
                bestData['timecourses']['data'].append((self.tc.shortnames[i], self.tc.shapes[i][0], 1.0E300))
            else:
                S = (Y- data[:, 1:])**2
                score = nansum(S)
                bestData['timecourses']['data'].append((self.tc.shortnames[i], self.tc.shapes[i][0], score))
        bestData['timecourses']['format'] = "%s\t%d\t%g"
        bestData['timecourses']['header'] = ['Name', 'Points', 'Score']
        
        return bestData

def test():
    m1 = Model("Glyoxalase system in L.infantum")
    m1.glo1 = react("HTA -> SDLTSH", rate = "V1*HTA/(Km1 + HTA)")
    m1.glo2 = react("SDLTSH -> "   , rate = "V2*SDLTSH/(Km2 + SDLTSH)")
    m1.V1  = 2.57594e-05
    m1.V1.uncertainty(0.00001, 0.0001)
    m1.Km1 = 0.252531
    m1.Km1.uncertainty(0.01, 1)
    m1.V2  = 2.23416e-05
    m1.V2.uncertainty(0.00001, 0.0001)
    m1.Km2 = 0.0980973
    m1.Km2.uncertainty(0.01, 1)
    m1.init = state(SDLTSH = 7.69231E-05, HTA = 0.1357)
    #print m1
    
    optSettings={'genomesize':80, 'generations':200}
    timecourses = timecourse.readTimeCourses(['TSH2a.txt', 'TSH2b.txt'], 'examples', (0,2,1))
    
    solver = DeODESolver(m1,optSettings, timecourses, None, None)
    
    time0 = time()
    
    solver.Solve()
    
    print "Optimization took %f s"% (time()-time0)


    print
    print '---------------------------------------------------------'
    print "Results for %s" % m1.title

    # generate report
    reportText = ""
    sections = [solver.optimum[s] for s in ['parameters', 'optimization', 'timecourses']]
    for section in sections:
        reportText += "%-20s --------------------------------\n" % section['name'].upper()
        if section['header']:
            reportText += '\t'.join(section['header'])+'\n'
        reportText += "\n".join([section['format'] % i for i in section['data']])
        reportText += '\n\n'
    print reportText


if __name__ == "__main__":
    test()
