#!/usr/bin/env python
# -*- coding:latin1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2009 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

from time import time
import de
from numpy import *
from scipy import integrate
from model import *
from modelparser import read_model
from analysis import *
import fim
import timecourse
import dyncriteria

#----------------------------------------------------------------------------
#         Class to perform DE optimization for ODE systems
#----------------------------------------------------------------------------

class DeODESolver(de.DESolver):
    """Overides energy function and report functions.
    
    The energy function solves ODEs and computes a least-squares score.
    Ticker functions are called on generation completion and when optimization finishes.
    """
    
    def __init__(self, model, optSettings, tcs, weights = None,
                    aGenerationTicker=None, anEndComputationTicker=None, 
                    dump_pars=False):
        self.model = model
        self.tc    = tcs
        self.generationTicker = aGenerationTicker
        self.endTicker        = anEndComputationTicker
        self.dump_pars        = dump_pars

        pars = model.uncertain
        mins = array([u.min for u in pars])
        maxs = array([u.max for u in pars])
        
        de.DESolver.__init__(self, len(pars), # number of parameters
                             int(optSettings['genomesize']),  # genome size
                             int(optSettings['generations']), # max number of generations
                             mins, maxs,              # min and max parameter values
                             "Best2Exp",              # DE strategy
                             0.7, 0.6, 0.0,           # DiffScale, Crossover Prob, Cut off Energy
                             True)                    # use class random number methods

        # cutoffEnergy is 1e-6 of deviation from data
        self.cutoffEnergy =  1.0e-6*sum([nansum(abs(tc.data)) for tc in self.tc])
        
        # scale times to maximum time in data
        scale = float(max([ (tc.t[-1]-tc.t[0]) for tc in self.tc]))
        
        self.calcDerivs = model.getdXdt(scale=scale, with_uncertain=True)
        self.salg=integrate._odepack.odeint
        
        # store initial values and (scaled) time points
        self.X0 = []
        self.times = []
        self.varindexes = []
        for data in self.tc:
            self.X0.append(copy(data[:, 0].T))
            t  = data.t
            t0 = t[0]
            times = (t-t0)/scale+t0  # this scales time points
            self.times.append(times)
        self.timecourse_scores = empty(len(self.tc))
        
        # find uncertain initial values
        self.mapinit2trial = []
        for iu, u in enumerate(self.model.uncertain):
            if u.name.startswith('init'):
                varname = u.name.split('.')[-1]
                ix = findWithNameIndex(varname, self.model.variables)
                self.mapinit2trial.append((ix,iu))
        
        self.criterium = dyncriteria.getCriteriumFunction(weights, self.tc)

        # open files to write parameter progression
        if self.dump_pars:
            self.parfilehandes = [open(par[0]+".par", 'w') for par in model.parameters]

    def computeSolution(self,i,trial):
        """Computes solution for timecourse i, given parameters trial."""
        
        y0 = copy(self.X0[i])
        # fill uncertain initial values
        for varindex, trialindex in self.mapinit2trial:
            y0[varindex] = trial[trialindex]
            
        t  = copy(self.times[i])
        #~ Y, infodict = integrate.odeint(self.calcDerivs, y0, t, full_output=True, printmessg=False)
        output = self.salg(self.calcDerivs, y0, t, (), None, 0, -1, -1, 0, None, 
                        None, None, 0.0, 0.0, 0.0, 0, 0, 0, 12, 5)
        if output[-1] < 0: return None
        #~ if infodict['message'] != 'Integration successful.':
            #~ return (1.0E300)
        return output[0]
        

    def externalEnergyFunction(self,trial):
        #if out of bounds flag with error energy
        for trialpar, minInitialValue, maxInitialValue in zip(trial, self.minInitialValue, self.maxInitialValue):
            if trialpar > maxInitialValue or trialpar < minInitialValue:
                return 1.0E300
        #set parameter values from trial
        self.model.set_uncertain(trial)
        
        #compute solutions and scores
        for i in range(len(self.tc)):
            Y = self.computeSolution(i,trial)
            if Y is not None:
                self.timecourse_scores[i]=self.criterium(Y, i)
            else:
                return (1.0E300)
        
        globalscore = self.timecourse_scores.sum()
        return globalscore

    def reportGeneration (self):
        if not self.generationTicker:
            print self.reportGenerationString()
        else:
            self.generationTicker(self.generation, float(self.bestEnergy))
        if self.dump_pars:
            for par in range(self.parameterCount):
                parvector = [str(self.population[k][par]) for k in range(self.populationSize)]
                print >>self.parfilehandes[par], " ".join(parvector)
            
    
    def reportFinal (self):
        if self.exitCode==0: outCode = -1 
        else: 
            outCode = self.exitCode
            self.generateOptimumData()
        if not self.endTicker:
            print self.reportFinalString()
        else:
            self.endTicker(outCode)
        if self.dump_pars:
            for par in self.parfilehandes:
                par.close()

    def generateOptimumData (self):
        #compute parameter standard errors, based on FIM-1
        #generate TC solutions
        best = {'parameters'       : {'name':"parameters"},
                'optimization'     : {'name':"D.E. optimization"}, 
                'timecourses'      : {'name':"timecourses"},
                'best timecourses' : {'name':"best timecourses"}}
        best['parameters']['format'] = "%s\t%12s +- %s"
        best['parameters']['header'] = None
        
        best['optimization']['data'] = [('Final Score', "%g"% self.bestEnergy),
                                        ('Generations', "%d"% self.generation),
                                        ('Exit by    ', "%s"% self.exitCodeStrings[self.exitCode])]
        best['optimization']['format'] = "%s\t%s"
        best['optimization']['header'] = None

        #TODO: Store initial solver parameters?

        #generate best time-courses
        best['timecourses']['data'] = []
        #best['best timecourses']['data'] = []

        varnames = [x.name for x in self.model.variables]
        pars = [self.model.uncertain[i].name for i in range(len(self.bestSolution))]
        parvalues = [value for value in self.bestSolution]
        parsdict = dict (zip(pars, parvalues))
        self.model.set_uncertain(self.bestSolution)
        
        sols = timecourse.Solutions()
        
        for (i,tc) in enumerate(self.tc):
            Y = self.computeSolution(i, self.bestSolution)
            if Y is not None:
                score =self.criterium(Y, i)
            else:
                score = 1.0E300
            nt = tc.ntimes
            sol = timecourse.SolutionTimeCourse (tc.t, Y.T, varnames)
            sols += sol
            #best['best timecourses']['data'].append(sol)
            best['timecourses']['data'].append((self.tc.shortnames[i], self.tc[i].shape[0], score))
            
            varnames = []
            varindexes=[]
            for iline,line in enumerate(tc.data):
                #count NaN
                yexp = line
                nnan = len(yexp[isnan(yexp)])
                if nnan >= nt-1: continue
                varnames.append(str(self.model.variables[iline].name))
                varindexes.append(iline)
        best['timecourses']['format'] = "%s\t%d\t%g"
        best['timecourses']['header'] = ['Name', 'Points', 'Score']
        best['best timecourses']['data']=sols
        
        consterror = [0.0 for i in range(len(varnames))]
        for ix, x in enumerate(varnames):
            for tc in self.tc:
                yexp = tc.data[varindexes[ix]]
                tpe = (max(yexp) - min(yexp))
                if tpe > consterror[ix]:
                    consterror[ix] = tpe
        
        consterror = [r * 0.05 for r in consterror] #assuming 5% error
        
        #print varnames
        #print consterror
        FIM1, invFIM1 = fim.computeFIM(self.model, parsdict, varnames, sols, consterror)
        #print "compute FIM concluded"
        STDerrors = []
        for i,p in enumerate(parsdict.keys()):
            STDerrors.append((p,invFIM1[i,i]**0.5))    
        perrors =[]
        for p in pars:
            for p2,v in STDerrors:
                if p2 == p:
                    perrors.append(v)
                    break
        best['parameters']['data'] = [(self.model.uncertain[i].name, "%g"%value, "%g"%perrors[i]) for (i,value) in enumerate(self.bestSolution)]
        
        self.optimum = best

def reportResults(solver):
    reportText = ""
    sections = [solver.optimum[s] for s in ['parameters', 'optimization', 'timecourses']]
    for section in sections:
        reportText += "--- %-20s -----------------------------\n" % section['name'].upper()
        if section['header']:
            reportText += '\t\t'.join(section['header'])+'\n'
        reportText += "\n".join([section['format'] % i for i in section['data']])
        reportText += '\n\n'
    return reportText


def test():
    m1 = read_model("""
title Glyoxalase system in L. Infantum

glx1 : HTA -> SDLTSH, V1*HTA/(Km1 + HTA)
glx2 : SDLTSH ->,     V2*SDLTSH/(Km2 + SDLTSH)

V1  = 2.57594e-05
find V1  in [0.00001, 0.0001]

Km1 = 0.252531
find Km1 in [0.01, 1]

V2  = 2.23416e-05
find V2 in [0.00001, 0.0001]

Km2 = 0.0980973
find Km2 in (0.01, 1)

init = state(SDLTSH = 7.69231E-05, HTA = 0.1357)
""")
    
    #print m1
    optSettings={'genomesize':80, 'generations':200}
    timecourses = timecourse.readTCs(['TSH2a.txt', 'TSH2b.txt'], '../models', intvarsorder=(0,2,1))
    
    solver = DeODESolver(m1,optSettings, timecourses)
    
    time0 = time()
    
    solver.Solve()
    
    print "Optimization took %f s"% (time()-time0)

    print
    print '---------------------------------------------------------'
    print "Results for %s\n" % m1.getData('title')
    print reportResults(solver)

    #--- an example with unknown initial values --------------------
    
    m2 = m1.clone()
    
    # Now, assume init.HTA is uncertain
    m2.init.HTA.uncertainty(0.05,0.25)
    # do not estimate Km1 and Km2 to help the analysis
    m2.Km1.uncertainty(None)
    m2.Km2.uncertainty(None)
    
    optSettings={'genomesize':60, 'generations':200}
    
    ## VERY IMPORTANT:
    ## only one time course can be used: 
    ## cannot fit one uncertain initial value to several timecourses!!!
    timecourses = timecourse.readTCs(['TSH2a.txt'], '../models', intvarsorder=(0,2,1))
    
    solver = DeODESolver(m2,optSettings, timecourses)
    
    time0 = time()
    
    solver.Solve()
    
    print "Optimization took %f s"% (time()-time0)

    print
    print '---------------------------------------------------------'
    print "Results for %s\n" % m2.getData('title')
    print reportResults(solver)

if __name__ == "__main__":
    test()
