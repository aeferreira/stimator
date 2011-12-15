#!/usr/bin/env python
# -*- coding:latin1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2010 Ant�nio Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

import de
from numpy import *
from scipy import integrate
from model import *
from dynamics import *
from modelparser import read_model
from analysis import *
import fim
import timecourse
import dyncriteria
import expcov

#----------------------------------------------------------------------------
#         Class to perform DE optimization for ODE systems
#----------------------------------------------------------------------------

class OptimumData(object):
    pass

class DeODESolver(de.DESolver):
    """Overides energy function and report functions.
    
    The energy function solves ODEs and computes a least-squares score.
    Ticker functions are called on generation completion and when optimization finishes.
    """
    
    def __init__(self, model, optSettings, tcs, weights = None,
                    aMsgTicker=None, anEndComputationTicker=None, 
                    dump_pars=False):
        self.model = model
        self.tc    = tcs
        self.endTicker        = anEndComputationTicker
        self.msgTicker        = aMsgTicker
        self.dump_pars        = dump_pars
        
        #reorder variables according to model
        self.tc.orderByModelVars(self.model)

        pars = uncertain(model)
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
        t0 = self.tc[0].t[0]
        
        self.calcDerivs = getdXdt(model, scale=scale, with_uncertain=True, t0=t0)
        self.salg=integrate._odepack.odeint
        
        # store initial values and (scaled) time points
        self.X0 = []
        self.times = []
        for data in self.tc:
            self.X0.append(copy(data[:, 0].T))
            t  = data.t
            times = (t-t0)/scale #+t0  # this scales time points
            self.times.append(times)
        self.timecourse_scores = empty(len(self.tc))
        
        # find uncertain initial values
        mapinit2trial = []
        for iu, u in enumerate(uncertain(self.model)):
            if get_name(u).startswith('init'):
                varname = get_name(u).split('.')[-1]
                ix = varnames(self.model).index(varname)
                mapinit2trial.append((ix,iu))
        self.trial_initindexes = array([j for (i,j) in mapinit2trial], dtype=int)
        self.vars_initindexes = array([i for (i,j) in mapinit2trial], dtype=int)
        
        self.criterium = dyncriteria.getCriteriumFunction(weights, self.tc)

        # open files to write parameter progression
        if self.dump_pars:
            self.parfilehandes = [open(par[0]+".par", 'w') for par in model.parameters]

    def computeSolution(self,i,trial):
        """Computes solution for timecourse i, given parameters trial."""
        
        y0 = copy(self.X0[i])
        # fill uncertain initial values
        y0[self.vars_initindexes] = trial[self.trial_initindexes]
            
##         t  = copy(self.times[i])
##         Y, infodict = integrate.odeint(self.calcDerivs, y0, t, full_output=True, printmessg=False)
        output = self.salg(self.calcDerivs, y0, self.times[i], (), None, 0, -1, -1, 0, None, 
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

    def reportInitial (self):
        msg = "\nSolving %s..."%self.model['title']
        if not self.msgTicker:
            print msg
        else:
            self.msgTicker(msg)

    def reportGeneration (self):
        msg = "%-4d: %f" % (self.generation, float(self.bestEnergy))
        if not self.msgTicker:
            print msg
        else:
            self.msgTicker(msg)
        if self.dump_pars:
            for par in range(self.parameterCount):
                parvector = [str(self.population[k][par]) for k in range(self.populationSize)]
                print >>self.parfilehandes[par], " ".join(parvector)


##     def reportGeneration (self):
##         if not self.generationTicker:
##             print self.reportGenerationString()
##         else:
##             self.generationTicker(self.generation, float(self.bestEnergy))
##         if self.dump_pars:
##             for par in range(self.parameterCount):
##                 parvector = [str(self.population[k][par]) for k in range(self.populationSize)]
##                 print >>self.parfilehandes[par], " ".join(parvector)
            
    
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
        best = OptimumData()
        best.optimization_score = self.bestEnergy
        best.optimization_generations = self.generation
        best.optimization_exit_by = self.exitCodeStrings[self.exitCode]

        #TODO: Store initial solver parameters?

        #generate best time-courses
        best.tcdata = []

        modelvnames = varnames(self.model)
        pars = [get_name(uncertain(self.model)[i]) for i in range(len(self.bestSolution))]
        parvalues = [value for value in self.bestSolution]
        parszip = zip(pars, parvalues)
        self.model.set_uncertain(self.bestSolution)
        
        sols = timecourse.Solutions()
        
        for (i,tc) in enumerate(self.tc):
            Y = self.computeSolution(i, self.bestSolution)
            if Y is not None:
                score =self.criterium(Y, i)
            else:
                score = 1.0E300
            nt = tc.ntimes
            sol = timecourse.SolutionTimeCourse (tc.t, Y.T, modelvnames)
            sols += sol
            best.tcdata.append((self.tc.shortnames[i], self.tc[i].ntimes, score))
            
            vnames = []
            varindexes=[]
            for iline,line in enumerate(tc.data):
                #count NaN
                yexp = line
                nnan = len(yexp[isnan(yexp)])
                if nnan >= nt-1: continue
                vnames.append(modelvnames[iline])
                varindexes.append(iline)
##             print vnames
##             print varindexes
            #print len(varindexes), varnames
        best.optimum_tcs=sols
        
        
        if not (fim.sympy_installed):
            best.parameters = [(get_name(self.model.uncertain[i]), "%g"%value, "0.0") for (i,value) in enumerate(self.bestSolution)]
        else:
            consterror = [0.0 for i in range(len(vnames))]
            for ix, x in enumerate(vnames):
                for tc in self.tc:
                    yexp = tc.data[varindexes[ix]]
                    tpe = (max(yexp) - min(yexp))
                    if tpe > consterror[ix]:
                        consterror[ix] = tpe
            consterror = expcov.constError_func([r * 0.05 for r in consterror]) #assuming 5% error
            
            #print consterror
            FIM1, invFIM1 = fim.computeFIM(self.model, parszip, vnames, sols, consterror)
            
            STDerrors = {}
            for i,p in enumerate(pars):
                STDerrors[p] =invFIM1[i,i]**0.5
            best.parameters = [(get_name(uncertain(self.model)[i]), "%g"%value, "%g"%STDerrors[get_name(uncertain(self.model)[i])]) for (i,value) in enumerate(self.bestSolution)]
        
        self.optimum = best

    def reportResults(self):
        headerformat = "--- %-20s -----------------------------\n"
        reportText = "\n"
        reportText += headerformat % 'PARAMETERS'
        reportText += "\n".join(["%s\t%12s +- %s" % i for i in self.optimum.parameters])
        reportText += '\n\n'
        reportText += headerformat % 'OPTIMIZATION'
        reportText += "%s\t%g\n" % ('Final Score', self.optimum.optimization_score)
        reportText += "%s\t%g\n" % ('generations', self.optimum.optimization_generations)
        reportText += "%s\t%s\n" % ('Exit by',     self.optimum.optimization_exit_by)
        reportText += '\n\n'
        reportText += headerformat % 'TIME COURSES'
        reportText += '\t\t'.join(['Name', 'Points', 'Score'])+'\n'
        reportText += "\n".join(["%s\t%d\t%g" % i for i in self.optimum.tcdata])
        reportText += '\n\n'
        
        return reportText

    def draw(self, figure):
        figure.clear()
        tcsubplots = []
        bestsols = self.optimum.optimum_tcs
        expsols = self.tc
        tcstats = self.optimum.tcdata
        ntc = len(bestsols)
        ncols = int(math.ceil(math.sqrt(ntc)))
        nrows = int(math.ceil(float(ntc)/ncols))
        for i in range(ntc):
            tcsubplots.append(figure.add_subplot(nrows,ncols,i+1))

        for i in range(ntc):
            subplot = tcsubplots[i]
            #subplot.set_xlabel("time")
            subplot.set_title("%s (%d pt) %g"% tcstats[i], fontsize = 12)
            expsol = expsols[i]
            symsol = bestsols[i]
            for line in range(len(expsol)):
                #count NaN and do not plot if they are most of the timecourse
                yexp = expsol[line]
                nnan = len(yexp[isnan(yexp)])
                if nnan >= expsol.ntimes-1: continue
                #otherwise plot lines
                ysim = symsol[line]
                subplot.plot(expsol.t,yexp, '-b')
                subplot.plot(expsol.t,ysim, '-r')

def test():
    m1 = read_model("""
title Glyoxalase system in L. Infantum

glx1 : HTA -> SDLTSH, V1*HTA/(Km1 + HTA)
#glx1 : HTA -> SDLTSH, V*HTA/(Km1 + HTA), V=2.57594e-05
glx2 : SDLTSH ->,     V2*SDLTSH/(Km2 + SDLTSH)

#find glx1.V  in [0.00001, 0.0001]
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
    timecourses = timecourse.readTCs(['TSH2a.txt', 'TSH2b.txt'], 'examples', names = ['SDLTSH', 'HTA'], verbose = True)
    #intvarsorder=(0,2,1), verbose=True)
    
    solver = DeODESolver(m1,optSettings, timecourses)
    solver.Solve()
    
    print solver.reportResults()

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
    timecourses = timecourse.readTCs(['TSH2a.txt'], 'examples', names = ['SDLTSH', 'HTA'], verbose = True)
    #, intvarsorder=(0,2,1), verbose=True)
    
    solver = DeODESolver(m2,optSettings, timecourses)
    solver.Solve()

    print solver.reportResults()

if __name__ == "__main__":
    test()
