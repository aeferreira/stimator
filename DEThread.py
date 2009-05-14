#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2009 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""

import thread
import time
import DESolver
import wx
from numpy import *
from scipy import integrate

DUMP_PARS_2FILES = False

##------------- Computing thread class

class DeOdeSolver(DESolver.DESolver):
    """Overides energy function and report functions.
    
    The energy function solves ODEs and computes a least-squares score.
    Report functions post events to the main window (nothing else)."""

    def setup(self, model, win, atimecoursecollection, anUpdateGenerationEvent, anEndComputationEvent):
        self.model = model
        self.win = win # the widget that receives the notifications
        self.tc = atimecoursecollection
        self.timecoursedata = self.tc.data
        self.UpdateGenerationEvent = anUpdateGenerationEvent
        self.EndComputationEvent = anEndComputationEvent
        
        # cutoffEnergy is 0.1% of deviation from data
        self.cutoffEnergy =  1.0e-6*sum([nansum(abs(tc[:,1:])) for tc in self.timecoursedata])
        
        # scale times to maximum time in data
        scale = float(max([ (tc[-1,0]-tc[0,0]) for tc in self.timecoursedata]))
        
        self.calcDerivs = model.scaled_dXdt(scale)
        
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
        
        self.model.set_unknown(trial)
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
        evt = self.UpdateGenerationEvent(generation = self.generation, energy = float(self.bestEnergy))
        wx.PostEvent(self.win, evt)
        # write parameters to files
        if DUMP_PARS_2FILES:
            for par in range(self.parameterCount):
                parvector = [str(self.population[k][par]) for k in range(self.populationSize)]
                print >>self.parfilehandes[par], " ".join(parvector)
            
    
    def reportFinal (self):
        if self.exitCode==0: outCode = -1 
        else: 
            outCode = self.exitCode
            self.win.bestData = self.generateBestData() # bestData is own by main window
        evt = self.EndComputationEvent(exitCode = outCode)
        wx.PostEvent(self.win, evt)
        # close parameter files
        if DUMP_PARS_2FILES:
            for par in self.parfilehandes:
                par.close()

    def generateBestData (self):
        timecoursedata = self.timecoursedata
        bestData = [{'section':"PARAMETERS"}, 
                    {'section':"D.E. OPTIMIZATION"}, 
                    {'section':"TIMECOURSES"},
                    {'section':"best timecourses"}]
        bestData[0]['data'] = [(self.model.unknown[i].name, "%g"%value) for (i,value) in enumerate(self.bestSolution)]
        bestData[0]['format'] = "%s\t%s"
        bestData[0]['header'] = None
        
        bestData[1]['data'] = [('Final Score', "%g"% self.bestEnergy),
                                ('Generations', "%d"% self.generation),
                                ('Exit by    ', "%s"% self.exitCodeStrings[self.exitCode])]
        bestData[1]['format'] = "%s\t%s"
        bestData[1]['header'] = None

        #TODO: Store initial solver parameters?

        #generate best time-courses
        bestData[2]['data'] = []
        bestData[3]['data'] = []
        besttimecoursedata = bestData[3]['data']
        for (i,data) in enumerate(timecoursedata):
            y0 = copy(self.X0[i])
            t  = self.times[i]
            Y, infodict = integrate.odeint(self.calcDerivs, y0, t, full_output=True, printmessg=False)
            besttimecoursedata.append(Y)
            if infodict['message'] != 'Integration successful.':
                bestData[2]['data'].append((self.tc.shortnames[i], self.tc.shapes[i][0], 1.0E300))
            else:
                S = (Y- data[:, 1:])**2
                score = nansum(S)
                bestData[2]['data'].append((self.tc.shortnames[i], self.tc.shapes[i][0], score))
        bestData[2]['format'] = "%s\t%d\t%g"
        bestData[2]['header'] = ['Name', 'Points', 'Score']
        
        return bestData
        

class CalcOptmThread:
    def __init__(self, win):
        self.win = win 

    def Start(self, model, optSettings, timecoursecollection, anUpdateGenerationEvent, anEndComputationEvent):
        mins = array([u.min for u in model.unknown])
        maxs = array([u.max for u in model.unknown])
        
        self.solver = DeOdeSolver(len(model.unknown), # number of parameters
                                 int(optSettings['genomesize']),  # genome size
                                 int(optSettings['generations']), # max number of generations
                                 mins, maxs,              # min and max parameter values
                                 "Best2Exp",              # DE strategy
                                 0.7, 0.6, 0.0,           # DiffScale, Crossover Prob, Cut off Energy
                                 True)                    # use class random number methods

        self.solver.setup(model,self.win,timecoursecollection, anUpdateGenerationEvent, anEndComputationEvent)
        
        self.keepGoing = self.running = True
        thread.start_new_thread(self.Run, ())

    def Stop(self):
        self.keepGoing = False

    def IsRunning(self):
        return self.running

    def Run(self):
        while self.keepGoing:
            self.solver.computeGeneration()
            if self.solver.exitCode !=0: self.Stop()

        self.solver.finalize()
        self.running = False

