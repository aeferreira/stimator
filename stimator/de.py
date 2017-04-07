# This module is based on
#    PythonEquations is a collection of equations expressed as Python classes
#    Copyright (C) 2007 James R. Phillips
#    2548 Vera Cruz Drive
#    Birmingham, AL 35235 USA
#    email: zunzun@zunzun.com
#

"""DE : Real-value optimization by Differential evolution."""

from __future__ import print_function, absolute_import, division
import stimator.utils as utils
import random
import time
import numpy as np
import scipy.optimize

class DESolver(object):

    def __init__(self, pars_count, 
                       pop_size,
                       min_values, 
                       max_values, 
                       deStrategy, 
                       diffScale, 
                       crossoverProb, 
                       cutoffEnergy, 
                       useClassRandomNumberMethods,
                       max_generations=200,
                       convergence_noimprovement=20):

        random.seed(3)
        np.random.seed(3)

        self.max_generations = max_generations
        self.convergence_noimprovement = convergence_noimprovement
        self.pars_count = pars_count
        self.pop_size = pop_size
        self.cutoffEnergy = cutoffEnergy
        
        self.useClassRandomNumberMethods = bool(useClassRandomNumberMethods)
        if self.useClassRandomNumberMethods:
            self.SetupClassRandomNumberMethods()

        # deStrategy is the name of the DE function to use
        #self.calcTrialSolution = eval('self.' + deStrategy)
        self.calcTrialSolution = getattr(self, deStrategy)

        self.scale = diffScale
        self.crossOverProbability = crossoverProb
        self.min_values = min_values
        self.max_values = max_values

        # a random initial population, returns numpy arrays directly
        # min_values and max_values must be scalars
        # or vectors of size pars_count
        self.pop = np.random.uniform(0.0, 1.0, size=(pop_size, pars_count))
        self.pop = self.pop*(max_values-min_values)+min_values

        # initial energies for comparison
        self.popEnergy = np.ones(self.pop_size) * 1.0E300

        self.bestSolution = np.zeros(self.pars_count)
        self.bestEnergy = 1.0E300
        self.generation = 0
        self.generationsWithNoImprovement = 0
        self.atSolution = False
        self.exitCode = 0
        
        self.elapsed = 0.0
    
    
    exitCodeStrings = (
    "not done",
    "Solution found by energy criterium",
    "Solution found by diversity criterium",
    "Hit max generations ",
    "Too many generations with no improvement",
    "Solution found by convergence criterium")

    def reportInitialString (self):
        return "Solving..."

    def reportGenerationString (self):
        return "%-4d: %f" % (self.generation, self.bestEnergy)

    def reportFinalString(self):
        code = DESolver.exitCodeStrings[self.exitCode]
        res = ['Done!',
               '%s in %d generations.' % (code, self.generation),
               'best score = %f' % self.bestEnergy,
               'best solution: %s' % self.bestSolution]
        res = '\n' + '\n'.join(res)
        took_msg = '\nOptimization took {:.3f} s ({})'.format
        res += took_msg(self.elapsed, utils.s2HMS(self.elapsed))
        return res

    def reportInitial(self):
        print(self.reportInitialString())

    def reportGeneration(self):
        print(self.reportGenerationString())

    def reportFinal(self):
        print(self.reportFinalString())

    def GetRandIntInPars(self):
        return random.randint(0, self.pars_count-1)

    def GetRandFloatIn01(self):
        return random.uniform(0.0, 1.0)
        
    def GetRandIntInPop(self):
        return random.randint(0, self.pop_size-1)


    # this class might normally be subclassed and this method overridden,
    # or the self.externalEnergyFunction(trial) set
    # and this method used as is
    def EnergyFunction(self, trial):
        try:
            energy = self.externalEnergyFunction(trial)
        except (ArithmeticError, FloatingPointError):
            energy = 1.0E300 # high energies for arithmetic exceptions

        # we will be "done" if the energy is less than or equal to the cutoff energy
        if energy <= self.cutoffEnergy:
            return energy, True
        else:
            return energy, False

    def computeGeneration(self):
        # TODO: parallelization here
        

        # TODO: this is for performance on non-parallelized hardware
        if self.generationsWithNoImprovement > self.convergence_noimprovement:
            self.exitCode = 4
            return
                
        # Hit max generations
        if self.generation >= self.max_generations:
            self.exitCode = 3
            return
                
        # no need to try another generation if we are done (energy criterium)
        if self.atSolution:
            self.exitCode = 1
            return
        
        if self.generation == 0: #compute energies for generation 0
            self.reportInitial()
            for candidate in range(self.pop_size):
                trialEnergy, self.atSolution = self.EnergyFunction(np.copy(self.pop[candidate]))
                self.popEnergy[candidate] = trialEnergy
                if trialEnergy < self.bestEnergy:
                    self.bestEnergy = trialEnergy
                    self.bestSolution = np.copy(self.pop[candidate])
            #initialize stopwatch
            self.elapsed = time.clock()
        
        if (self.popEnergy.ptp()/self.popEnergy.mean()) < 1.0E-2:
            self.exitCode = 5
            return
           
            
        #print '==============================generation', self.generation
                
        for candidate in range(self.pop_size):
            if self.atSolution: break
                
            self.calcTrialSolution(candidate)
            trialEnergy, self.atSolution = self.EnergyFunction(self.trialSolution)
            
            if trialEnergy < self.popEnergy[candidate]:
                # New low for this candidate
                self.popEnergy[candidate] = trialEnergy
                self.pop[candidate] = np.copy(self.trialSolution)

                # Check if all-time low
                if trialEnergy < self.bestEnergy:
                    self.bestEnergy = trialEnergy
                    self.bestSolution = np.copy(self.trialSolution)
                    #self.bestSolution = self.trialSolution
                    self.generationsWithNoImprovement = 0
            
            #print self.pop[candidate],'=', self.popEnergy[candidate]
            
            # no need to try another candidate if we are done
            if self.atSolution:
                # it is possible for self.EnergyFunction()
                # to return self.atSolution == True even if we are not at
                # the best energy.  Just in case, copy the current values
                self.bestEnergy = trialEnergy
                self.bestSolution = np.copy(self.trialSolution)
                break # from candidate loop
            
        self.reportGeneration()
        if not self.atSolution:
            self.generation +=1
            self.generationsWithNoImprovement += 1
        return
        
    def finalize(self):
        # try to polish the best solution using scipy.optimize.fmin.
        if self.exitCode == 0:
            self.exitCode = -1
        if self.exitCode > 0:
            print ('refining last solution ...')
            self.bestSolution = scipy.optimize.fmin(self.externalEnergyFunction,
                                                    self.bestSolution,
                                                    disp=0) 
            self.bestEnergy, self.atSolution = self.EnergyFunction(self.bestSolution)
        self.elapsed = time.clock() - self.elapsed
        self.reportFinal()

    def run(self):
        while self.exitCode == 0:
            self.computeGeneration()
        self.finalize()

    
    # DE models
    def Best1Exp(self, candidate):
        r1,r2 = self.SelectSamples(candidate, 2)
        n = self.GetRandIntInPars()

        self.trialSolution = np.copy(self.pop[candidate])
        i = 0
        while(1):
            k = self.GetRandFloatIn01()
            if k >= self.crossOverProbability or i == self.pars_count:
                break
            self.trialSolution[n] = self.bestSolution[n] + self.scale * (self.pop[r1][n] - self.pop[r2][n])
            n = (n + 1) % self.pars_count
            i += 1


    def Rand1Exp(self, candidate):
        r1,r2,r3 = self.SelectSamples(candidate, 3)
        n = self.GetRandIntInPars()
        
        self.trialSolution = np.copy(self.pop[candidate])
        i = 0
        while(1):
            k = self.GetRandFloatIn01()
            if k >= self.crossOverProbability or i == self.pars_count:
                break
            self.trialSolution[n] = self.pop[r1][n] + self.scale * (self.pop[r2][n] - self.pop[r3][n])
            n = (n + 1) % self.pars_count
            i += 1

    def genIndxOfGenesToXover(self):
        #TODO this must be some discrete classic distribution random sample
        n = self.GetRandIntInPars()
        indx = np.zeros(self.pars_count,dtype=int)
        for i in range(self.pars_count):
            k = self.GetRandFloatIn01()
            if k >= self.crossOverProbability:
                break
            indx[n]=1
            n = (n + 1) % self.pars_count
        print (indx)
        return indx
        

    def Best2Exp(self, candidate):
        r1,r2,r3,r4 = self.SelectSamples(candidate, 4)
        self.trialSolution = np.copy(self.pop[candidate])
        n = self.GetRandIntInPars()
        #~ indx = self.genIndxOfGenesToXover()
        #~ print self.trialSolution
        #~ self.trialSolution[indx] = self.bestSolution[indx] + self.scale * (self.pop[r1][indx] + self.pop[r2][indx] - self.pop[r3][indx] - self.pop[r4][indx])
        for i in range(self.pars_count):
            #popn = self.pop[:,n]
            k = self.GetRandFloatIn01()
            if k >= self.crossOverProbability:
                break
            self.trialSolution[n] = self.bestSolution[n] + self.scale * (self.pop[r1][n] + self.pop[r2][n] - self.pop[r3][n] - self.pop[r4][n])
            #self.trialSolution[n] = self.bestSolution[n] + self.scale * (popn[r1] + popn[r2] - popn[r3] - popn[r4])
            n = (n + 1) % self.pars_count
            
    def RandToBest1Exp(self, candidate):
        r1,r2 = self.SelectSamples(candidate, 2)
        n = self.GetRandIntInPars()

        self.trialSolution = np.copy(self.pop[candidate])
        i = 0
        while(1):
            k = self.GetRandFloatIn01()
            if k >= self.crossOverProbability or i == self.pars_count:
                break
            self.trialSolution[n] += self.scale * (self.bestSolution[n] - self.trialSolution[n]) + self.scale * (self.pop[r1][n] - self.pop[r2][n])
            n = (n + 1) % self.pars_count
            i += 1


    def Rand2Exp(self, candidate):
        r1,r2,r3,r4,r5 = self.SelectSamples(candidate, 5)
        n = self.GetRandIntInPars()

        self.trialSolution = np.copy(self.pop[candidate])
        i = 0
        while(1):
            k = self.GetRandFloatIn01()
            if k >= self.crossOverProbability or i == self.pars_count:
                break
            self.trialSolution[n] = self.pop[r1][n] + self.scale * (self.pop[r2][n] + self.pop[r3][n] - self.pop[r4][n] - self.pop[r5][n])
            n = (n + 1) % self.pars_count
            i += 1

    def Best1Bin(self, candidate):
        r1,r2 = self.SelectSamples(candidate, 2)
        n = self.GetRandIntInPars()

        self.trialSolution = np.copy(self.pop[candidate])
        i = 0
        while(1):
            k = self.GetRandFloatIn01()
            if k >= self.crossOverProbability or i == self.pars_count:
                break
            self.trialSolution[n] = self.bestSolution[n] + self.scale * (self.pop[r1][n] - self.pop[r2][n])
            n = (n + 1) % self.pars_count
            i += 1


    def Rand1Bin(self, candidate):
        r1,r2,r3 = self.SelectSamples(candidate, 3)
        n = self.GetRandIntInPars()

        self.trialSolution = np.copy(self.pop[candidate])
        i = 0
        while(1):
            k = self.GetRandFloatIn01()
            if k >= self.crossOverProbability or i == self.pars_count:
                break
            self.trialSolution[n] = self.pop[r1][n] + self.scale * (self.pop[r2][n] - self.pop[r3][n])
            n = (n + 1) % self.pars_count
            i += 1


    def RandToBest1Bin(self, candidate):
        r1,r2 = self.SelectSamples(candidate, 2)
        n = self.GetRandIntInPars()

        self.trialSolution = np.copy(self.pop[candidate])
        i = 0
        while(1):
            k = self.GetRandFloatIn01()
            if k >= self.crossOverProbability or i == self.pars_count:
                break
            self.trialSolution[n] += self.scale * (self.bestSolution[n] - self.trialSolution[n]) + self.scale * (self.pop[r1][n] - self.pop[r2][n])
            n = (n + 1) % self.pars_count
            i += 1


    def Best2Bin(self, candidate):
        r1,r2,r3,r4 = self.SelectSamples(candidate, 4)
        n = self.GetRandIntInPars()

        self.trialSolution = np.copy(self.pop[candidate])
        i = 0
        while(1):
            k = self.GetRandFloatIn01()
            if k >= self.crossOverProbability or i == self.pars_count:
                break
            self.trialSolution[n] = self.bestSolution[n] + self.scale * (self.pop[r1][n] + self.pop[r2][n] - self.pop[r3][n] - self.pop[r4][n])
            n = (n + 1) % self.pars_count
            i += 1


    def Rand2Bin(self, candidate):
        r1,r2,r3,r4,r5 = self.SelectSamples(candidate, 5)
        n = self.GetRandIntInPars()

        self.trialSolution = np.copy(self.pop[candidate])
        i = 0
        while(1):
            k = self.GetRandFloatIn01()
            if k >= self.crossOverProbability or i == self.pars_count:
                break
            self.trialSolution[n] = self.pop[r1][n] + self.scale * (self.pop[r2][n] + self.pop[r3][n] - self.pop[r4][n] - self.pop[r5][n])
            n = (n + 1) % self.pars_count
            i += 1

    def SelectSamples(self, candidate, n):
        """Select n different members of population which are different from candidate."""
        
        s = random.sample(list(range(self.pop_size)),n)
        while candidate in s:
            s = random.sample(list(range(self.pop_size)),n)
        
#        universe = range(0,candidate)+range(candidate+1, self.pop_size-1)
#        s = random.sample(universe, n)
        return s

    def SetupClassRandomNumberMethods(self):
        np.random.seed(3) # this yields same results each time run() is run
        self.nonStandardRandomCount = self.pop_size * self.pars_count * 3
        if self.nonStandardRandomCount < 523: # set a minimum number of random numbers
            self.nonStandardRandomCount = 523
            
        self.ArrayOfRandomIntegersBetweenZeroAndParameterCount = \
        np.random.random_integers(0, self.pars_count-1, size=(self.nonStandardRandomCount))
        self.ArrayOfRandomRandomFloatBetweenZeroAndOne = \
        np.random.uniform(size=(self.nonStandardRandomCount))
        self.ArrayOfRandomIntegersBetweenZeroAndPopulationSize \
        = np.random.random_integers(0, self.pop_size-1, size=(self.nonStandardRandomCount))
        self.randCounter1 = 0
        self.randCounter2 = 0
        self.randCounter3 = 0


    def GetClassRandomIntegerBetweenZeroAndParameterCount(self):
        self.randCounter1 += 1
        if self.randCounter1 >= self.nonStandardRandomCount:
            self.randCounter1 = 0
        return self.ArrayOfRandomIntegersBetweenZeroAndParameterCount[self.randCounter1]

    def GetClassRandomFloatBetweenZeroAndOne(self):
        self.randCounter2 += 1
        if self.randCounter2 >= self.nonStandardRandomCount:
            self.randCounter2 = 0
        return self.ArrayOfRandomRandomFloatBetweenZeroAndOne[self.randCounter2]
        
    def GetClassRandomIntegerBetweenZeroAndPopulationSize(self):
        self.randCounter3 += 1
        if self.randCounter3 >= self.nonStandardRandomCount:
            self.randCounter3 = 0
        return self.ArrayOfRandomIntegersBetweenZeroAndPopulationSize[self.randCounter3]

