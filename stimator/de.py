# This module is based on
#    PythonEquations is a collection of equations expressed as Python classes
#    Copyright (C) 2007 James R. Phillips
#    2548 Vera Cruz Drive
#    Birmingham, AL 35235 USA
#    email: zunzun@zunzun.com
#

"""DE : Real-value optimization by Differential evolution."""

from __future__ import print_function, absolute_import, division
from . import utils
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
                       cutoff_score,
                       max_generations=200,
                       convergence_noimprovement=20,
                       useClassRandomNumberMethods=False):

        np.random.seed(3)

        self.max_generations = max_generations
        self.convergence_noimprovement = convergence_noimprovement
        self.pars_count = pars_count
        self.pop_size = pop_size
        self.cutoff_score = cutoff_score
        
        self.useClassRandomNumberMethods = bool(useClassRandomNumberMethods)
        if self.useClassRandomNumberMethods:
            self.SetupClassRandomNumberMethods()

        # deStrategy is the name of the DE function to use
        self.calcTrialSolution = getattr(self, deStrategy)

        # DE hyperparameters
        self.scale = diffScale
        self.crossOverProbability = crossoverProb
        
        # bounds
        self.min_values = min_values
        self.max_values = max_values

        # a random initial population, returns numpy arrays directly
        # min_values and max_values must be scalars
        # or vectors of size pars_count
        self.pop = np.random.uniform(0.0, 1.0, size=(pop_size, pars_count))
        self.pop = self.pop * (max_values - min_values) + min_values

        # initial scores for comparison
        self.scores = np.ones(self.pop_size) * float('inf')

        self.best = np.empty(self.pars_count)
        self.best_score = float('inf')
        self.generation = 0
        self.generationsWithNoImprovement = 0
        self.atSolution = False
        self.exitCode = 0

    exitCodeStrings = (
    "not done",
    "Solution found by energy criterium",
    "Solution found by diversity criterium",
    "Hit max generations ",
    "Too many generations with no improvement",
    "Solution found by convergence criterium")

    def reportInitial(self):
        print("Solving...")
        #initialize stopwatch
        self.start_time = time.clock()

    def reportGeneration(self):
        print("{:-4d}: {:f}".format(self.generation, self.best_score))

    def reportFinal(self):
        ttime = time.clock() - self.start_time
        code = DESolver.exitCodeStrings[self.exitCode]
        res = ['Done!',
               '%s in %d generations.' % (code, self.generation),
               'best score = %f' % self.best_score,
               'best solution: %s' % self.best]
        res = '\n' + '\n'.join(res)
        took_msg = '\nOptimization took {:.3f} s ({})'.format
        res += took_msg(ttime, utils.s2HMS(ttime))
        print(res)

    # this class might normally be subclassed and this method overridden,
    # or the self.external_score_function(individual) set
    # and this method used as is
    def score_function(self, individual):
        try:
            score = self.external_score_function(individual)
        except (ArithmeticError, FloatingPointError):
            score = float('inf') # high energies for computational exceptions

        # we will be "done" if the score is
        # less than or equal to the cut-off score
        if score <= self.cutoff_score:
            return score, True
        else:
            return score, False

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
                
        
        if self.generation == 0: #compute energies for generation 0
            self.reportInitial()
            for i in range(self.pop_size):
                score, self.atSolution = self.score_function(np.copy(self.pop[i]))
                self.scores[i] = score
                if score < self.best_score:
                    self.best_score = score
                    self.best = np.copy(self.pop[i])
            self.reportGeneration()
            if not self.atSolution:
                self.generationsWithNoImprovement += 1
        
        # no need to try another generation if we are done (energy criterium)
        if self.atSolution:
            self.exitCode = 1
            return

        # return if whole population "converged"
        if (self.scores.ptp()/self.scores.mean()) < 1.0E-2:
            self.exitCode = 5
            return

        self.generation +=1

        for i in range(self.pop_size):
            if self.atSolution: break
                
            new_i = self.calcTrialSolution(i)
            if new_i is not None:
                score, self.atSolution = self.score_function(new_i)
                
                if score < self.scores[i]:
                    # New low for this individual
                    self.pop[i], self.scores[i] = new_i, score

                    # Check if all-time low
                    if score < self.best_score:
                        self.best, self.best_score = new_i, score
                        self.generationsWithNoImprovement = 0

            # no need to try another i if we are done
            if self.atSolution:
                # it is possible for self.score_function()
                # to return self.atSolution == True even if we are not at
                # the best score. Copy the current values
                self.best, self.best_score = new_i, score
            
        self.reportGeneration()
        if not self.atSolution:
            self.generationsWithNoImprovement += 1
        return
        
    def finalize(self):
        # try to polish the best solution using scipy.optimize.fmin.
        if self.exitCode == 0:
            self.exitCode = -1
        if self.exitCode > 0:
            print ('refining last solution ...')
            self.best = scipy.optimize.fmin(self.external_score_function,
                                            self.best, disp=0) 
            self.best_score, self.atSolution = self.score_function(self.best)
        self.reportFinal()

    def run(self):
        while self.exitCode == 0:
            self.computeGeneration()
        self.finalize()

    
    # DE models

    def get_pars2change(self):
        pars = []
        n = np.random.randint(self.pars_count)
        for i in range(self.pars_count):
            d = np.random.uniform()
            if d >= self.crossOverProbability:
                return pars
            pars.append(n)
            n = (n + 1) % self.pars_count
        return pars


    def Best1Exp(self, i):
        r1,r2 = self.SelectSamples(i, 2)
        pop = self.pop
        new_i = np.copy(pop[i])
        change = self.get_pars2change()
        if len(change) == 0:
            return None
        for n in change:
            new_i[n] = self.best[n] + self.scale * (pop[r1][n] - pop[r2][n])
        return new_i

    def Rand1Exp(self, i):
        r1,r2,r3 = self.SelectSamples(i, 3)
        pop = self.pop
        new_i = np.copy(pop[i])
        change = self.get_pars2change()
        if len(change) == 0:
            return None
        for n in change:
            new_i[n] = pop[r1][n] + self.scale * (pop[r2][n] - pop[r3][n])
        return new_i
    
    def Best2Exp(self, i):
        r1,r2,r3,r4 = self.SelectSamples(i, 4)
        pop = self.pop
        new_i = np.copy(pop[i])
        change = self.get_pars2change()
        if len(change) == 0:
            return None
        for n in change:
            new_i[n] = self.best[n] + self.scale * (pop[r1][n] + pop[r2][n] - pop[r3][n] - pop[r4][n])
        return new_i

    def RandToBest1Exp(self, i):
        r1,r2 = self.SelectSamples(i, 2)
        pop = self.pop
        new_i = np.copy(pop[i])
        change = self.get_pars2change()
        if len(change) == 0:
            return None
        for n in change:
            new_i[n] = self.scale * (self.best[n] - new_i[n]) + self.scale * (pop[r1][n] - pop[r2][n])
        return new_i

    def Rand2Exp(self, i):
        r1,r2,r3,r4,r5 = self.SelectSamples(i, 5)
        pop = self.pop
        new_i = np.copy(pop[i])
        change = self.get_pars2change()
        if len(change) == 0:
            return None
        for n in change:
            new_i[n] = pop[r1][n] + self.scale * (pop[r2][n] + pop[r3][n] - pop[r4][n] - pop[r5][n])
        return new_i

    def Best1Bin(self, i):
        r1,r2 = self.SelectSamples(i, 2)
        new_i = np.copy(self.pop[i])
        change = self.get_pars2change()
        if len(change) == 0:
            return None
        for n in change:
            new_i[n] = self.best[n] + self.scale * (self.pop[r1][n] - self.pop[r2][n])
        return new_i

    def Rand1Bin(self, i):
        r1,r2,r3 = self.SelectSamples(i, 3)
        pop = self.pop
        new_i = np.copy(pop[i])
        change = self.get_pars2change()
        if len(change) == 0:
            return None
        for n in change:
            new_i[n] = pop[r1][n] + self.scale * (pop[r2][n] - pop[r3][n])
        return new_i

    def RandToBest1Bin(self, i):
        r1,r2 = self.SelectSamples(i, 2)
        new_i = np.copy(self.pop[i])
        change = self.get_pars2change()
        if len(change) == 0:
            return None
        for n in change:
            new_i[n] += self.scale * (self.best[n] - new_i[n]) + self.scale * (self.pop[r1][n] - self.pop[r2][n])
        return new_i

    def Best2Bin(self, i):
        r1,r2,r3,r4 = self.SelectSamples(i, 4)
        new_i = np.copy(self.pop[i])
        change = self.get_pars2change()
        if len(change) == 0:
            return None
        for n in change:
            new_i[n] = self.best[n] + self.scale * (self.pop[r1][n] + self.pop[r2][n] - self.pop[r3][n] - self.pop[r4][n])
        return new_i

    def Rand2Bin(self, i):
        r1,r2,r3,r4,r5 = self.SelectSamples(i, 5)
        new_i = np.copy(self.pop[i])
        change = self.get_pars2change()
        if len(change) == 0:
            return None
        for n in change:
            new_i[n] = self.pop[r1][n] + self.scale * (self.pop[r2][n] + self.pop[r3][n] - self.pop[r4][n] - self.pop[r5][n])
        return new_i

    def SelectSamples(self, i, n):
        """Select n different indexes of pop which are different from i."""
        #the_set = [k for k in range(self.pop_size) if k != i]
        #the_set = list(range(i)) + list(range(i + 1, self.pop_size))
        the_set = range(self.pop_size)
        s = np.random.choice(the_set, n, replace=False)
        if i in s:
            s = np.random.choice(the_set, n, replace=False)
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

