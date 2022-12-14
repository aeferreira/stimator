from __future__ import print_function, absolute_import, division
from time import time
import numpy as np
from scipy import integrate

import stimator.timecourse as timecourse
import stimator.dynamics as dynamics
import stimator.utils as utils
from stimator.de import DESolver
from stimator.moo.rmc import remove_most_crowded
from stimator.moo.sorting import MOOSorter

def dominance(vec1, vec2):
    """Compute Pareto dominance relationship."""
    d_result = 0
    for vo,vn in zip(vec1, vec2):
        d = vn-vo
        if d <= 0 and d_result <=0:
            d_result -= 1
        elif d >= 0 and d_result >=0:
            d_result += 1
        else:
            return 0
    return d_result

## def dominance(vec1, vec2):
##     """Compute Pareto dominance relationship."""
##     size = len(vec1)
##     d = vec2 <= vec1
##     if np.all(d): return -size
##     d = vec2 >= vec1
##     if np.all(d): return size
##     return 0

def nondominated_solutions(energies):
    """Returns the indexes of non-dominated solutions in a population."""
    nondominated = []
    for i in range(len(energies)):
        is_dominated = False
        for k in range(len(energies)):
            d_result = 0
            for j in range(len(energies[i])):
                d = energies[i][j] - energies[k][j]
                if d <= 0 and d_result <=0:
                    d_result -= 1
                elif d >= 0 and d_result >=0:
                    d_result += 1
                else:
                    d_result = 0
                    break
            if d_result > 0:
                is_dominated = True
                break
        if not is_dominated:
            nondominated.append(i)
    return nondominated

def dominance_delta(old_energies, new_energies):
    """Compute dominance relationship between solutions of two generations."""
    return [dominance(o,n) for o,n in zip(old_energies, new_energies)]


class ModelSolver(object):
    """Driver to compute timecourses, given a model and a trial vector.

    'model' is a model object
    't0' and 'tf' are floats
    'npoints' is a positive integer.
    'observed' is a list of variable names to output.
    'optnames'are variable names. Their initial values will be optimized
    in the experimental design."""

    def __init__(self, model, t0, npoints, tf, optnames, observed):
        self.t0 = t0
        self.npoints = npoints
        if tf is None:
            tf = 1.0
        self.tf = tf

        self.vector = np.copy(dynamics.init2array(model))

        self.names = model.varnames

        self.model = model
        self.optvars = utils.listify(optnames)
        self.observed = utils.listify(observed)

        #check that the names exist as variables in the model
        for name in self.optvars+self.observed:
            if not (name in self.names):
                raise AttributeError('%s is not a variable in model'%name)

        self.optvars_indexes = np.array([self.names.index(name) for name in self.optvars])
        self.obsvars_indexes = np.array([self.names.index(name) for name in self.observed])

        self.times = np.linspace (self.t0, self.tf, self.npoints)                
        
        # scale times for output
        scale = float(self.times[-1] - self.t0)
        #scale = 1.0
        
        self.f = dynamics.getdXdt(model, scale=scale, t0=self.t0)
        self.t = np.copy((self.times-self.t0)/scale) # this scales time points

    def solve(self, trial):
        """Returns the solution for a particular trial of initial values."""
        #salg=integrate._odepack.odeint

        for value,indx in zip(trial, self.optvars_indexes):
            self.vector[indx] = value
        y0 = np.copy(self.vector)
        
        
        solver = integrate.odeint
        output = solver(self.f, y0, self.t,
                        args=(),
                        Dfun=None,
                        col_deriv=0,
                        full_output=True,
                        ml=None,
                        rtol=None,
                        mu=None,
                        atol=None,
                        tcrit=None, 
                        h0=0.0, 
                        hmax=0.0,
                        hmin=0.0,
                        ixpr=0,
                        mxstep=0,
                        mxhnil=0,
                        mxordn=12,
                        mxords=5)#, tfirst=False)
        out_message = output[1]['message'].strip()
        if out_message != 'Integration successful.':
            #print('Solution failed:', out_message)
            return None

        Y = output[0]
        title = self.model.metadata.get('title', '')
        Y = np.copy(Y.T)

        sol = timecourse.SolutionTimeCourse(self.times, Y, self.names, title)
        sol = sol.copy(names=self.observed)
        return sol


class GDE3Solver(DESolver):
    """ 
    Adaptation of DESolver for multiobjective optimization.
    Based on 
    Kukkonen, S. and Lampinen, J. (2005) GDE3: 
    The third step of generalized differential evolution, 
    2005 IEEE Congress on Evolutionary Computation. 
    """

    def __init__(self, models, 
                 toOpt, objectiveFunction, observed, 
                 npoints, t0, tf, 
                 pop_size, deStrategy, 
                 diffScale, crossoverProb, cutoffEnergy, 
                 useClassRandomNumberMethods=False, 
                 dif='0',
                 dump_generations=None,
                 max_generations=200):

        self.models = models
        self.npoints = npoints
        self.t0 = t0
        self.tf = tf
        self.observed = observed

        self.nmodels = len(models)

        self.toOpt = toOpt
        self.toOptKeys = []
        self.opt_maxs = []
        self.opt_mins = []
        for n, min_v, max_v in toOpt:
            self.toOptKeys.append(n)
            self.opt_maxs.append(max_v)
            self.opt_mins.append(min_v)
        self.opt_maxs = np.array(self.opt_maxs)
        self.opt_mins = np.array(self.opt_mins)

        #self.toOptKeys = self.toOpt.keys()
        self.objFunc = objectiveFunction
        self.dump_generations = dump_generations

        DESolver.__init__(self, 
                          len(toOpt), pop_size,
                          self.opt_mins, self.opt_maxs, 
                          deStrategy, diffScale, crossoverProb, 
                          cutoffEnergy, 
                          max_generations=max_generations,
                          useClassRandomNumberMethods=useClassRandomNumberMethods)

        self.deltaT = (tf - t0) / npoints

        self.objFuncList = []

        for m in self.models:
            self.objFuncList.append(ModelSolver(m, 
                                                self.t0, 
                                                self.npoints, 
                                                self.tf, 
                                                self.toOptKeys, 
                                                self.observed))

        self.dif = dif

        # threshold for improvement based on 5 % of new solution count
        self.roomForImprovement = int(round(0.05 * pop_size))

        #counter of number of generations with only non-dominated solutions
        self.fullnondominated = 0

        str2distance = {'extKL'   :timecourse.extendedKLdivergence,
                        'L2_midpoint_weights':timecourse.L2_midpoint_weights,
                        'KL'      :timecourse.KLdivergence,
                        'L2'      :timecourse.L2}
        if self.objFunc not in str2distance:
            raise ("%s is not an implemented divergence function"%self.objFunc)

        if self.objFunc in ('L2_midpoint_weights','L2'):  #symetric measures
            self.trueMetric = True
        else:
            self.trueMetric = False

        self.distance_func = str2distance[self.objFunc]

        # generate list of pairs of model indexes. Each pair corresponds to a different comparison
        self.model_indexes = []
        for i in range(self.nmodels-1):
            for j in range(i+1, self.nmodels):
                self.model_indexes.append((i,j))
        # if not a true metric, include also the symetric pairs
        if not self.trueMetric:
            self.model_indexes.extend([(j,i) for (i,j) in self.model_indexes])

        # working storage arrays
        self.new_generation_energies = [[] for i in range(self.pop_size)]
        self.new_population = np.empty((self.pop_size,len(self.toOpt)))
        self.gen_times = []
        self.calcTrialSolution = self.Rand1Bin

    def Rand1Bin(self, i):
        """This function overrides the base class (de.py).
           The intent is to aggressively generate a new solution from
           current population."""
        r1,r2,r3 = self.SelectSamples(i, 3)
        new_i = np.copy(self.pop[i])
        for n in range(self.pars_count):
            new_i[n] = self.pop[r1][n] + self.scale * (self.pop[r2][n] - self.pop[r3][n])
        return new_i

    def score_function(self, trial):
        #compute solution for each model, using trial vector
        sols = [s.solve(trial) for s in self.objFuncList]
        if self.distance_func is not None:
            return self.distance_func(sols, self.deltaT, self.model_indexes)
        return sols

    def computeGeneration(self):
        # Hit max generations with no improvement
        if self.generationsWithNoImprovement > 20:
            self.exitCode = 4
            return
        # Hit max generations
        if self.generation >= self.max_generations:
            self.exitCode = 3
            return
        # Hit several generations with only non-dominate solutions
        if self.fullnondominated >= 4:
            self.exitCode = 6
            return

        # generation 0: initialization
        if self.generation == 0:
            print ('------------------------------------\nGeneration 0')
            self.gen_times = []
            self.fullnondominated = 0
            if self.dump_generations is not None:
                self.dumpfile = open('generations.txt', 'w')

            time0 = time()
            self.start_time = time0

            # compute initial energies
            # base class DESolver.__init__()  populates the initial population 
            
            self.population_energies = []
            for trial in self.pop:
                energies = self.score_function(trial)
                if self.dif in '+-':
                    difs = []
                    for i in range(self.nmodels):
                        for j in range(i+1,self.nmodels):
                            res = energies[i] - energies[j]
                            if self.dif == '-':
                                res = -1.0 * res
                            difs.append(res)
                    self.population_energies.append(difs)
                else:
                    self.population_energies.append(energies)

##             print '------- BEGIN CONTROL ----------------------------'
##             np.set_printoptions(precision=14)
##             for i in range(30):
##                 print i, self.pop[i], self.population_energies[i]
##             print '------- END CONTROL ------------------------------'

            timeElapsed = time() - time0
            print ('generation took {} s'.format(timeElapsed))
            self.gen_times.append(timeElapsed)
            if self.dump_generations is not None:
                print >> self.dumpfile, self.generation_string('0')
        
        else: # generation >= 1
            time0 = time()
            print ('------------------------------------')
            print ('\nGeneration {}'.format(self.generation))

            print ("Generating new candidates...")
            self.new_generation_energies = []
            for p in range(self.pop_size):
                # generate new solutions.,reject those out-of-bounds or repeated

                while True:
                    # force a totally new solution
                    ltrial = self.calcTrialSolution(p)
##                     if np.all(ltrial == self.pop[p]):
##                         print 'SOL REPEATED'
##                         continue
                    # check if out of bounds
                    inbounds = np.all(np.logical_and(ltrial <= self.opt_maxs, ltrial >= self.opt_mins))
                    if inbounds: break
                    else: continue
                
                #Handle new solution
                self.new_population[p,:] = ltrial
                
                # compute trialSolution energies
                trialEnergies = self.score_function(ltrial)
                if self.dif in '+-':
                    difs = []
                    for i in range(self.nmodels):
                        for j in range(i+1,self.nmodels):
                            dif = trialEnergies[i] - trialEnergies[j]
                            if self.dif == '-':
                                dif = -1.0 * dif
                            difs.append(dif)
                    self.new_generation_energies.append(difs)
                else:
                    self.new_generation_energies.append(trialEnergies)

            timeElapsed1 = time() - time0
            
            energyComparison = dominance_delta(self.new_generation_energies, 
                                               self.population_energies)

            #working structures for sorting
            objectives = [[]]
            working_sols = [np.array([0.0])]
            
            n_keys = 0
            newBetterSols = 0
            for i in range(len(energyComparison)):
                if energyComparison[i] > 0:
                    objectives.append(self.new_generation_energies[i])
                    working_sols.append(np.copy(self.new_population[i]))
                    n_keys += 1
                    newBetterSols += 1
                elif energyComparison[i] < 0:
                    objectives.append(self.population_energies[i])
                    working_sols.append(np.copy(self.pop[i]))
                    n_keys += 1
                else: #energyComparison[i] == 0
                    objectives.append(self.population_energies[i])
                    working_sols.append(np.copy(self.pop[i]))
                    objectives.append(self.new_generation_energies[i])
                    working_sols.append(np.copy(self.new_population[i]))
                    n_keys += 2
            
            print ('New dominant solutions: {}'.format(newBetterSols))
            print
            sortingtime = time()
            print ("Sorting solutions...")

            # sort solutions by dominance
            sorter = MOOSorter(objectives, indexes = list(range(1,n_keys+1)))
            nondominated_waves = sorter.get_non_dominated_fronts()
            flengths = [len(i) for i in nondominated_waves]
            print ('Front lengths: {}'.format(flengths))
            
            # rebuild current population to pop_size
            self.pop = []
            self.population_energies = []
            
            fronts = []    # holds indexes of (used) non-dominated waves
            
            remaining = []
            for ndf in nondominated_waves:
                excess = len(self.pop) + len(ndf) - self.pop_size
                if len(self.pop) < self.pop_size and excess >= 0:
                    # create a set of solutions to exactly complete pop to pop_size
                    remaining = remove_most_crowded(objectives,
                                                    indexes = ndf,
                                                    knumber = 3, 
                                                    remove_n = excess)
                    break
                elif len(self.pop) == self.pop_size:
                    # do nothing, pop complete
                    break
                else:
                    # just copy  the whole ndf 'front' into pop
                    fronts.append(ndf)
                    for s in ndf:
                        self.pop.append(working_sols[s])
                        self.population_energies.append(objectives[s])
            
            # use the (trimmed)  last front to complete pop to self.pop_size 
            fronts.append(remaining)
            for s in remaining:
                self.pop.append(working_sols[s])
                self.population_energies.append(objectives[s])
            
            timeElapsed2 = time() - sortingtime

            flengths = [len(i) for i in fronts]
            print ('Used front lengths:', flengths, '= %d'%sum(flengths))
            
##             print '------- BEGIN CONTROL ----------------------------'
##             np.set_printoptions(precision=14)
##             for i in fronts[0]:
##                 print i, working_sols[i], objectives[i]
##             print '------- END CONTROL ------------------------------'

            n_nondominated = flengths[0]
            print ('%d non-dominated solutions'%(n_nondominated))
            if n_nondominated == len(self.pop):
                self.fullnondominated += 1
            else:
                self.fullnondominated = 0

            #print 'room for improvement', self.roomForImprovement
            
            if ((not self.trueMetric) and newBetterSols <= self.roomForImprovement and len(fronts) == 1) or (self.trueMetric and self.nmodels == 2 and newBetterSols <= self.roomForImprovement):
                self.generationsWithNoImprovement += 1
            else:
                self.generationsWithNoImprovement = 0
            print ('generations with no improvement:', self.generationsWithNoImprovement)
            
            timeElapsed = time() - time0
            print
            print ("Generation %d finished, took %6.3f s" % (self.generation, timeElapsed))
            print ('%6.3f'% timeElapsed1, 's generating new pop')
            print ('%6.3f'% timeElapsed2, 's in non-dominant sorting')
            self.gen_times.append(timeElapsed)
            if self.dump_generations is not None:
                if self.generation in self.dump_generations:
                    print >> self.dumpfile, self.generation_string(self.generation)

        self.generation += 1
        return
    
    exitCodeStrings = (
    "not done",
    "Solution found by energy criterium",
    "Solution found by diversity criterium",
    "Hit max generations",
    "Too many generations with no improvement",
    "Solution found by convergence criterium",
    "Too many generations with only non-dominant solutions")

    def generation_string(self, generation):
        generation = str(generation)
        res = 'generation %s -------------------------\n'%generation
        for s,o in zip(self.pop, self.population_energies):
            sstr = ' '.join([str(i) for i in s])
            ostr = ' '.join([str(i) for i in o])
            res = res + '%s %s\n'%(sstr, ostr)
        return res
        
    def finalize(self):
        if self.exitCode == 0:
            self.exitCode = -1
        ttime = time() - self.start_time
        print ('=============================================')
        print ("Finished!")
        print (GDE3Solver.exitCodeStrings[self.exitCode])
        print()
        print ('%d generations'%(self.generation))
        print ("Total time: %g s (%s)"% (ttime, utils.s2HMS(ttime)))
        print()
        if self.dump_generations is not None:
            if self.generation-1 not in self.dump_generations:
                print(self.generation_string(self.generation-1), file=self.dumpfile)
            self.dumpfile.close()
