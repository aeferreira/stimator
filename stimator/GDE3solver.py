# Project S-timator

from time import time
import numpy, timecourse

from de import DESolver
import dynamics
import utils
import moosorting

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
##     if numpy.all(d): return -size
##     d = vec2 >= vec1
##     if numpy.all(d): return size
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

def energies_dominance_delta(old_energies, new_energies):
    """Compute dominance relationship between solutions of two generations."""
    return [dominance(o,n) for o,n in zip(old_energies, new_energies)]


class ModelSolver(object):
    
    """Driver to compute timecourses, given a model and a trial vector.
    
    'model' is a S-timator model object;
    't0' and 'tf' are floats;
    'npoints' is a positive integer.
    'observerd' is a list of variable names to output.
    'optnames'are variable names. Their initial values will be optimized
    in the experimental design."""
    
    def __init__(self, model, t0, npoints, tf, optnames, observed):
        self.t0 = t0
        self.npoints = npoints
        self.tf = tf
        
        self.vector = numpy.copy(dynamics.state2array(model,"init"))

        self.modelvarnames = model().varnames
        
        self.model = model
        self.optvars = listify(optnames)
        self.observed = listify(observed)
        
        #check that the names exist as variables in the model
        for name in self.optvars+self.observed:
            if not (name in self.modelvarnames):
                raise AttributeError('%s is not a variable in model'%name)
        
  
        self.optvars_indexes = numpy.array([self.modelvarnames.index(name) for name in self.optvars])
        self.obsvars_indexes = numpy.array([self.modelvarnames.index(name) for name in self.observed])
            
    def solve(self, trial):
        """Returns the solution for self.model for a particular trial of initial values."""
        self.trial = trial # trial is a list of values
                    
        for value,indx in zip(self.trial, self.optvars_indexes):
            self.vector[indx] = value

        return dynamics.solve(self.model, tf = self.tf, npoints = self.npoints, t0 = self.t0, 
                                          initial = self.vector).copy(self.observed)
 
# helper to transform string arguments in lists:
def listify(arguments):
    if isinstance(arguments, list) or isinstance(arguments, tuple):  return [a.strip() for a in arguments]
    if isinstance(arguments, str) or isinstance(arguments, unicode): 
        arguments = arguments.split()
        return [a.strip() for a in arguments]


class GDE3Solver(DESolver):

    """ 
    Adaptation of DESolver for multiobjective optimization.
    Based on 
    Kukkonen, S. and Lampinen, J. (2005) GDE3: 
    The third step of generalized differential evolution, 
    2005 IEEE Congress on Evolutionary Computation. 
    """

    def __init__(self, models, toOpt, objectiveFunction, observed, npoints, t0, tf, populationSize, maxGenerations, deStrategy, diffScale, crossoverProb, 
                 cutoffEnergy, 
                 useClassRandomNumberMethods, 
                 dif = '0', dump_generations = None):
        
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
        self.opt_maxs = numpy.array(self.opt_maxs)
        self.opt_mins = numpy.array(self.opt_mins)
        
        #self.toOptKeys = self.toOpt.keys()
        self.objFunc = objectiveFunction
        self.dump_generations = dump_generations
        
        DESolver.__init__(self, 
                          len(toOpt), populationSize, maxGenerations, 
                          self.opt_mins, self.opt_maxs, 
                          deStrategy, diffScale, crossoverProb, 
                          cutoffEnergy, useClassRandomNumberMethods)

        self.deltaT = (tf - t0)/npoints
                
        self.objFuncList = []

        for m in self.models:
            self.objFuncList.append(ModelSolver(m, self.t0, self.npoints, self.tf, 
                                                self.toOptKeys, self.observed))
        
        self.dif = dif
        
        # threshold for improvement based on 5 % of new solution count
        self.roomForImprovement = int(round(0.05 * populationSize))
        
        #counter of number of generations with only  non-dominated solutions
        self.fullnondominated = 0
        
        str2distance = {'extKL'   :timecourse.extendedKLdivergence,
                        'kremling':timecourse.kremling,
                        'KL'      :timecourse.KLdivergence,
                        'L2'      :timecourse.L2}
        if self.objFunc not in str2distance.keys():
            raise ("%s is not an implemented divergence function" % self.objFunc)
        
        if self.objFunc in ('kremling','L2'):  #symetric measures
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
        self.population_energies = []
        self.new_generation_energies = [[] for i in range(self.populationSize)]
        self.new_population = numpy.empty((self.populationSize,len(self.toOpt)))
        self.gen_times = []

    def EnergyFunction(self, trial):
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
        if self.generation >= self.maxGenerations:
            self.exitCode = 3
            return
        # Hit several generations with only non-dominate solutions
        if self.fullnondominated >= 4:
            self.exitCode = 6
            return
        
        ## # no need to try another generation if we are done (energy criterium)
        ## if self.atSolution:
            ## self.exitCode = 1
            ## return
        # generation 0: initialization
        if self.generation == 0:
            print '------------------------------------\nGeneration 0'
            self.gen_times = []
            self.fullnondominated = 0
            if self.dump_generations is not None:
                self.dumpfile = open('generations.txt', 'w')
            
            time0 = time()
            self.elapsed = time0

            # compute initial energies
            # base class DESolver.__init__()  populates the initial population 
            
            self.population_energies = []
            for trial in self.population:
                energies = self.EnergyFunction(trial)
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
            print self.firstpop_string(20, 'after regeneration 0 before dump')

            timeElapsed = time() - time0
            print 'generation took', timeElapsed, 's'
            self.gen_times.append(timeElapsed)
            if self.dump_generations is not None:
                print >> self.dumpfile, self.generation_string('0')
            print self.firstpop_string(20, 'after regeneration 0 after dump')
        
        else: # generation >= 1
            time0 = time()
            print '------------------------------------\nGeneration %d'% self.generation

            already_aslists = [list(i) for i in self.population]
            
            print "Generating new candidates...",
            for p in range(self.populationSize):
                # generate new solutions.,reject those out-of-bounds or repeated
                insideBoundaries = False
                while insideBoundaries == False:
                    # force a totally new solution
                    self.calcTrialSolution(p)
                    ltrial = list(self.trialSolution)
                    if ltrial not in already_aslists:
                        already_aslists.append(ltrial)
                    else:
                        continue
                    # check if out of bounds
                    inbounds = numpy.all(numpy.logical_and(self.trialSolution <= self.opt_maxs, self.trialSolution >= self.opt_mins))
                    if inbounds: break
                    else: continue
                
                self.new_population[p,:] = self.trialSolution
                
                # compute trialSolution energies
                trialEnergies = self.EnergyFunction(self.trialSolution)
                if self.dif in '+-':
                    difs = []
                    for i in range(self.nmodels):
                        for j in range(i+1,self.nmodels):
                            dif = trialEnergies[i] - trialEnergies[j]
                            if self.dif == '-':
                                dif = -1.0 * dif
                            difs.append(dif)
                    self.new_generation_energies[p] = difs
                else:
                    self.new_generation_energies[p] = trialEnergies

            timeElapsed = time() - time0
            print 'done, took %6.3f'% timeElapsed, 's'
            
            energyComparison = energies_dominance_delta(self.new_generation_energies, self.population_energies)
            #print 'Dominance comparison with previous generation:'
            #print energyComparison

            #working structures for sorting
            self.objectives = {}
            working_sols = [numpy.array([0.0])]
            
            n_keys = 1
            newBetterSols = 0
            for i in range(len(energyComparison)):
                if energyComparison[i] > 0:
                    self.objectives[n_keys] = self.new_generation_energies[i]
                    working_sols.append(numpy.copy(self.new_population[i]))
                    n_keys += 1
                    newBetterSols += 1
                elif energyComparison[i] < 0:
                    self.objectives[n_keys] = self.population_energies[i]
                    working_sols.append(numpy.copy(self.population[i]))
                    n_keys += 1
                else: #energyComparison[i] == 0
                    self.objectives[n_keys] = self.new_generation_energies[i]
                    working_sols.append(numpy.copy(self.new_population[i]))
                    n_keys += 1
                    self.objectives[n_keys] = self.population_energies[i]
                    working_sols.append(numpy.copy(self.population[i]))
                    n_keys += 1
            
            print '%%%%%%%%%%%%%%%%% ini OBJS'
            for k in self.objectives:
                print k, ' ---->', self.objectives[k]
            print '%%%%%%%%%%%%%%%%% end OBJS'
            
            print 'New dominant solutions: %d' %(newBetterSols)
            print
            sortingtime = time()
            print "Sorting solutions..."

            # sort solutions by dominance
            sorter = moosorting.MOOSorter(self.objectives)
            
            nondominated_waves = sorter.get_non_dominated_fronts()

            flengths = [len(i) for i in nondominated_waves]
            print 'Waves lengths:', flengths
            
            # rebuild current population to populationSize
            self.population = []
            self.population_energies = []
            
            fronts = []    # holds indexes of (used) non-dominated waves
            
            for ndf in nondominated_waves:
                tempObjDic = {}
                excess = len(self.population) + len(ndf) - self.populationSize
                if len(self.population) < self.populationSize and excess >= 0:
                    # create a set of solutions to exactly complete pop to populationSize
                    for k in ndf:
                        tempObjDic[k] = self.objectives[k]
                    tempObjDic = removeMostCrowded(tempObjDic, 3, excess)
                    break
                elif len(self.population) == self.populationSize:
                    # do nothing, pop complete
                    break
                else:
                    # just copy  the whole 'wave' of non-dominated solutions into pop
                    fronts.append(ndf)
                    for k in ndf:
                        self.population.append(working_sols[k])
                        self.population_energies.append(self.objectives[k])
            
            
            # use the (trimmed)  tempObjDic to complete pop to self.populationSize 
            # and record last front
            fronts.append(tempObjDic.keys())
            for i in tempObjDic.keys():
                self.population.append(working_sols[i])
                self.population_energies.append(self.objectives[i])
            
            print self.firstpop_string(20, 'after regeneration')
            
            timeElapsed = time() - sortingtime
            print 'done, took %6.3f'% timeElapsed, 's'

            flengths = [len(i) for i in fronts]
            print 'Front lengths:', flengths, ' = ', sum(flengths)
            #nondominated_indxs = nondominated_solutions(self.population_energies)
            #n_nondominated = len(nondominated_indxs)
            n_nondominated = flengths[0]
            print '%d non-dominated solutions'%(n_nondominated)
            if n_nondominated == len(self.population):
                self.fullnondominated += 1
            else:
                self.fullnondominated = 0

            #print 'room for improvement', self.roomForImprovement
            
            if ((not self.trueMetric) and newBetterSols <= self.roomForImprovement and len(fronts) == 1) or (self.trueMetric and self.nmodels == 2 and newBetterSols <= self.roomForImprovement):
                self.generationsWithNoImprovement += 1
            else:
                self.generationsWithNoImprovement = 0
            timeElapsed = time() - time0
            print
            print "Generation %d finished, took %6.3f s" % (self.generation, timeElapsed)
            print 'generations with no improvement:', self.generationsWithNoImprovement
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

    def firstpop_string(self, nelems, msg = ""):
        res = 'CONTROL %s -------------------------\n'%msg
        for s,o in zip(self.population[:nelems], self.population_energies[:nelems]):
            sstr = ' '.join([str(i) for i in s])
            ostr = ' '.join([str(i) for i in o])
            res = res + '%s %s\n'%(sstr, ostr)
        return res

    def generation_string(self, generation):
        generation = str(generation)
        res = 'generation %s -------------------------\n'%generation
        for s,o in zip(self.population, self.population_energies):
            sstr = ' '.join([str(i) for i in s])
            ostr = ' '.join([str(i) for i in o])
            res = res + '%s %s\n'%(sstr, ostr)
        return res
        
    def finalize(self):
        if self.exitCode == 0:
            self.exitCode = -1
        ttime = self.elapsed = time() - self.elapsed
        print '============================================='
        print "Finished!"
        print GDE3Solver.exitCodeStrings[self.exitCode]
        print
        print '%d generations'%(self.generation)
        print "Total time: %g s (%s)"% (ttime, utils.s2HMS(ttime))
        print
        if self.dump_generations is not None:
            if self.generation-1 not in self.dump_generations:
                print >> self.dumpfile, self.generation_string(self.generation-1)
            self.dumpfile.close()

def removeMostCrowded(x, knumber = 3, remove_n = 1, verbose = False):

    """
    Removes a point with the smaller crowiding distance measured as
    the sum of the distances of the points to their k-nearest neighbours.
    x is a dictionary which values are the point coordinates as a list.
    """
    
    n_points = len(x)
    if n_points == 0 :
        return x
    if remove_n < 1:
        return x
    keys = list(x.keys())
    n_objs = len(x[keys[0]])
    
    points = range(n_points) #holds current indexes of points still present
    
    #compute matrix of objectives
    obj_matrix = numpy.empty((n_points, n_objs))
    for i,k in enumerate(keys):
        obj_matrix[i, :] = x[k]

    # find indexes of extreme points (max and min in each dimension)
    
    maxima = numpy.amax(obj_matrix, axis=0)
    minima = numpy.amin(obj_matrix, axis=0)
    
    extremes = []
    for i in range(n_objs):
        dimextremes = []
        for j in points:
            v = obj_matrix[j][i]
            if v == minima[i] or v == maxima[i]:
                dimextremes.append(j)                
        for j in dimextremes:
            if j not in extremes:
                extremes.append(j)
    extremes.sort()
    
    if verbose:
        print 'extreme points', [keys[i] for i in extremes]

    #compute distances
    # TODO: use numpy function
    distanceMatrix = []
    
    for i in range(n_points):
        distances = []
        for j in range(n_points):
            if j < i:
                distances.append(distanceMatrix[j][i])
            elif j ==i:
                distances.append(0.0)
            elif j > i:
                temp = numpy.array(obj_matrix[i]) - numpy.array(obj_matrix[j])
                d = numpy.sqrt(numpy.dot(temp,temp))
                distances.append(d)
        distanceMatrix.append(distances)
    
    #sort distances for each point
    distances = []
    for i in points:
        dd = zip(list(distanceMatrix[i]), range(n_points))
        dd.sort()
        distances.append(dd)
    
    #compute k shortest (note: position 0 after sorting is always 0.0)
    last_removed = None
    
    # loop to remove each point
    for n in range(remove_n):
        if len(points) == 0:
            if verbose:
                print '\nNo more points to remove'
            return x

        if verbose:
            print '---------------------------\nRemoving point #%d' % (n+1), ':'
            print '\nremaining points' #, points
            print [keys[k] for k in points]
        
        
        if last_removed is not None:
            #remove reference to last removed point in list of distances
            for i in points:
                indx_last_remove = -1
                d = distances[i]
                for j in range(len(d)):
                    if d[j][1] == last_removed:
                        indx_last_remove = j
                        break
                del(distances[i][indx_last_remove])
        
        if verbose:
            print '\nk-distances'
            for i in points:
                dd = distances[i][1:knumber+1]
                print keys[i],
                print ', '.join([ "(%-5.3g to %s)"% (t1, keys[t2]) for (t1,t2) in dd])
            print
        
        ksums = [sum([d[0] for d in distances[i][1:knumber+1]]) for i in points]
        
        #find and remove most crowded
        distancesAndKeys = []
        
        allextremes = True
        for i in points:
            if i not in extremes:
                allextremes = False
                break

        if not allextremes:
            for i,k in enumerate(points):
                if k in extremes:
                    distancesAndKeys.append((10**300, k))
                else:
                    distancesAndKeys.append((ksums[i], k))
        else:
            for i,k in enumerate(points):
                distancesAndKeys.append((ksums[i], k))      
        
        last_removed = min(distancesAndKeys)[1]
        points.remove(last_removed)
        mc_key = keys[last_removed]
        
        if verbose:
            mcv = x[mc_key]
            print 'Point to remove:', mc_key, mcv
        
        del x[mc_key]
    
    return x
