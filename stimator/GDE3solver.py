# Project S-timator

from time import time
import numpy, timecourse

from de import DESolver
import dynamics
import utils

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
        
        self.dominance_rdepth = 0
        
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
            timeElapsed = time() - time0
            print 'generation took', timeElapsed, 's'
            self.gen_times.append(timeElapsed)
            if self.dump_generations is not None:
                print >> self.dumpfile, self.generation_string('0')
        
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
            self.objectives = [[]]
            working_sols = [numpy.array([0.0])]
            
            n_keys = 0
            newBetterSols = 0
            for i in range(len(energyComparison)):
                if energyComparison[i] > 0:
                    self.objectives.append(self.new_generation_energies[i])
                    working_sols.append(numpy.copy(self.new_population[i]))
                    n_keys += 1
                    newBetterSols += 1
                elif energyComparison[i] < 0:
                    self.objectives.append(self.population_energies[i])
                    working_sols.append(numpy.copy(self.population[i]))
                    n_keys += 1
                else: #energyComparison[i] == 0
                    self.objectives.append(self.population_energies[i])
                    working_sols.append(numpy.copy(self.population[i]))
                    self.objectives.append(self.new_generation_energies[i])
                    working_sols.append(numpy.copy(self.new_population[i]))
                    n_keys += 2
            self.dom_dict = {}
            for indx in range(1,n_keys+1): # (1 based indexation)
                self.dom_dict[indx] = []
            
            print 'New dominant solutions: %d' %(newBetterSols)
            print
            sortingtime = time()
            print "Sorting solutions..."

            # sort solutions by dominance
            self.getDominanceTree(self.dom_dict.keys())
            nondominated_waves = self.ndf2list()
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

        
        #self.reportFinal()

    #------------------------------------------------------------------------------------------------------------------------------------
    #This code is an adaptation of the non-dominated sorting algorithm with delayed insertion published in
    #Fang et al (2008) An Efficient Non-dominated Sorting Method for Evolutionary Algorithms, Evolutionary Computation 16(3):355-384

    def depthspaces(self):
        return ' '*(self.dominance_rdepth-1)*4
    
    def getDominanceTree(self, nodeList, verbose = False):
        """ This function applies the 'Divide and Conquer' method recursively generating
        the dominance tree (divides) and returning the product of the mergeDominanceTree function (conquers).
        This version doesn't use delayed insertion of dominated node yet."""
        self.dominance_rdepth += 1
        if verbose:
            spcs = self.depthspaces()
            print spcs, '---------------------------------'
            print spcs, 'ENTERING getDominanceTree()'
            print spcs, 'nodes', nodeList
        size = len(nodeList)
        if size > 1:
            leftTree  = self.getDominanceTree(nodeList[:size/2], verbose)
            rightTree = self.getDominanceTree(nodeList[ size/2:], verbose)
            if verbose:
                print spcs, '@@  left tree', leftTree
                print spcs, '@@ right tree', rightTree
                print spcs, '@@ dom dict', self.dom_dict
            res = self.merge_Dom_trees(leftTree, rightTree)
            if verbose:
                print spcs, '--> return from merge', res
                print spcs, '--> dom dict', self.dom_dict
                print spcs, '**************************'
            self.dominance_rdepth -= 1
            return res
        else:
            if verbose:
                print spcs, '.single node return', nodeList
                print spcs, '**************************'
            self.dominance_rdepth -= 1
            return nodeList

    def merge_Dom_trees(self, left_tree, right_tree):
        """ This function merges (conquers) the dominance trees recursively (not using delayed insertion yet). """
        if len(left_tree) == 0:
            if len(right_tree) < 2:
                return right_tree
        left = Node(left_tree[0])
        right = Node(right_tree[0])
        one_dominates = False 
        while left.valid() and right.valid(): #and left_tree != [] and right_tree != []:
            one_dominates = False
            d = dominance(self.objectives[right.indx], self.objectives[left.indx])
            #rightDelayedInsertionList and leftDelayedInsertionList are [] for now.
            leftDelayedInsertionList = []
            rightDelayedInsertionList = []
            if d < 0: 
                one_dominates = True
                if rightDelayedInsertionList != []:
                    merge_Dom_trees([self.dom_dict[right.indx][0]], rightDelayedInsertionList)
                node = right.indx
                right.move_and_remove(right_tree) 
                for k in self.dom_dict:
                    if node in self.dom_dict[k]:
                        self.dom_dict[k].remove(node)
                self.dom_dict[left.indx] = self.merge_Dom_trees(self.dom_dict[left.indx], [node])
                if left.valid() and not right.valid() and right_tree != []: #New!
                    right.indx = right_tree[0]
            elif d > 0:
                one_dominates = True
                if leftDelayedInsertionList != []:
                    self.merge_Dom_trees([self.dom_dict[left.indx][0]], leftDelayedInsertionList)
                node = left.indx
                left.move_and_remove(left_tree) 
                for k in self.dom_dict:
                    if node in self.dom_dict[k]:
                        self.dom_dict[k].remove(node)
                self.dom_dict[right.indx] = self.merge_Dom_trees(self.dom_dict[right.indx], [node])
                if left_tree == []:
                    left_tree = right_tree
                    right_tree = []
                else:
                    right.indx = right_tree[0] #This is necessary for the new element from left_tree to be compared with every element os right_tree
            else:
                if rightDelayedInsertionList != []:
                    self.dom_dict[right.indx] = self.merge_Dom_trees(self.dom_dict[right.indx], rightDelayedInsertionList)
                right.point_to_next_sibling(right_tree)
                if not right.valid():
                    right.indx = right_tree[0]
                    left.point_to_next_sibling(left_tree)
        if not one_dominates:
            left.indx = left_tree[0]
            right.indx = right_tree[0]
        while right_tree != [] and right.valid():
            left_tree.append(right.indx)
            right.move_and_remove(right_tree)
            if not right.valid() and len(right_tree) > 0:
                left_tree.extend(right_tree)
        #This block merges cousins, i.e. the nodes at the same non-dominance level which are children of different non-dominant sibling solutions
        #I'm not sure if this should be done every time or just in the final.
        if self.dom_dict != {}:
            firstFrontToBeCompared = [0, 0]
            grandChildren = []
            #This finds and merges cousins which are children of dominated solutions
            for i in left_tree:
                if self.dom_dict[i] != []:
                    grandChildren.append(self.dom_dict[i])
                    if len(grandChildren) == 2:
                        grandChildren = [self.merge_Dom_trees(grandChildren[0], grandChildren[1])]
            #Continues to merge the cousins which are children of non-dominated parents
            while len(firstFrontToBeCompared) > 1:
                firstFrontToBeCompared = []
                for i in left_tree:
                    nonDominantSolution = True
                    for k in self.dom_dict:
                        if i in self.dom_dict[k]:
                            nonDominantSolution = False
                    if nonDominantSolution and self.dom_dict[i]!=[]:
                        firstFrontToBeCompared.append(i)
                    if len(firstFrontToBeCompared) > 1:
                        #Next two lines are necessary to avoid errors inside merge function
                        cousin1 = self.dom_dict[firstFrontToBeCompared[0]]
                        cousin2 = self.dom_dict[firstFrontToBeCompared[1]]
                        self.merge_Dom_trees(cousin1, cousin2)
                        break
        return left_tree

    def ndf2list(self):
        """Organizes a list of nondominated fronts from self.dom_dict"""
        nonDominatedFronts = []
        while len(self.dom_dict) > 0:
            allvls = []
            for v in self.dom_dict.values():
                allvls.extend(v)
            nonDominated = [k for k in self.dom_dict if k not in allvls]
            nonDominatedFronts.append(nonDominated)
            for k in nonDominated:
                del self.dom_dict[k]
        return nonDominatedFronts

    #------------------------------------------------------------------------------------------------------------------------------------
    #End of non-dominated sorted algorithm methods


class Node:
    """Node object, which basically contains an index that can be incremented 
    by method point_to_next_sibling to iterate over a list of sibling nodes in a tree.
    """

    def __init__(self, nodeIndex):
        """Initializes Node object with an index"""
        self.indx = nodeIndex

    def point_to_next_sibling(self, tree):
        """Increments the index to the next sibling of the node, if not past the end of the 'tree' list"""
        indx = tree.index(self.indx)
        if indx < len(tree)-1:
            self.indx = tree[indx+1]
        else:
            self.indx = -1
    
    def valid(self):
        return self.indx != -1
        
    def move_and_remove(self, tree):
        node = self.indx
        self.point_to_next_sibling(tree)
        tree.remove(node)
        
        

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

#________________________________________________________________
#Tests for the non-dominated sorted algorithm methods


if __name__ == "__main__":

    class FangEtAL_test(GDE3Solver):
        def __init__(self, report = False):
            data = [[182.08, 100.13, 192.21],
[187.53, 246.16, 203.2],
[197.15, 201.57, 318.86],
[47.48, 74.96, 22.69],
[37.05, 304.83, 381.19],
[126.88, 54.58, 144.17],
[101.77, 49.18, 111.91],
[37.47, 18.63, 446.57]]
            self.dominance_rdepth = 0
            n_nodes = len(data)
            self.n_objectives = len(data[0])
            print 'Test initialized'
            print 'Nodes: %d  Objectives: %d' % (n_nodes, self.n_objectives)
            #dict keys are integers that begin at 1.
            self.dom_dict = {}
            self.objectives = {}
            for i_node in range(n_nodes):
                self.dom_dict[i_node+1] = []
                self.objectives[i_node+1] = data[i_node]
            keys = self.objectives.keys()

            print 'objectives dict: %d keys with %d elements' % (len(self.objectives), len(self.objectives[1]))
            print 'self.dom_dict and objectiveDic created.'
            print '\nComputing nondominated_waves...'
            self.getDominanceTree(keys, verbose = True)
            print 'DOMINANCE DICT'
            for k in self.dom_dict:
                print k, '>', self.dom_dict[k]
            nondominated_waves = self.ndf2list()
            print '%d nondominated_waves:'% len(nondominated_waves)
            for wave in nondominated_waves:
                print wave
            print '\nTesting non-dominance between solutions in the same front...',
            for wave in nondominated_waves:
                if len(wave) == 0:
                    print '\n FAILED: empty front found!'
                    return
                if len(wave) == 1:
                    pass
                else:
                    for p in range(len(wave)-1):
                        for r in range(p+1, len(wave)):
                            d = dominance(self.objectives[wave[r]], self.objectives[wave[p]])
                            if d != 0:
                                print '\n FAILED:'
                                print 'Domination relationship in front', wave, 'between nodes', wave[p], 'and', wave[r],'. Test not passed.\n\n'
                                return
            print 'passed.'
            if len(nondominated_waves) == 1:
                print 'Only non-dominated solutions - no test between different fronts'
                return
            else:
                print 'Testing dominance relationship between solutions in different fronts...',
                #Solution in rFront must be dominated by at least one solution in pFront and cannot dominate any solution in pFront.
                for pFront in range(len(nondominated_waves)-1):
                    for rFront in range(pFront+1, len(nondominated_waves)):
                        for down in nondominated_waves[rFront]:
                            for up in nondominated_waves[pFront]:
                                d = dominance(self.objectives[down], self.objectives[up])
                                if d == 1:
                                    print '\n FAILED:'
                                    print 'Solution in front', pFront, ', (', up, ') is dominated by solution in front', rFront, ', (', down, '). Test not passed.\n\n'
                                    return
                print 'passed.'

    class ndsaTest(GDE3Solver):

        def __init__(self, n_nodes, n_objectives, report = False):
            self.dominance_rdepth = 0
            self.n_nodes = n_nodes
            self.n_objectives = n_objectives
            print 'Test initialized'
            print 'Nodes: %d  Objectives: %d' % (n_nodes, n_objectives)
            numpy.random.seed(2)
            #dict keys are integers that begin at 1.
            self.dom_dict = {}
            self.objectives = {}
            for i_node in range(self.n_nodes):
                self.dom_dict[i_node+1] = []
                self.objectives[i_node+1] = []
                for i_objective in range(self.n_objectives):
                    self.objectives[i_node+1].append(numpy.random.rand())
            keys = self.objectives.keys()

            print 'objectives dict: %d keys with %d elements' % (len(self.objectives), len(self.objectives[1]))
            print 'self.dom_dict and objectiveDic created.'
            print '\nComputing nondominated_waves...'
            self.getDominanceTree(keys, verbose = True)
            print 'DOMINANCE DICT'
            for k in self.dom_dict:
                print k, '>', self.dom_dict[k]
            nondominated_waves = self.ndf2list()
            print '%d nondominated_waves:'% len(nondominated_waves)
            for wave in nondominated_waves:
                print wave
            print '\nTesting non-dominance between solutions in the same front...',
            for wave in nondominated_waves:
                if len(wave) == 0:
                    print '\n FAILED: empty front found!'
                    return
                if len(wave) == 1:
                    pass
                else:
                    for p in range(len(wave)-1):
                        for r in range(p+1, len(wave)):
                            d = dominance(self.objectives[wave[r]], self.objectives[wave[p]])
                            if d != 0:
                                print '\n FAILED:'
                                print 'Domination relationship in front', wave, 'between nodes', wave[p], 'and', wave[r],'. Test not passed.\n\n'
                                return
            print 'passed.'
            if len(nondominated_waves) == 1:
                print 'Only non-dominated solutions - no test between different fronts'
                return
            else:
                print 'Testing dominance relationship between solutions in different fronts...',
                #Solution in rFront must be dominated by at least one solution in pFront and cannot dominate any solution in pFront.
                for pFront in range(len(nondominated_waves)-1):
                    for rFront in range(pFront+1, len(nondominated_waves)):
                        for down in nondominated_waves[rFront]:
                            for up in nondominated_waves[pFront]:
                                d = dominance(self.objectives[down], self.objectives[up])
                                if d == 1:
                                    print '\n FAILED:'
                                    print 'Solution in front', pFront, ', (', up, ') is dominated by solution in front', rFront, ', (', down, '). Test not passed.\n\n'
                                    return
                print 'passed.'
    
##     i = 50
##     j = 2
##     ndsaTest(i, j)
    FangEtAL_test()
    
