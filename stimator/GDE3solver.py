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

            #working structures for sorting (1 based indexation)
            self.objectives = [[]]
            working_sols = [numpy.array([0.0])]
            self.dom_dict = {}
            
            indx = 1
            newBetterSols = 0
            for i in range(len(energyComparison)):
                if energyComparison[i] > 0:
                    self.objectives.append(self.new_generation_energies[i])
                    self.dom_dict[indx] = []
                    working_sols.append(numpy.copy(self.new_population[i]))
                    indx += 1
                    newBetterSols += 1
                elif energyComparison[i] < 0:
                    self.objectives.append(self.population_energies[i])
                    self.dom_dict[indx] = []
                    working_sols.append(numpy.copy(self.population[i]))
                    indx += 1
                else: #energyComparison[i] == 0
                    self.objectives.append(self.population_energies[i])
                    self.dom_dict[indx] = []
                    working_sols.append(numpy.copy(self.population[i]))
                    indx += 1
                    self.objectives.append(self.new_generation_energies[i])
                    self.dom_dict[indx] = []
                    working_sols.append(numpy.copy(self.new_population[i]))
                    indx += 1
            print 'New dominant solutions: %d' %(newBetterSols)
            print
            sortingtime = time()
            print "Sorting solutions...",

            # sort solutions by dominance
            self.getDominanceTree(self.dom_dict.keys())
            nondominated_waves = self.ndf2list()
            
            # rebuild current population to populationSize
            self.population = []
            self.population_energies = []
            
            fronts = []    # holds indexes of (used) non-dominated waves
            
            full = False
            while not full:
                for ndf in nondominated_waves:
                    tempObjDic = {}
                    if len(self.population) + len(ndf) > self.populationSize and len(self.population) < self.populationSize:
                        # create a set of solutions to exactly complete pop to populationSize
                        for k in ndf:
                            tempObjDic[k] = self.objectives[k]
                        while len(self.population) + len(tempObjDic) != self.populationSize:
                            tempObjDic = removeMostCrowded(tempObjDic, 3)
                        full = True
                        break
                    elif len(self.population) == self.populationSize:
                        # do nothing, pop complete
                        full = True
                        break
                    else:
                        # just copy  another 'wave' of non-dominated solutions into pop, because
                        # len(ndf) + len(pop) <= populationSize
                        fronts.append(ndf)
                        for k in ndf:
                            self.population.append(working_sols[k])
                            self.population_energies.append(self.objectives[k])
            
            
            # use the (trimmed)  tempObjDic to complete pop to self.populationSize 
            # and complete last front
            fronts.append(tempObjDic.keys())
            for i in tempObjDic.keys():
                self.population.append(working_sols[i])
                self.population_energies.append(self.objectives[i])
            
            timeElapsed = time() - sortingtime
            print 'done, took %6.3f'% timeElapsed, 's'

            flengths = [len(i) for i in fronts]
            print 'Front lengths:', flengths
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

    def getDominanceTree(self, nodeList):
        """ This function applies the 'Divide and Conquer' method recursively generating
        the dominance tree (divides) and returning the product of the mergeDominanceTree function (conquers).
        This version doesn't use delayed insertion of dominated node yet."""
        size = len(nodeList)
        if size > 1:
            leftTree  = self.getDominanceTree(nodeList[:size/2])
            rightTree = self.getDominanceTree(nodeList[ size/2:])
            return self.mergeDominanceTrees(leftTree, rightTree)
        else:
            return nodeList

    def mergeDominanceTrees(self, leftTree, rightTree):
        """ This function merges (conquers) the dominance trees recursively (not using delayed insertion yet). """
        leftNode = Node(leftTree[0])
        rightNode = Node(rightTree[0])
        switch = 0 #switch indicates at the end of the loop which conditional of the loop is used. If 'else' is used switch does not change.
        while leftNode.nodeIndex != -1 and rightNode.nodeIndex != -1: #and leftTree != [] and rightTree != []:
            switch = 0
            dominanceTestResult = dominance(self.objectives[rightNode.nodeIndex], self.objectives[leftNode.nodeIndex])
            #rightDelayedInsertionList and leftDelayedInsertionList are [] for now.
            leftDelayedInsertionList = []
            rightDelayedInsertionList = []
            if dominanceTestResult < 0: 
                switch = 1
                if rightDelayedInsertionList != []:
                    mergeDominanceTrees([self.dom_dict[rightNode.nodeIndex][0]], rightDelayedInsertionList)
                tempNode = rightNode.nodeIndex
                rightNode.nextSiblingNode(rightTree) 
                for i in self.dom_dict.keys():
                    if tempNode in self.dom_dict[i]:
                        self.dom_dict[i].remove(tempNode)
                if len(self.dom_dict[leftNode.nodeIndex]) > 0:
                    self.dom_dict[leftNode.nodeIndex] = self.mergeDominanceTrees(self.dom_dict[leftNode.nodeIndex], [tempNode])
                else:
                    self.dom_dict[leftNode.nodeIndex].append(tempNode)
                if tempNode in rightTree:
                    rightTree.remove(tempNode)
                if leftNode.nodeIndex != -1 and rightNode.nodeIndex == -1 and rightTree != []: #New!
                    rightNode.nodeIndex = rightTree[0]
            elif dominanceTestResult > 0:
                switch = 2
                if leftDelayedInsertionList != []:
                    self.mergeDominanceTrees([self.dom_dict[leftNode.nodeIndex][0]], leftDelayedInsertionList)
                tempNode = leftNode.nodeIndex
                leftNode.nextSiblingNode(leftTree) 
                for i in self.dom_dict.keys():
                    if tempNode in self.dom_dict[i]:
                        self.dom_dict[i].remove(tempNode)
                if len(self.dom_dict[rightNode.nodeIndex]) > 0:
                    self.dom_dict[rightNode.nodeIndex] = self.mergeDominanceTrees(self.dom_dict[rightNode.nodeIndex], [tempNode])
                else:
                    self.dom_dict[rightNode.nodeIndex].append(tempNode)
                if tempNode in leftTree:
                    leftTree.remove(tempNode)
                if leftTree == []:
                    leftTree = rightTree
                    rightTree = []
                else:
                    rightNode.nodeIndex = rightTree[0] #This is necessary for the new element from leftTree to be compared with every element os rightTree
            else:
                #No switch assignment
                if rightDelayedInsertionList != []:
                    self.dom_dict[rightNode.nodeIndex] = self.mergeDominanceTrees(self.dom_dict[rightNode.nodeIndex], rightDelayedInsertionList)
                rightNode.nextSiblingNode(rightTree)
                if rightNode.nodeIndex == -1:
                    rightNode.nodeIndex = rightTree[0]
                    leftNode.nextSiblingNode(leftTree)
        if switch == 0:
            leftNode.nodeIndex = leftTree[0]
            rightNode.nodeIndex = rightTree[0]
        while rightTree != [] and rightNode.nodeIndex != -1:
            leftTree.append(rightNode.nodeIndex)
            tempNode = rightNode.nodeIndex
            rightNode.nextSiblingNode(rightTree)
            rightTree.remove(tempNode)
            if rightNode.nodeIndex == -1 and rightTree !=[]:
                leftTree = leftTree + rightTree
        #This block merges cousins, i.e. the nodes at the same non-dominance level which are children of different non-dominant sibling solutions
        #I'm not sure if this should be done every time or just in the final.
        if self.dom_dict != {}:
            firstFrontToBeCompared = [0, 0]
            grandChildren = []
            #This finds and merges cousins which are children of dominated solutions
            for i in leftTree:
                if self.dom_dict[i] != []:
                    grandChildren.append(self.dom_dict[i])
                    if len(grandChildren) == 2:
                        grandChildren = [self.mergeDominanceTrees(grandChildren[0], grandChildren[1])]
            #Continues to merge the cousins which are children of non-dominated parents
            while len(firstFrontToBeCompared) > 1:
                firstFrontToBeCompared = []
                for i in leftTree:
                    nonDominantSolution = True
                    for j in self.dom_dict.keys():
                        if i in self.dom_dict[j]:
                            nonDominantSolution = False
                    if nonDominantSolution == True and self.dom_dict[i]!=[]:
                        firstFrontToBeCompared.append(i)
                    if len(firstFrontToBeCompared) > 1:
                        #Next two lines are necessary to avoid errors inside merge function
                        cousin1 = self.dom_dict[firstFrontToBeCompared[0]] #The 'list' function is necessary for the merging to work
                        cousin2 = self.dom_dict[firstFrontToBeCompared[1]] #If 'list' is not used an instance of a dictionary object passes into the merging function and that raises problems
                        self.mergeDominanceTrees(cousin1, cousin2)
                        break
        return leftTree

    def ndf2list(self):
        """Organizes a list of nondominated fronts from self.dom_dict"""
        nonDominatedFronts = []
        while len(self.dom_dict) > 0:
            nonDominated = []
            for k in self.dom_dict:
                dominated = False
                for j in self.dom_dict.values():
                    if k in j:
                        dominated = True
                        break
                if not dominated:
                    nonDominated.append(k)
            nonDominatedFronts.append(nonDominated)
            for k in nonDominated:
                del self.dom_dict[k]
        return nonDominatedFronts

    #------------------------------------------------------------------------------------------------------------------------------------
    #End of non-dominated sorted algorithm methods


class Node:
    """Node object, which basically contains an index that can be incremented 
    by method nextSiblingNode to iterate over a list of sibling nodes in a tree.
    """

    def __init__(self, nodeIndex):
        """This function creates the 'Node' object."""
        self.nodeIndex = nodeIndex

    def nextSiblingNode(self, tree):
        """ This function increments the index of the 'Node' object and returns the sibling next to the current node."""
        #Each node should have a sibling list of which itself is part and the lists for every subling must be in the same order. This is not done yet.
        if tree.index(self.nodeIndex) < len(tree)-1:
            self.nodeIndex = tree[tree.index(self.nodeIndex)+1]
        else:
            self.nodeIndex = -1


def removeMostCrowded(x, knumber, pop_removed=False, distance_fn=None):

    """
    Removes a point with the smaller crowiding distance measured as
    the sum of the distances of the points to their k-nearest neighbours.
    x is a dictionary which values are the lists of objective values for each solution.
    distance_fn
    is an optional function that takes two points and returns the
    distance between them.  If distance_fn is None (the default), the
    Euclidean distance is used.
    This function was written by extensive modification of the Calculate function
    in the kNN module of the Biopython package.
    """
    if len(x) == 0 :
        return x
    if len(x) < 3:
        x.popitem()
        return x
    else:
        keys = x.keys()
        keys.sort()
        n_objs = len(x[keys[0]])
        extremes = []
        for i in range(n_objs):
            maximum = x[keys[0]][i]
            tempMax = [keys[0]]
            minimum = x[keys[0]][i]
            tempMin = [keys[0]]
            for j in keys[1:]:
                v = x[j][i]
                if v < minimum:
                    tempMin = [j]
                    minimum = v
                elif v == minimum and j not in tempMin:
                    tempMin.append(j)                
                elif v > maximum:
                    tempMax = [j]
                    maximum = v
                elif v == maximum and j not in tempMax:
                    tempMax.append(j)                
            for k in tempMin:
                if k not in extremes:
                    extremes.append(k)
            for k in tempMax:
                if k not in extremes:
                    extremes.append(k)
        for i in keys:
            x[i] = numpy.array(x[i])
        lista =[]
        distanceMatrix = []
        for i in range(len(keys)):
            distances = []
            for j in range(len(keys)):
                if j < i:
                    distances.append(distanceMatrix[j][i])
                elif j ==i:
                    distances.append(0)
                elif j > i:
                    temp = x[keys[i]] - x[keys[j]]
                    dist = numpy.sqrt(numpy.dot(temp,temp))
                    distances.append(dist)
            distanceMatrix.append(numpy.copy(distances))
            distances.sort()
            if knumber < len(distances):
                lista.append(distances[1:knumber+1])
            else:
                lista.append(distances[1:])
        distancesAndKeys = []
        tolook = []
        extremes.sort()
        if keys != extremes:
            for i in range(len(lista)):
                if keys[i] in extremes:
                    distancesAndKeys.append((10**300, keys[i]))
                    tolook.append(('inf', keys[i]))
                else:
                    distancesAndKeys.append(((sum(lista[i])), keys[i]))
                    tolook.append(((sum(lista[i])), keys[i]))
        else:
            for i in range(len(lista)):
                distancesAndKeys.append(((sum(lista[i])), keys[i]))        
        if pop_removed == True:
            print 'distancesAndKeys minimum', x[min(distancesAndKeys)[1]]
        del x[min(distancesAndKeys)[1]]
        for i in x.keys(): #TODO confirm if this object should be returned as list or Arrays
            x[i] = list(x[i]) #dto be used by other functions in next generations
        return x


#________________________________________________________________
#Tests for the non-dominated sorted algorithm methods

if __name__ == "__main__":

    class ndsaTest(GDE3Solver):

        def __init__(self, n_nodes, n_objectives, report = False):
            self.n_nodes = n_nodes
            self.n_objectives = n_objectives
            print 'Test initialized'
            print 'Nodes: %d  Objectives: %d' % (n_nodes, n_objectives)
            numpy.random.seed(2)
            #This must begin in 1.
            self.dom_dict = {}
            self.objectives = {}
            for i_node in range(self.n_nodes):
                self.dom_dict[i_node+1] = []
                self.objectives[i_node+1] = []
                for i_objective in range(self.n_objectives):
                    self.objectives[i_node+1].append(numpy.random.rand())
            keys = self.objectives.keys()
            print 'objectives dict: %d keys with %d elements' % (len(self.objectives), len(self.objectives[1]))
            print 'self.dom_dict and objectiveDic created; entering getDominanceTree'
            self.getDominanceTree(keys)
            print 'Entering nondominated_waves'
            nondominated_waves = self.ndf2list()
            print '%d nondominated_waves'% len(nondominated_waves)
            print '\nTesting non-dominance between solutions in the same front...',
            for k in nondominated_waves:
                if len(k) == 1:
                    d =0
                elif len(k) > 1:
                    p = 0
                    while p < len(k)-1:
                        r = p + 1
                        while r < len(k):
                            d = dominance(self.objectives[k[r]], self.objectives[k[p]])
                            if d != 0:
                                print '\n FAILED:'
                                print '\n\n\nNumber of solutions', self.n_nodes, 'number of objectives', self.n_objectives
                                print 'Domination relationship in front', k, 'between nodes', k[p], 'and', k[r],'. Test not passed.\n\n'
                                break
                            else:
                                r += 1
                        if d != 0:
                            break
                        else:
                            p += 1
                    if d != 0:
                        break
                if d != 0:
                    return
                elif len(k) == 0:
                    print '\n FAILED: empty front found!'
                    return
            print 'passed.'
            if len(nondominated_waves) > 1:
                print '\nTesting dominance relationship between solutions in different fronts...'
                #Solution in rFront must be dominated by at least one solution in pFront and cannot dominate any solution in pFront.
                pFront = 0
                rFront = pFront + 1
                while pFront < len(nondominated_waves)-1:
                    for down in nondominated_waves[rFront]:
                        print 'Assigning totalDominance for the first time in the test'
                        totalDominance = 0
                        for up in nondominated_waves[pFront]:
                            d = dominance(self.objectives[down], self.objectives[up])
                            if d == 1:
                                print '\n\n\nNumber of solutions', self.n_nodes, 'number of objectives', self.n_objectives
                                print 'Solution in front', pFront, ', (', up, ') is dominated by solution in front', rFront, ', (', down, '). Test not passed.\n\n'
                                break
                            totalDominance = totalDominance + d
                        if d == 1 or totalDominance < 1:
                            break
                    if d == 1 or totalDominance < 1:
                        break
                    pFront += 1
                    rFront += 1
                if d == 1 or totalDominance < 1:
                    return #pFront, up, rFront, down
                print 'Non-dominance test between solutions in different fronts passed.'
            if len(nondominated_waves) == 1:
                print 'Only non-dominated solutions - no test between different fronts'

    counter = 0
    for i in range(10, 60):
        for j in range(2, 5):
            print '--------------------------------------------------------'
            print 'test %d'% (counter +1)
            if (i == 0 and j == 0):
                print '**************************'
                ndsaTest(i, j, report = False)
                print '**************************'
            else:
                ndsaTest(i, j)
            counter += 1
    print 'Tests finished successfully!'
    
