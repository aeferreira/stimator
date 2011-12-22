# -*- coding: ISO-8859-1 -*-


#    PythonEquations is a collection of equations expressed as Python classes
#    Copyright (C) 2007 James R. Phillips
#    2548 Vera Cruz Drive
#    Birmingham, AL 35235 USA
#    email: zunzun@zunzun.com
#
#    Version info: $  $
#
#    License: BSD-style (see LICENSE.txt in main source directory)


import numpy, tcdif
from de import DESolver
from numpy import transpose, array, float64, zeros, empty
from time import time
import random


def dominanceComparison(energyListOld, energyListNew):
    """ This function tests the dominance relationship the solutions of the new and the previous generations."""
    dominanceList = []
    for i in range(len(energyListOld)):
        dominanceTestResult = 0
        for j in range(len(energyListOld[i])):
            d = energyListNew[i][j] - energyListOld[i][j]
            if d <= 0 and dominanceTestResult <=0:
                dominanceTestResult -= 1
            elif d >= 0 and dominanceTestResult >=0:
                dominanceTestResult += 1
            else:
                dominanceTestResult = 0
                break
        dominanceList.append(dominanceTestResult)
    return dominanceList

class Node:
    """ This class creates a 'Node' object, which basically contains an index that can be incremented 
    by method nextSiblingNode to iterate over a list of sibling nodes in a tree. """

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

class GDE3Solver(DESolver):

    """ This class is an adaptation of DESolver for multiobjective optimizations.
        This code is based on 
        Kukkonen, S. and Lampinen, J. (2005) GDE3: The third step of generalized differential evolution, The 2005 IEEE Congress on Evolutionary Computation. """

    def __init__(self, models, absoluteMeasurementError, toOpt, 
                 objectiveFunction, 
                 observed, 
                 npoints, t0, tf, 
                 populationSize, 
                 maxGenerations, 
                 biasedCurveNumber, 
                 biasStandardDeviation, 
                 deStrategy, 
                 diffScale, 
                 crossoverProb, 
                 cutoffEnergy, 
                 useClassRandomNumberMethods, 
                 dif = None,
                 keep_track = False):
        
        self.keep_track = keep_track
        self.models = models
        self.npoints = npoints
        self.t0 = t0
        self.tf = tf
        self.observed = observed

        self.nmodels = len(models)
        
        self.toOpt = toOpt
        self.toOptKeys = self.toOpt.keys()
        self.objFunc = objectiveFunction

        self.toOptBounding = [toOpt[i] for i in self.toOptKeys]
        self.toOptBounding = numpy.transpose(self.toOptBounding)
        maxs = numpy.array(self.toOptBounding)[1]
        mins = numpy.array(self.toOptBounding)[0]

        DESolver.__init__(self, 
                          len(toOpt), populationSize, maxGenerations, 
                          mins, maxs, 
                          deStrategy, 
                          diffScale, 
                          crossoverProb, 
                          cutoffEnergy, 
                          useClassRandomNumberMethods)

        self.deltaT = (tf - t0)/npoints
        
        if self.objFunc in ('AIC', 'AICc', 'AICu', 'criterionA', 'modCriterionA', 'criterionD', 'criterionE', 'modCriterionE'):
            #Initialize bias for error simulation
            self.bias = tcdif.generateBias(biasedCurveNumber, len(self.toOptKeys), biasStandardDeviation)
            self.measurementErrors = tcdif.generateUniformMeasurementError(absoluteMeasurementError, npoints)
        else:
            self.bias = None
            self.measurementErrors = None
            
        
        self.objFuncList = []

        for m in self.models:
            self.objFuncList.append(tcdif.Objective(m, self.t0, self.npoints, self.tf, self.objFunc, self.toOptKeys, self.observed, self.bias, self.measurementErrors))
        
        self.dif = dif
        
        # threshold for improvement based on 5 % of new solution count
        self.roomForImprovement = int(round(0.05 * populationSize))
        
        if self.objFunc in ('kremling','L2'):  #symetric measures
            self.trueMetric = True
        else:
            self.trueMetric = False

        str2distance = {'KL'      :tcdif.KLDiscrepancies,
                        'kremling':tcdif.kremling,
                        'KLs'     :tcdif.KLs,
                        'L2'      :tcdif.L2}
        self.distance_func = str2distance.get(self.objFunc, None)
        if self.distance_func is not None:
            self.model_indexes = []
            for i in range(self.nmodels-1):
                for j in range(i+1, self.nmodels):
                    self.model_indexes.append((i,j))
            if self.objFunc in ['KLs', 'KL']:
                self.model_indexes.extend([(j,i) for (i,j) in self.model_indexes])
        
        # working storage arrays
        self.newGenerationEnergyList = [[] for i in range(self.populationSize)]
        
        # storage arrays
        self.fronts = [[]]    # holds fronts created in current generation
        self.frontObj = [[]]  # holds objectives for current generation
        self.ftimes = []

    def EnergyFunction(self, trial):
        trialDic = {}
        for i in range(len(self.toOptKeys)):
            trialDic[self.toOptKeys[i]] = trial[i]
        energies = [f(trialDic) for f in self.objFuncList]
        
        if self.distance_func is not None:
            energies = self.distance_func(energies, self.deltaT, self.model_indexes)
        return energies, False

    def computeGeneration(self):
        # TODO: parallelization here
        # TODO: this is for performance on non-parallelized hardware
        if self.generationsWithNoImprovement > 20:
            self.exitCode = 4
            return
        # Hit max generations
        if self.generation >= self.maxGenerations:
            self.exitCode = 3
            return
        # no need to try another generation if we are done (energy criterium)
        if self.atSolution:
            self.exitCode = 1
            return
        if self.generation == 0:      #compute energies for generation 0
            print '\n\nComputing generation 0.\n'
            self.ftimes = []
            if self.keep_track:
                self.frontsfile = open('fronts.txt', 'w')
                self.objsfile = open('objectives.txt', 'w')
            
            time0 = time()

            self.generationEnergyList = []
            for candidate in range(self.populationSize):
                trialEnergies, self.atSolution = self.EnergyFunction(self.population[candidate]) #This argument must include the optimization candidates AND the fixed initial values
                #print 'energies in computeGeneration', trialEnergies
                if self.dif == '+':
                    difs = []
                    for i in range(self.nmodels):
                        for j in range(i+1,self.nmodels):
                            difs.append(trialEnergies[i] - trialEnergies[j])
                    self.generationEnergyList.append(difs)
                elif self.dif == '-':
                    difs = []
                    for i in range(self.nmodels):
                        for j in range(i+1,self.nmodels):
                            difs.append(trialEnergies[j] - trialEnergies[i])
                    self.generationEnergyList.append(difs)
                else:
                    self.generationEnergyList.append(trialEnergies)
            timeElapsed = time() - time0
            print 'generation took', timeElapsed, 's'
            self.ftimes.append(timeElapsed)
        
        else: # generation >= 1
            time0 = time()
            print '\n\ngeneration', self.generation, '\n'
            self.oldGeneration = self.population
            global objectiveDic
            objectiveDic = {}
            global dominanceDic
            dominanceDic = {}
            solutionDic = {}
            
            self.newPopulation = []
            
            newPopulationList = []
            oldPopulationList = [list(i) for i in self.population]
            
            for candidate in range(self.populationSize):
                # generate new solutions.,reject those out-of-bounds or repeated
                insideBoundaries = False
                insideBoundariesCounter = 0
                while insideBoundaries == False:
                    self.calcTrialSolution(candidate)
                    ltrial = list(self.trialSolution)
                    if ltrial not in newPopulationList and ltrial not in oldPopulationList:
                        newPopulationList.append(ltrial)
                    else:
                        continue
                    for z in range(len(self.trialSolution)):
                        if self.trialSolution[z] <= self.toOptBounding[1][z] and self.trialSolution[z] >= self.toOptBounding[0][z]:
                            insideBoundariesCounter += 1
                            if insideBoundariesCounter == len(self.trialSolution):
                                insideBoundaries = True
                        else:
                            insideBoundariesCounter = 0
                            break
                
                self.newPopulation.append(self.trialSolution)
                
                trialEnergies, self.atSolution = self.EnergyFunction(self.newPopulation[candidate])

                if self.dif == '+':
                    difs = []
                    for i in range(self.nmodels):
                        for j in range(i+1,self.nmodels):
                            difs.append(trialEnergies[i] - trialEnergies[j])
                    self.newGenerationEnergyList[candidate] = difs
                elif self.dif == '-':
                    difs = []
                    for i in range(self.nmodels):
                        for j in range(i+1,self.nmodels):
                            difs.append(trialEnergies[j] - trialEnergies[i])
                    self.newGenerationEnergyList[candidate] = difs
                else:
                    self.newGenerationEnergyList[candidate] = trialEnergies
              
##             print self.generation,': len new energy list', len(self.newGenerationEnergyList)
            
            print "finished generating new candidates"
            energyComparison = dominanceComparison(self.newGenerationEnergyList, self.generationEnergyList)
            print 'dominance comparison with previous generation:'
            print energyComparison
            
            dicIndex = 1
            newBetterSols = 0

            for i in range(len(energyComparison)):
                if energyComparison[i] > 0:
                    objectiveDic[dicIndex] = numpy.copy(self.newGenerationEnergyList[i])
                    dominanceDic[dicIndex] = []
                    solutionDic[dicIndex] = numpy.copy(self.newPopulation[i])
                    dicIndex = dicIndex +1
                    newBetterSols += 1
                if energyComparison[i] < 0:
                    objectiveDic[dicIndex] = numpy.copy(self.generationEnergyList[i])
                    dominanceDic[dicIndex] = []
                    solutionDic[dicIndex] = numpy.copy(self.population[i])
                    dicIndex = dicIndex +1
                if energyComparison[i] == 0:
                    objectiveDic[dicIndex] = numpy.copy(self.generationEnergyList[i])
                    dominanceDic[dicIndex] = []
                    solutionDic[dicIndex] = numpy.copy(self.population[i])
                    dicIndex = dicIndex +1
                    objectiveDic[dicIndex] = numpy.copy(self.newGenerationEnergyList[i])
                    dominanceDic[dicIndex] = []
                    solutionDic[dicIndex] = numpy.copy(self.newPopulation[i])
                    dicIndex = dicIndex +1
            
            self.getDominanceTree(objectiveDic.keys())
            nonDominatedFrontsOut = self.orgndf(dominanceDic)
            
            #print 'objectiveDic a', objectiveDic
            
            self.population = []
            self.generationEnergyList = []
            
            self.fronts = []    # holds iterations of fronts in this generation
            self.frontObj = []  # holds objective values for each front
            
            full = False
            while not full:
                for i in range(len(nonDominatedFrontsOut)):
                    if len(self.population) + len(nonDominatedFrontsOut[i]) > self.populationSize and len(self.population) < self.populationSize:
                        tempFront = nonDominatedFrontsOut[i]
                        tempObjDic = {}
                        for k in tempFront:
                            tempObjDic[k] = objectiveDic[k]
                        while len(self.population) + len(tempObjDic.keys()) != self.populationSize:
                            tempObjDic = removeMostCrowded(tempObjDic, 3)
                        full = True
                        break
                    elif len(self.population) == self.populationSize:
                        full = True
                        break
                    else:
                        self.fronts.append([])
                        self.frontObj.append([])
                        for k in nonDominatedFrontsOut[i]:
                            self.population.append(solutionDic[k])
                            self.generationEnergyList.append(numpy.copy(objectiveDic[k]))
                            self.fronts[-1].append(numpy.copy(solutionDic[k]))
                            self.frontObj[-1].append(numpy.copy(objectiveDic[k]))
                        tempObjDic = {}
            
##             print 'fronts0:',[len(i) for i in self.fronts]
##             print 'pop0:',len(self.population)
            if self.fronts ==[]:
                self.fronts.append([])
            if self.frontObj ==[]:
                self.frontObj.append([])
            
            # complete pop to self.populationSize and complete last front
            for i in tempObjDic.keys():
                self.population.append(numpy.copy(solutionDic[i]))
                self.generationEnergyList.append(numpy.copy(objectiveDic[i]))
                self.fronts[-1].append(numpy.copy(solutionDic[i]))
                self.frontObj[-1].append(numpy.copy(objectiveDic[i]))
            
##             print 'fronts:',[len(i) for i in self.fronts]
##             print 'pop:',len(self.population)

            print "\ngeneration %d finished"%(self.generation)
            print 'number of new better solutions', newBetterSols
            print 'room for improvement', self.roomForImprovement
            
            if ((not self.trueMetric) and newBetterSols <= self.roomForImprovement and len(self.fronts) == 1) or (self.trueMetric and self.nmodels == 2 and newBetterSols <= self.roomForImprovement):
                self.generationsWithNoImprovement += 1
            else:
                self.generationsWithNoImprovement = 0
            timeElapsed = time() - time0
            print 'generation took', timeElapsed, 's'
            self.ftimes.append(timeElapsed)
            
            if self.keep_track:
                for front in self.fronts:
                    print >> self.frontsfile, [list(elem) for elem in front]
                for objs in self.frontObj:
                    print >> self.objsfile, [list(elem) for elem in objs]
        
        self.generation += 1
        print 'generations with no improvement:', self.generationsWithNoImprovement
        return

    def finalize(self):
        if self.keep_track:
            self.frontsfile.close()
            self.objsfile.close()

    #------------------------------------------------------------------------------------------------------------------------------------
    #This code is an adaptation of the non-dominated sorting algorithm with delayed insertion published in
    #Fang et al (2008) An Efficient Non-dominated Sorting Method for Evolutionary Algorithms, Evolutionary Computation 16(3):355-384
    #This algorithm is so far implemented without delayed insertion. 
    #Related to that or not, in some situations comparisons are repeated - for example use test3/seed(2)/2 objectives and 10 solutions.
    #Modified for the last time on April 14th, 2009

    def getDominanceTree(self, nodeList):
        """ This function applies the 'Divide and Conquer' method recursively generating
        the dominance tree (divides) and returning the product of the mergeDominanceTree function (conquers).
        This version doesn't use delayed insertion of dominated node yet. """
        size = len(nodeList)
        if size > 1:
            leftTree = self.getDominanceTree(nodeList[:size/2])
            rightTree = self.getDominanceTree(nodeList[size/2:])
        else:
            return nodeList
        return self.mergeDominanceTrees(leftTree, rightTree)

    def mergeDominanceTrees(self, leftTree, rightTree):
        """ This function merges (conquers) the dominance trees recursively (not using delayed insertion yet). """
        leftNode = Node(leftTree[0])
        rightNode = Node(rightTree[0])
        switch = 0 #switch indicates at the end of the loop which conditional of the loop is used. If 'else' is used switch does not change.
        while leftNode.nodeIndex != -1 and rightNode.nodeIndex != -1: #and leftTree != [] and rightTree != []:
            switch = 0 #switch indicates at the end of the loop which conditional of the loop is used. If 'else' is used switch does not change.
            dominanceTestResult = self.dominanceTest(leftNode.nodeIndex, rightNode.nodeIndex)
            #rightDelayedInsertionList and leftDelayedInsertionList are [] for now just for the sake of simplicity.
            leftDelayedInsertionList = []
            rightDelayedInsertionList = []
            if dominanceTestResult < 0: 
                switch = 1
                if rightDelayedInsertionList != []:
                    mergeDominanceTrees([dominanceDic[rightNode.nodeIndex][0]], rightDelayedInsertionList)
                tempNode = rightNode.nodeIndex
                rightNode.nextSiblingNode(rightTree) 
                for i in dominanceDic.keys():
                    if tempNode in dominanceDic[i]:
                        dominanceDic[i].remove(tempNode)
                if len(dominanceDic[leftNode.nodeIndex]) > 0:
                    dominanceDic[leftNode.nodeIndex] = self.mergeDominanceTrees(dominanceDic[leftNode.nodeIndex], [tempNode])
                else:
                    dominanceDic[leftNode.nodeIndex].append(tempNode)
                if tempNode in rightTree:
                    rightTree.remove(tempNode)
                if leftNode.nodeIndex != -1 and rightNode.nodeIndex == -1 and rightTree != []: #New!
                    rightNode.nodeIndex = rightTree[0]
            elif dominanceTestResult > 0:
                switch = 2
                if leftDelayedInsertionList != []:
                    self.mergeDominanceTrees([dominanceDic[leftNode.nodeIndex][0]], leftDelayedInsertionList)
                tempNode = leftNode.nodeIndex
                leftNode.nextSiblingNode(leftTree) 
                for i in dominanceDic.keys():
                    if tempNode in dominanceDic[i]:
                        dominanceDic[i].remove(tempNode)
                if len(dominanceDic[rightNode.nodeIndex]) > 0:
                    dominanceDic[rightNode.nodeIndex] = self.mergeDominanceTrees(dominanceDic[rightNode.nodeIndex], [tempNode])
                else:
                    dominanceDic[rightNode.nodeIndex].append(tempNode)
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
                    dominanceDic[rightNode.nodeIndex] = self.mergeDominanceTrees(dominanceDic[rightNode.nodeIndex], rightDelayedInsertionList)
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
        if dominanceDic != {}:
            firstFrontToBeCompared = [0, 0]
            grandChildren = []
            #This finds and merges cousins which are children of dominated solutions
            for i in leftTree:
                if dominanceDic[i] != []:
                    grandChildren.append(dominanceDic[i])
                    if len(grandChildren) == 2:
                        grandChildren = [self.mergeDominanceTrees(grandChildren[0], grandChildren[1])]
            #Continues to merge the cousins which are children of non-dominated parents
            while len(firstFrontToBeCompared) > 1:
                firstFrontToBeCompared = []
                for i in leftTree:
                    nonDominantSolution = True
                    for j in dominanceDic.keys():
                        if i in dominanceDic[j]:
                            nonDominantSolution = False
                    if nonDominantSolution == True and dominanceDic[i]!=[]:
                        firstFrontToBeCompared.append(i)
                    if len(firstFrontToBeCompared) > 1:
                        #Next two lines are necessary to avoid errors inside merge function
                        cousin1 = dominanceDic[firstFrontToBeCompared[0]] #The 'list' function is necessary for the merging to work
                        cousin2 = dominanceDic[firstFrontToBeCompared[1]] #If 'list' is not used an instance of a dictionary object passes into the merging function and that raises problems
                        self.mergeDominanceTrees(cousin1, cousin2)
                        break
        return leftTree

    def orgndf(self, dominanceDic):
        """ This function organizes the nondominated fronts found by the function getDominanceTree in a list. """
        nonDominatedFronts = []
        while dominanceDic != {}:
            nonDominated = []
            for i in dominanceDic.keys():
                dominated = False
                for j in dominanceDic.values():
                    if i in j:
                        dominated = True
                        break
                if dominated == False:
                    nonDominated.append(i)
            nonDominatedFronts.append(nonDominated)
            for k in nonDominated:
                del dominanceDic[k]
        return nonDominatedFronts

    def dominanceTest(self, solution1, solution2):
        """ This function tests the dominance relationship between two solutions,
        given as keys in, objectiveDic, the dictionary containing as values the lists with the objective functions evaluations. """
        dominanceTestResult = 0
        for i,j in zip(objectiveDic[solution1],objectiveDic[solution2]):
            d = i - j
            if d <= 0 and dominanceTestResult <=0:
                dominanceTestResult += -1
            elif d >= 0 and dominanceTestResult >=0:
                dominanceTestResult += 1
            else:
                dominanceTestResult = 0
                break
        return dominanceTestResult
    #------------------------------------------------------------------------------------------------------------------------------------
    #End of non-dominated sorted algorithm methods

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
    
    if len(x) < 3:
        x.popitem()
        return x
    else:
        keys = x.keys()
        keys.sort()
        extremes = []
        for i in range(len(random.choice(x.values()))):
            maximum = x[keys[0]][i]
            tempMax = [keys[0]]
            minimum = x[keys[0]][i]
            tempMin = [keys[0]]
            for j in keys:
                if x[j][i] < minimum:
                    tempMin = []
                    minimum = x[j][i]
                    tempMin.append(j)
                elif x[j][i] == minimum and j not in tempMin:
                    tempMin.append(j)                
                elif x[j][i] > maximum:
                    tempMax = []
                    maximum = x[j][i]
                    tempMax.append(j)
                elif x[j][i] == maximum and j not in tempMax:
                    tempMax.append(j)                
            for k in tempMin:
                if k not in extremes:
                    extremes.append(k)
            for k in tempMax:
                if k not in extremes:
                    extremes.append(k)
        for i in keys:
            x[i] = array(x[i])
        lista =[]
        #Use a fast implementation of the Euclidean distance
        temp = numpy.zeros(len(random.choice(x.values())))
        # Predefining temp allows reuse of this array, making this
        # function about twice as fast.
        distanceMatrix = []
        for i in range(len(keys)):
            distances = []  # list of (distance, index)
            for j in range(len(keys)):
                if j < i:
                    distances.append(distanceMatrix[j][i])
                elif j ==i:
                    distances.append(0)
                elif j > i:
                    temp[:] = x[keys[i]] - x[keys[j]]
                    dist = numpy.sqrt(numpy.dot(temp,temp))
                    distances.append(dist)
            distanceMatrix.append(numpy.copy(distances))
            distances.sort()
            if knumber < len(distances):
                lista.append(distances[1:knumber+1])
            if knumber >= len(distances):
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
        for i in x.keys(): #Confirmar se este objecto deve ser retornado com listas ou arrays
            x[i] = list(x[i]) #de modo a ser usado pelas outras funções na geração seguinte.
        return x


#________________________________________________________________
#Tests for the non-dominated sorted algorithm methods

if __name__ == "__main__":

    from numpy import random
    from numpy.random import rand, seed

    class ndsaTest(GDE3Solver):

        def __init__(self, numberOfNodes, numberOfObjectives, printOption = False):
            self.numberOfNodes = numberOfNodes
            self.numberOfObjectives = numberOfObjectives
            print 'Test initialized'
            seed(2)
            #Attention! Because these variables are defined as global, they must be reset every time a test starts or each generation during an optimization with an EA. 
            global dominanceDic
            global objectiveDic #This must begin in 1.
            dominanceDic = {}
            objectiveDic = {}
            print 'self.numberOfNodes', numberOfNodes
            print 'self.numberOfObjectives', self.numberOfObjectives
            for nodeNumber in range(self.numberOfNodes):
                dominanceDic[nodeNumber+1] = []
                objectiveDic[nodeNumber+1] = []
                for objectiveNumber in range(self.numberOfObjectives):
                    objectiveDic[nodeNumber+1].append(rand())
            print 'dominanceDic and objectiveDic created; entering getDominanceTree'
            print 'objectiveDic', objectiveDic
            chaves = objectiveDic.keys()
            self.getDominanceTree(chaves)
            print self.numberOfNodes, 'nodes and ', self.numberOfObjectives, 'objectives'
            print 'Entering nonDominatedFrontsOut'
            nonDominatedFrontsOut = self.orgndf(dominanceDic)
            print 'objectiveDic', objectiveDic
            print 'nonDominatedFrontsOut', nonDominatedFrontsOut
            print 'Testing non-dominance between solutions in the same front...'
            for k in nonDominatedFrontsOut:
                if len(k) == 1:
                    dominance =0
                elif len(k) > 1:
                    p = 0
                    while p < len(k)-1:
                        r = p + 1
                        while r < len(k):
                            dominance = self.dominanceTest(k[p], k[r])
                            if dominance != 0:
                                print '\n\n\nNumber of solutions', self.numberOfNodes, 'number of objectives', self.numberOfObjectives
                                print 'Domination relationship in front', k, 'between nodes', k[p], 'and', k[r],'. Test not passed.\n\n'
                                break
                            else:
                                r += 1
                        if dominance != 0:
                            break
                        else:
                            p += 1
                    if dominance != 0:
                        break
                if dominance != 0:
                    return
                elif len(k) == 0:
                    print 'Empty front found!'
                    return
            print 'Non-dominance test between solutions in the same front passed.'
            if len(nonDominatedFrontsOut) > 1:
                print 'Testing dominance relationship between solutions in different fronts...'
                #Solution in rFront must be dominated by at least one solution in pFront and cannot dominate any solution in pFront.
                pFront = 0
                rFront = pFront + 1
                while pFront < len(nonDominatedFrontsOut)-1:
                    for down in nonDominatedFrontsOut[rFront]:
                        print 'Assigning totalDominance for the first time in the test'
                        totalDominance = 0
                        for up in nonDominatedFrontsOut[pFront]:
                            dominance = self.dominanceTest(up, down)
                            if dominance == 1:
                                print '\n\n\nNumber of solutions', self.numberOfNodes, 'number of objectives', self.numberOfObjectives
                                print 'Solution in front', pFront, ', (', up, ') is dominated by solution in front', rFront, ', (', down, '). Test not passed.\n\n'
                                break
                            totalDominance = totalDominance + dominance
                        if dominance == 1 or totalDominance < 1:
                            break
                    if dominance == 1 or totalDominance < 1:
                        break
                    pFront += 1
                    rFront += 1
                if dominance == 1 or totalDominance < 1:
                    return #pFront, up, rFront, down
                print 'Non-dominance test between solutions in different fronts passed.'
            if len(nonDominatedFrontsOut) == 1:
                print 'Only non-dominated solutions - no test between different fronts'
            print 'Tests finished successfully!\n\n\n'

    counter = 0
    for i in range(10, 60):
        for j in range(2, 5):
            if (i == 0 and j == 0):
                print '\n\n\n\n\n'
                ndsaTest(i, j, printOption = False)                    #(numberOfNodes, numberOfObjectives)
                print '\n\n\n\n\n'
            else:
                ndsaTest(i, j)                    #(numberOfNodes, numberOfObjectives)
            counter += 1
            print 'test', counter
    print 'Tests finished successfully!'
    
