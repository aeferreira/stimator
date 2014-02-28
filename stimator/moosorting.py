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


class MOOSorter(object):

    def __init__(self):
        self.dom_dict = {}
        self.dominance_rdepth = 0
        
    def depthspaces(self):
        return ' '*(self.dominance_rdepth-1)*4

    #------------------------------------------------------------------------------------------------------------------------------------
    #This code is an adaptation of the non-dominated sorting algorithm with delayed insertion published in
    #Fang et al (2008) An Efficient Non-dominated Sorting Method for Evolutionary Algorithms, Evolutionary Computation 16(3):355-384

    
    def getDominanceTree(self, nodeList, verbose = False):
        """ This function applies the 'Divide and Conquer' method recursively generating
        the dominance tree (divides) and returning the product of the mergeDominanceTree function (conquers).
        This version doesn't use delayed insertion of dominated node yet."""
        self.dominance_rdepth += 1
        if verbose:
            spcs = self.depthspaces()
            print spcs, '---------------------------------'
            print spcs, 'IN getDominanceTree()'
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
            self.dominance_rdepth -= 1
            return res
        else:
            if verbose:
                print spcs, '.single node return', nodeList
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
            if d < 0: 
                one_dominates = True
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
##         #This block merges cousins, i.e. the nodes at the same non-dominance level which are children of different non-dominant siblings
##         #I'm not sure if this should be done every time or just in the final.
##         if self.dom_dict != {}:
##             firstFrontToBeCompared = [0, 0]
##             grandChildren = []
##             #This finds and merges cousins which are children of dominated solutions
##             for i in left_tree:
##                 if self.dom_dict[i] != []:
##                     grandChildren.append(self.dom_dict[i])
##                     if len(grandChildren) == 2:
##                         grandChildren = [self.merge_Dom_trees(grandChildren[0], grandChildren[1])]
##             #Continues to merge the cousins which are children of non-dominated parents
##             while len(firstFrontToBeCompared) > 1:
##                 firstFrontToBeCompared = []
##                 for i in left_tree:
##                     nonDominantSolution = True
##                     for k in self.dom_dict:
##                         if i in self.dom_dict[k]:
##                             nonDominantSolution = False
##                     if nonDominantSolution and self.dom_dict[i]!=[]:
##                         firstFrontToBeCompared.append(i)
##                     if len(firstFrontToBeCompared) > 1:
##                         #Next two lines are necessary to avoid errors inside merge function
##                         cousin1 = self.dom_dict[firstFrontToBeCompared[0]]
##                         cousin2 = self.dom_dict[firstFrontToBeCompared[1]]
##                         self.merge_Dom_trees(cousin1, cousin2)
##                         break
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
        

#________________________________________________________________
#Tests for the non-dominated sorted algorithm methods

from scipy.spatial.distance import squareform, pdist
if __name__ == "__main__":

    class FangEtAL_test(MOOSorter):
        def __init__(self, report = False):
            MOOSorter.__init__(self)
            
            data = [[182.08, 100.13, 192.21],[187.53, 246.16, 203.2],
                    [197.15, 201.57, 318.86],[47.48, 74.96, 22.69],
                    [37.05, 304.83, 381.19], [126.88, 54.58, 144.17],
                    [101.77, 49.18, 111.91], [37.47, 18.63, 446.57]]
            
            print 'distance matrix'
            print squareform(pdist(data, dominance))
            print '##########################################'
            print

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
            print '======================================================'
            print 'Comparing to published results...'
            published =[[4, 5, 7, 8],[6],[1],[2, 3]]
            if nondominated_waves == published:
                print 'PASSED'
            else:
                print 'FAILED'
            print '======================================================'

    FangEtAL_test()
    
