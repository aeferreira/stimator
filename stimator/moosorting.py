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

    def __init__(self, obj_dict = {}):
        self.objectives = obj_dict
        self.keys = list(obj_dict.keys())
        self.dom_dict = {}
        for k in self.keys:
            self.dom_dict[k] = []
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
            print spcs, 'getDominanceTree(%s)'%str(nodeList)
        size = len(nodeList)
        if size > 1:
            leftTree  = self.getDominanceTree(nodeList[:size/2], verbose)
            rightTree = self.getDominanceTree(nodeList[ size/2:], verbose)
            if verbose:
                print spcs, '@@  left tree', leftTree
                print spcs, '@@ right tree', rightTree
            res = self.merge_Dom_trees(leftTree, rightTree)
            if verbose:
                print spcs, '--> node list returned from merge', res
                print spcs, '--> dom dict', self.dom_dict
            self.dominance_rdepth -= 1
            return res
        else:
            if verbose:
                print spcs, '.single node return', nodeList
            self.dominance_rdepth -= 1
            return nodeList

    def merge_Dom_trees(self, left_tree, right_tree):
        """ This function merges (conquers) the dominance trees recursively 
        (not using delayed insertion yet). """
        if len(left_tree) == 0:
            if len(right_tree) < 2:
                return right_tree
        left = Node(left_tree[0])
        right = Node(right_tree[0])
        while left.valid() and right.valid():
            d = dominance(self.objectives[right.indx], self.objectives[left.indx])
            if d < 0: 
                node = right.indx
                right.move_and_remove(right_tree) 
                for k in self.dom_dict:
                    if node in self.dom_dict[k]:
                        self.dom_dict[k].remove(node)
                self.dom_dict[left.indx] = self.merge_Dom_trees(self.dom_dict[left.indx], [node])
                if not right.valid() and len(right_tree) > 0:
                    right.indx = right_tree[0]
                    left.point_to_next_sibling(left_tree)
            elif d > 0:
                node = left.indx
                left.move_and_remove(left_tree) 
                for k in self.dom_dict:
                    if node in self.dom_dict[k]:
                        self.dom_dict[k].remove(node)
                self.dom_dict[right.indx] = self.merge_Dom_trees(self.dom_dict[right.indx], [node])
                right.indx = right_tree[0]
            else:
                right.point_to_next_sibling(right_tree)
                if not right.valid():
                    right.indx = right_tree[0]
                    left.point_to_next_sibling(left_tree)
        left_tree.extend(right_tree)
        return left_tree

    def get_non_dominated_fronts(self, nodeList, verbose = False):
        fronts = []
        top_nodes = nodeList
        while len(top_nodes) > 0:
            top_nodes = self.getDominanceTree(top_nodes, verbose)
            fronts.append(top_nodes)
            if verbose:
                print '\nnon-dominated FRONT :', top_nodes
            children = []
            for node in top_nodes:
                children.extend(self.dom_dict[node])
            if verbose:
                print 'CHILDREN :', children
            top_nodes = children
        return fronts
        
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

def pprint_dominance_matrix(keys, values, f):
    keys = ['%3s'%str(k) for k in keys]
    header = '    ' + ' '.join(keys)
    print header
    for i, k in enumerate(keys):
        line = ''
        for j in range(len(values)):
            if i == j:
                c = '0'
            else:
                d = f(values[i], values[j])
                if d < 0:
                    c = '<'
                elif d > 0:
                    c = '>'
                else:
                    c = '0'
            line = line + '%2s' % c
        line = ' '.join(line)
        print '%2s' % k, line
        
if __name__ == "__main__":

    class FangEtAL_test(MOOSorter):
        def __init__(self):
            
            data = [[182.08, 100.13, 192.21],[187.53, 246.16, 203.2],
                    [197.15, 201.57, 318.86],[47.48, 74.96, 22.69],
                    [37.05, 304.83, 381.19], [126.88, 54.58, 144.17],
                    [101.77, 49.18, 111.91], [37.47, 18.63, 446.57]]

            n_nodes = len(data)
            n_objectives = len(data[0])
            print '======================================================'
            print 'Test initialized: example from'
            print 'Fang et al (2008) An Efficient Non-dominated Sorting Method'
            print 'for Evolutionary Algorithms, Evol. Comput. 16(3):355-384'
            print '\nNodes: %d  Objectives: %d' % (n_nodes, n_objectives)
            #dict keys are integers that begin at 1.
            objectives = {}
            for i_node in range(n_nodes):
                objectives[i_node+1] = data[i_node]
            
            MOOSorter.__init__(self, objectives)
            
            print 'dominance matrix'
            pprint_dominance_matrix(self.keys, self.objectives.values(), dominance)
            print '------------------------------------------------'
            print
            print '\nComputing nondominated fronts...'
            
            fronts = self.get_non_dominated_fronts(self.keys, verbose = True)
            
            print
            print 'final DOMINANCE dict'
            for k in self.dom_dict:
                print k, '>', self.dom_dict[k]
            print '\n%d non-dominated fronts:'% len(fronts)
            for front in fronts:
                print front
            print '======================================================'
            print 'Comparing to published results...',
            published =[[4, 5, 7, 8],[6],[1],[2, 3]]
            if fronts == published:
                print 'PASSED'
            else:
                print 'FAILED'
            print
    
    class ndsaTest(MOOSorter):

        def __init__(self, n_nodes, n_objectives, report=False):

            print '======================================================'
            print 'Random objectives test initialized'
            print 'Nodes: %d  Objectives: %d' % (n_nodes, n_objectives)
            numpy.random.seed(2)
            #dict keys are integers that begin at 1.
            objectives = {}
            for i_node in range(n_nodes):
                objectives[i_node+1] = []
                for i_objective in range(n_objectives):
                    objectives[i_node+1].append(numpy.random.rand())
            
            MOOSorter.__init__(self, objectives)

            print 'objectives dict: %d keys with %d elements' % (len(self.objectives), len(self.objectives[1]))
            print 'dominance matrix'
            pprint_dominance_matrix(self.keys, self.objectives.values(),dominance)
            print '----------------------------------------------------'
            
            fronts = self.get_non_dominated_fronts(self.keys, verbose = report)
            
            print
            print 'final DOMINANCE dict'
            for k in self.dom_dict:
                print k, '>', self.dom_dict[k]
            print '\n%d non-dominated fronts:'% len(fronts)
            for front in fronts:
                print front
            print '======================================================'
            print '\nTesting non-dominance between solutions in the same front...',
            for front in fronts:
                if len(front) == 0:
                    print '\nFAILED: empty front found!'
                    return
                if len(front) == 1:
                    continue
                for p1 in range(len(front)-1):
                    for p2 in range(p1+1, len(front)):
                        d = dominance(self.objectives[front[p2]], self.objectives[front[p1]])
                        if d != 0:
                            print '\nFAILED:'
                            w1 = front[p1]
                            w2 = front[p2]
                            print 'Domination relationship in front', front, 'between nodes', w1, 'and', w2,'.\n\n'
                            return
            print 'passed.'
            if len(fronts) == 1:
                print 'Only non-dominated solutions - no test between different fronts'
                return
            print 'Testing dominance relationship between different fronts...',
            #Solution in rFront must be dominated by at least one solution in pFront and cannot dominate any solution in pFront.
            for pFront in range(len(fronts)-1):
                rFront = pFront + 1
                one_dominates = False
                for down in fronts[rFront]:
                    for up in fronts[pFront]:
                        d = dominance(self.objectives[down], self.objectives[up])
                        if d > 0:
                            print '\nFAILED:'
                            print 'Solution', up, ' is dominated by solution ', down, '.\n\n'
                            return
                        if d < 0:
                            one_dominates = True
                if not one_dominates:
                    print '\nFAILED:'
                    print 'No solution in front', fronts[pFront], ' dominates any solution in front ', fronts[rFront], '.\n\n'
                    return                            
            print 'passed.'

    FangEtAL_test()
    ndsaTest(15, 2, report=True)
    ndsaTest(50, 2)
    
