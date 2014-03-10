
import numpy

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


class MOOSorter(object):

    def __init__(self, obj_dict = {}):
        self.objectives = obj_dict
        self.keys = list(obj_dict.keys())
        self.dom_dict = {}
        for k in self.keys:
            self.dom_dict[k] = []
        self.dominance_rdepth = 0
    
    def remove_refs2node(self, node):
        for k in self.dom_dict:
            if node in self.dom_dict[k]:
                self.dom_dict[k].remove(node)

    #------------------------------------------------------------------------------------------------------------------------------------
    #This code is based on the non-dominated sorting algorithm with delayed insertion by
    #Fang et al (2008) An Efficient Non-dominated Sorting Method for Evolutionary Algorithms,
    #Evolutionary Computation 16(3):355-384

    
    def getDominanceTree(self, nodeList, verbose = False):
        """Applies the 'Divide and Conquer' method to recursively generate
        the dominance tree (divides) and returns the result 
        of mergeDominanceTree() (conquers).
        """

        if verbose:
            self.enter_getDominanceTree(nodeList)

        size = len(nodeList)
        if size > 1:
            leftTree  = self.getDominanceTree(nodeList[:size/2], verbose)
            rightTree = self.getDominanceTree(nodeList[ size/2:], verbose)
            if verbose:
                self.print_lr_trees(leftTree, rightTree)
            res = self.merge_Dom_trees(leftTree, rightTree)
        else:
            res = nodeList
        if verbose:
            self.exit_getDominanceTree(res)
        return res

    class NodeListIterator(object):
        """NodeListIterator supports iteration of a list (of necessarily different elements)
        with possible removal of elements during the iteration process.
        current() is the basic method to retrieve the current element.
        reset() is used to point current() to first element of the list and
        remove() deletes current(), leaving it pointing to the next element in the list.
        valid() tests the past-the-end condition.
        """

        def __init__(self, nodelist):
            """Initializes Node object with an index"""
            self.node_list = nodelist
            self.reset()
        
        def reset(self):
            if len(self.node_list) < 1:
                self.__valid = False
            else:
                self.node = self.node_list[0]
                self.__valid = True

        def current(self):
            if not self.__valid:
                raise IndexError
            return self.node

        def next(self):
            """Increments the index to the next sibling of the node, if not past the end of the 'tree' list"""
            if not self.__valid:
                raise IndexError
            pos = self.node_list.index(self.node)
            if pos < len(self.node_list)-1:
                self.node = self.node_list[pos+1]
            else:
                self.__valid = False
        
        def valid(self):
            return self.__valid
            
        def remove(self):
            if not self.__valid:
                raise IndexError
            node = self.node
            self.next()
            self.node_list.remove(node)
        
    def merge_Dom_trees(self, left_tree, right_tree):
        """Merges (conquers) two dominance trees recursively 
        (not using delayed insertion yet)."""
        
        if len(left_tree) == 0:
            if len(right_tree) < 2:
                return right_tree
        left = self.NodeListIterator(left_tree)
        right = self.NodeListIterator(right_tree)
        while left.valid() and right.valid():
            left_value = left.current()
            right_value = right.current()
            d = dominance(self.objectives[right_value], self.objectives[left_value])
            if d < 0: 
                right.remove()
                self.remove_refs2node(right_value)
                self.dom_dict[left_value] = self.merge_Dom_trees(self.dom_dict[left_value], [right_value])
                if not right.valid() and len(right_tree) > 0:
                    right.reset()
                    left.next()
            elif d > 0:
                left.remove()
                self.remove_refs2node(left_value)
                self.dom_dict[right_value] = self.merge_Dom_trees(self.dom_dict[right_value], [left_value])
                right.reset()
            else:
                right.next()
                if not right.valid():
                    right.reset()
                    left.next()
        left_tree.extend(right_tree)
        return left_tree

    def get_non_dominated_fronts(self, verbose = False):
        """Driver function to obtain a set of non-dominated fronts from sorting algorithm.
        The function applies sucessive rounds of sorting, applying a divide and conquer approach.
        After the generation of the first dominance tree, the first front is extracted as the top elements.
        In the next rounds, the same process is applied to the children of the top elements.
        When no more children are found, the function exits.
        """
        top_nodes = self.keys # start with the whole node list
        fronts = []
        while len(top_nodes) > 0:
            top_nodes = self.getDominanceTree(top_nodes, verbose)
            fronts.append(top_nodes)
            children = []
            for n in top_nodes:
                children.extend(self.dom_dict[n])
            if verbose:
                print '\nNon-dominated FRONT :', top_nodes
                print 'CHILDREN :', children
            top_nodes = children
        return fronts
        
    #------------------------------------------------------------------------------------------------------------------------------------
    #Pretty print functions used in testing

    def depthspaces(self):
        return ' '*(self.dominance_rdepth-1)*4
    
    def enter_getDominanceTree(self, nodeList):
        self.dominance_rdepth += 1
        spcs = self.depthspaces()
        print spcs, '---------------------------------'
        if len(nodeList) == 1:
            print spcs, 'getDominanceTree(%s) --> %s'% (str(nodeList), str(nodeList))
        else:
            print spcs, 'getDominanceTree(%s)'% str(nodeList)
            

    def exit_getDominanceTree(self, res):
        spcs = self.depthspaces()
        if len(res) > 1:
            print spcs, '--> node list returned from merge', res
            print spcs, '--> dom dict', self.dom_dict
        self.dominance_rdepth -= 1
    
    def print_lr_trees(self, leftTree, rightTree):
        spcs = self.depthspaces()
        print spcs, '@@  left tree', leftTree
        print spcs, '@@ right tree', rightTree

    def pprint_dominance_matrix(self, f):
        keys = ['%3s'%str(k) for k in self.keys]
        values = self.objectives.values()
        print '    ' + ' '.join(keys)
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

        
#---------------------------------------------------------------------------
# Tests for the non-dominated sorting 
#---------------------------------------------------------------------------

def FangEtAL_test():
        
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

    objectives = {}
    for i in range(n_nodes):
        objectives[i+1] = data[i]
    
    sorter = MOOSorter(objectives)
    
    print 'dominance matrix'
    sorter.pprint_dominance_matrix(dominance)
    print '------------------------------------------------'
    print
    print '\nComputing nondominated fronts...'
    
    fronts = sorter.get_non_dominated_fronts(verbose = True)
    
    print
    print 'final DOMINANCE dict'
    for k in sorter.dom_dict:
        print k, '>', sorter.dom_dict[k]
    print '\n%d non-dominated fronts:'% len(fronts)
    for front in fronts:
        print front
    print '======================================================'
    print 'Comparing with published results...',
    published =[[4, 5, 7, 8],[6],[1],[2, 3]]
    if fronts == published:
        print 'PASSED'
    else:
        print 'FAILED'
    print


def random_objs_test(n_nodes, n_objectives, report=False):

    print '======================================================'
    print 'Random objectives test initialized'
    print 'Nodes: %d  Objectives: %d' % (n_nodes, n_objectives)

    numpy.random.seed(2)
    objectives = {}
    for i_node in range(n_nodes):
        objectives[i_node+1] = []
        for i_objective in range(n_objectives):
            objectives[i_node+1].append(numpy.random.rand())
    
    sorter = MOOSorter(objectives)

    print 'dominance matrix'
    sorter.pprint_dominance_matrix(dominance)
    print '----------------------------------------------------'
    
    fronts = sorter.get_non_dominated_fronts(verbose = report)
    
    print
    print 'final DOMINANCE dict'
    for k in sorter.dom_dict:
        print k, '>', sorter.dom_dict[k]
    print '\n%d non-dominated fronts:'% len(fronts)
    for front in fronts:
        print front
    print '======================================================'
    print '\nTesting non-dominance between solutions in the same front...',
    #Every solution in each front can not dominate any other solution in the same front.
    for front in fronts:
        if len(front) == 0:
            print '\nFAILED: empty front found!'
            return
        if len(front) == 1:
            continue
        for p1 in range(len(front)-1):
            for p2 in range(p1+1, len(front)):
                d = dominance(objectives[front[p2]], objectives[front[p1]])
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
    
    #Every solution in rFront must be dominated by at least one solution in pFront 
    #and 
    #can not dominate any solution in pFront.
    
    for pFront in range(len(fronts)-1):
        rFront = pFront + 1
        one_dominates = False
        for down in fronts[rFront]:
            for up in fronts[pFront]:
                d = dominance(objectives[down], objectives[up])
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


if __name__ == "__main__":

    FangEtAL_test()
    random_objs_test(15, 2, report=True)
    random_objs_test(50, 2)
    
