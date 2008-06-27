
import sys
sys.path.append('..')
import DESolver, numpy

class TestSolver(DESolver.DESolver):

    #some test data
    xData = numpy.array([5.357, 9.861, 5.457, 5.936, 6.161, 6.731])
    yData = numpy.array([0.376, 7.104, 0.489, 1.049, 1.327, 2.077])


    def externalEnergyFunction(self, trial):
        # inverse exponential with offset, y = a * exp(b/x) + c
        predicted = trial[0] * numpy.exp(trial[1] / self.xData) + trial[2]

        # sum of squared error
        error = predicted - self.yData
        return numpy.sum(error*error)



if __name__ == '__main__':
    solver = TestSolver(3, 600, 600, -10, 10, "Rand2Exp", 0.7, 0.6, 0.1, True)
    while solver.exitCode ==0:
        solver.computeGeneration()
    solver.finalize()
    #solver.Solve()
    
    print
    predicted = solver.bestSolution[0] * numpy.exp(solver.bestSolution[1] / solver.xData) + solver.bestSolution[2]
    print "yData yfit"
    for i in range(len(solver.yData)):
        print solver.yData[i], predicted[i]
