from __future__ import print_function, absolute_import, division

import time

from numpy import array, nansum, fabs, copy, empty, linspace, isnan
from scipy import integrate
import numpy as np

import stimator.de as de
from stimator.dynamics import getdXdt, init2array, ModelSolver, solve, get_outputs_decl
import stimator.fim as fim
import stimator.timecourse as timecourse
import stimator.plots as plots


# ----------------------------------------------------------------------------
#         Optimum criteria for TC similarity
# ----------------------------------------------------------------------------

def get_matching_indexes(repnames, tc):
    # find indexers to match a solution and tc
    # to create compatible subsets for objetive functions
    # return the subsets
    sol_indexes, tc_indexes = [], []
    ntimes = tc.ntimes

    for i, yexp in enumerate(tc):
        # count NaN, keep only columns with sufficient number of points
        nnan = len(yexp[np.isnan(yexp)])
        if nnan >= ntimes - 1:
            continue
        colname = tc.names[i]
        if colname not in repnames:
            continue
        tc_indexes.append(i)
        sol_indexes.append(repnames.index(colname))

    tc_indexes = array(tc_indexes, int)
    sol_indexes = array(sol_indexes, int)
    return sol_indexes, tc_indexes


def getRangeVars(tcs, varnames):
    ranges = [0.0 for i in range(len(varnames))]
    for ix, x in enumerate(varnames):
        for tc in tcs:
            yexp = tc[x]
            tpe = (max(yexp) - min(yexp))
            ranges[ix] = max(ranges[ix], tpe)
    return ranges

def get_criterium(repnames, tc, weights=None):
    """Returns a function to compute the objective function given a solution and a timecourse.

    The function has signature
    criterium(Y,i)
    Y is the predicted timecourse, for a given set of parameters.
    i is the index of the timecourse.
    The function returns a float.

    tc is a Solutions object holding ('experimental') timecourse data,
    each timecourse has shape (nvars, ntimes).

    weights can be:

    None         : no weighting (simple LSquares, S = sum((Ypred-Yexp)**2))
    all others are weighted LSquares, S = (Ypred-Yexp).T * W * (Ypred-Yexp)
    'demo'       : demo weighting  W = 1/j with j = 1,...,nvars
    """
    sol_indexes, tc_indexes = get_matching_indexes(repnames, tc)

    def unweighted_criterium(sol, tc):
        sol_subset, tc_subset = sol[sol_indexes], tc[tc_indexes]
        d = (sol_subset - tc_subset)
        return np.sum(d * d)

    if weights is None:
        return unweighted_criterium
    else:
        return unweighted_criterium

        #     return np.sum(d.T * W[i] * d)

    # TODO: weights not implemented
    return None


# ----------------------------------------------------------------------------
#         Class to perform DE optimization for ODE systems
# ----------------------------------------------------------------------------

class OptimumData(object):
    """Object that holds optimum solution data."""

    def __init__(self, optimizer):
        self.optimizer = optimizer

    def info(self):
        optimum = self
        headerformat = "--- %-20s -----------------------------\n"
        res = "\n" + (headerformat % 'PARAMETERS')
        res += "\n".join(["%s\t%12g +- %g" % i for i in optimum.parameters])
        res += '\n\n'
        res += headerformat % 'OPTIMIZATION'
        res += "%s\t%g\n" % ('Final Score', optimum.optimization_score)
        res += "%s\t%d\n" % ('generations', optimum.optimization_generations)

        res += "%s\t%d\n" % ('max generations', optimum.max_generations)
        res += "%s\t%d\n" % ('population size', optimum.pop_size)
        res += "%s\t%s\n" % ('Exit by',     optimum.optimization_exit_by)
        res += '\n\n'

        res += headerformat % 'TIME COURSES'
        res += '\t\t'.join(['Name', 'Points', 'Score'])+'\n'
        res += "\n".join(["%s\t%d\t%g" % i for i in optimum.tcdata])
        res += '\n\n'
        return res

    def print_info(self):
        print (self.info())
    
    def __str__(self):
        return self.info()

    def plot(self, **kwargs):
        return plots.plot_estim_optimum(self, **kwargs)

    def plot_generations(self, **kwargs):
        return plots.plot_generations(self, **kwargs)

class DeODEOptimizer(de.DESolver):
    """Overides energy function and report functions.

    The energy function solves ODEs and computes a least-squares score.
    Ticker functions are called on completion of a generation and when
    optimization finishes.
    """

    def __init__(self, model, optSettings, tcs, weights=None,
                 aMsgTicker=None,
                 anEndComputationTicker=None,
                 dump_generations=None,
                 dump_predictions=False,
                 initial='init',
                 max_generations=200,
                 convergence_noimprovement=20):
        #self.varnames = model.varnames
        self.model = model.copy()
        self.model_solvers = []
        self.tc = tcs
        self.endTicker = anEndComputationTicker
        self.msgTicker = aMsgTicker
        self.dump_predictions = dump_predictions
        self.dump_generations = dump_generations

        # reorder variables according to model
        self.tc.order_by_modelvars(self.model)

        pars = model.with_bounds
        self.par_names = [p.name for p in self.model.with_bounds]
        self.varnames = list(self.model.varnames)
        mins = array([u.bounds.lower for u in pars])
        maxs = array([u.bounds.upper for u in pars])

        if optSettings.get('pop_size', None) is None:
            optSettings['pop_size'] = optSettings['genomesize']
        if optSettings.get('max_generations', None) is None:
            optSettings['max_generations'] = optSettings['generations']
        max_generations=optSettings['max_generations']

        de.DESolver.__init__(self, len(pars),  # number of parameters
                             int(optSettings['pop_size']),  # pop size
                             mins, maxs,  # min and max parameter values
                             "Best2Exp",  # DE strategy
                             0.7, 0.6, 0.0, # DiffScale, p crossover, Cut-off S
                             max_generations=max_generations,
                             convergence_noimprovement=convergence_noimprovement)

        # cutoffEnergy is 1e-6 of deviation from data
        self.cutoffEnergy = 1.0e-6 * sum([nansum(fabs(tc.data)) for tc in self.tc])

        # create one ModelSolver per timecorse
        self.criterium = []

        # overide initial values of solver with first time-course point
        for tc in self.tc:
            X0 = []
            for xname in model.varnames:
                x0value = tc[xname][0] if xname in tc.names else self.model.get_init(xname)
                X0.append(x0value)
            X0 = array(X0, dtype=float)
            ms = ModelSolver(self.model,
                            times=tc.t.copy(),
                            initial=X0,
                            title=tc.title,
                            changing_pars=self.par_names)
            self.model_solvers.append(ms)
            solnames = ms.solutions_names()
            self.criterium.append(get_criterium(solnames, tc, weights))

        self.timecourse_scores = empty(len(self.tc))

    def computeSolution(self, i, trial, dense=None):
        """Computes solution for timecourse i, given parameters trial."""

        sol = self.model_solvers[i].solve(par_values=trial)
        return sol

    def external_score_function(self, trial):
        # if out of bounds flag with error energy
        for p, minInitialValue, maxInitialValue in zip(trial, self.min_values, self.max_values):
            if p > maxInitialValue or p < minInitialValue:
                return float('inf')
        # set parameter values from trial
        #self.model.set_uncertain(trial)

        # compute solutions and scores
        for i in range(len(self.tc)):
            sol = self.model_solvers[i].solve(par_values=trial)
            #sol = self.computeSolution(i, trial)
            if sol is not None:
                self.timecourse_scores[i] = self.criterium[i](sol, self.tc[i])
            else:
                return float('inf')

        globalscore = self.timecourse_scores.sum()
        return globalscore

    def reportInitial(self):
        msg = "\nSolving %s..." % self.model.metadata.get('title', '')
        #initialize stopwatch
        self.start_time = time.time()
        if self.dump_generations is not None:
            self.dumpfile = open('generations.txt', 'w')
        if not self.msgTicker:
            print (msg)
        else:
            self.msgTicker(msg)

    def reportGeneration(self):
        msg = "%-4d: %f" % (self.generation, float(self.best_score))
        if not self.msgTicker:
            print (msg)
        else:
            self.msgTicker(msg)
        if self.dump_generations is not None:
            print (self.generation_string(self.generation), file=self.dumpfile)

    def reportFinal(self):
        if self.exitCode <= 0:
            outCode = -1
        else:
            outCode = self.exitCode
            self.generate_optimum()
        if not self.endTicker:
            de.DESolver.reportFinal(self)
        else:
            self.endTicker(outCode)
        if self.dump_generations is not None:
            print (self.generation_string(self.generation), file=self.dumpfile)
            self.dumpfile.close()

    def generation_string(self, generation):
        generation = str(generation)
        # find if objectives is iterable
        isiter = hasattr(self.scores[0], '__contains__')
        res = 'generation %s -------------------------\n' % generation
        for s, o in zip(self.pop, self.scores):
            sstr = ' '.join([str(i) for i in s])
            if isiter:
                ostr = ' '.join([str(i) for i in o])
            else:
                ostr = str(o)
            res = res + '%s %s\n' % (sstr, ostr)
        return res

    def generate_optimum(self):
        # compute parameter standard errors, based on FIM-1
        # generate TC solutions
        best = OptimumData(self)
        best.optimization_score = self.best_score
        best.optimization_generations = self.generation
        best.optimization_exit_by = self.exitCodeStrings[self.exitCode]
        best.max_generations = self.max_generations
        best.pop_size = self.pop_size

        # TODO: Store initial solver parameters?

        # generate best time-courses

        parameters = list(zip(self.par_names, [x for x in self.best]))

        sols = timecourse.Solutions()
        best.tcdata = []

        for (i, tc) in enumerate(self.tc):
            sol = self.computeSolution(i, self.best)
            if sol is not None:
                score = self.criterium[i](sol, tc)
            else:
                score = 1.0E300
            sols += sol
            best.tcdata.append((self.tc[i].title, tc.ntimes, score))

        best.optimum_tcs = sols

        if not (fim.SYMPY_INSTALLED):
            best.parameters = [(p, v, 0.0) for (p, v) in parameters]
        else:
            commonvnames = self.tc.get_common_full_vars()
            commonvnames = set(commonvnames).intersection(set(self.model.varnames))
            if len(commonvnames) == 0:
                commonvnames = self.model.varnames
                # consterror = getRangeVars(best.optimum_tcs, commonvnames)
                consterror = None
            else:
                commonvnames = list(commonvnames)
                # print(commonvnames)
                consterror = getRangeVars(self.tc, commonvnames)
                # assume 5% of range
                consterror = timecourse.constError_func([r * 0.05 for r in consterror])

            try:
                FIM1, invFIM1 = fim.computeFIM(self.model,
                                            parameters,
                                            sols,
                                            consterror,
                                            commonvnames)
                best.parameters = [(self.par_names[i],
                                    value,
                                    invFIM1[i, i]**0.5)
                                    for (i, value) in enumerate(self.best)]
            except:
                best.parameters = [(p, v, 0.0) for (p, v) in parameters]

        sols = timecourse.Solutions()
        for (i, tc) in enumerate(self.tc):
            sol = self.computeSolution(i, self.best, dense=True)
            # ts = linspace(tc.t[0], tc.t[-1], 500)

            # sol = timecourse.SolutionTimeCourse(ts, Y.T,
            #                                     self.varnames,
            #                                     title=tc.title)
            sols += sol

        best.optimum_dense_tcs = sols

        if self.dump_generations is not None:
            best.generations_exist = True
        else:
            best.generations_exist = False

        self.optimum = best
        # self.generate_fitted_sols()

        if self.dump_predictions:
            fnames = ['pred_' + self.tc[i].title for i in range(len(self.tc))]
            best.optimum_tcs.write_to(fnames, verbose=True)


def s_timate(model, timecourses=None, opt_settings=None,
             tc_dir=None,
             names=None,
             verbose_readingTCs=True,
             **kwargs):

    # create a default dict of optimizer settings,
    # then update with .metadata['optSettings']
    # finally, update with argument opt_settings
    optSettings = {'pop_size': 80,
                   'max_generations': 200,
                   'optimizer': 'DeODEOptimizer'}
    if model.metadata.get('optSettings', None) is not None:
        optSettings.update(model.metadata['optSettings'])
    if opt_settings is not None:
        optSettings.update(opt_settings)

    # timecourses argument is used to indicate time-course files
    # if it is None, then use model.metadata['timecourses']
    if timecourses is None:
        timecourses = model  # use model as source in readTCs

    tcs = timecourse.readTCs(timecourses,
                             filedir=tc_dir,
                             names=names,
                             verbose=verbose_readingTCs)

    optimizer = DeODEOptimizer(model, optSettings, tcs, **kwargs)
    optimizer.run()
    return optimizer.optimum


def test():
    from stimator import read_model
    m1 = read_model("""
title Glyoxalase system in L. Infantum

glx1 : HTA -> SDLTSH, V1*HTA/(Km1 + HTA)
#glx1 : HTA -> SDLTSH, V*HTA/(Km1 + HTA), V=2.57594e-05
glx2 : SDLTSH ->,     V2*SDLTSH/(Km2 + SDLTSH)

#find glx1.V  in [0.00001, 0.0001]
find V1  in [0.00001, 0.0001]

Km1 = 0.252531
find Km1 in [0.01, 1]

V2  = 2.23416e-05
find V2 in [0.00001, 0.0001]

Km2 = 0.0980973
find Km2 in (0.01, 1)

init : (SDLTSH = 7.69231E-05, HTA = 0.1357)

timecourse TSH2a.txt
timecourse TSH2b.txt
""")

    # print m1
    # get tcdir
    import pathlib
    this_file = pathlib.Path(__file__)
    tcdir = pathlib.Path(this_file.parents[0], 'examples', 'timecourses')
    #print(tcdir)

    optimum = s_timate(m1, tc_dir=tcdir, #timecourses=['TSH2a.txt', 'TSH2b.txt'],
                       names=['SDLTSH', 'HTA'],
                       dump_generations=True) 
    # convergence_noimprovement=40)
    # ... intvarsorder=(0,2,1) ...

    print(optimum)
    optimum.plot()
    optimum.plot_generations(pars=['V2', 'Km1'], fig_size=(9,6))


    # --- an example with transformations --------------------

    mtransf = read_model("""
title Glyoxalase system in L. Infantum

glx1 : HTA -> SDLTSH, V1*HTA/(Km1 + HTA)
#glx1 : HTA -> SDLTSH, V*HTA/(Km1 + HTA), V=2.57594e-05
glx2 : SDLTSH ->,     V2*SDLTSH/(Km2 + SDLTSH)

#find glx1.V  in [0.00001, 0.0001]
find V1  in [0.00001, 0.0001]

Km1 = 0.252531
find Km1 in [0.01, 1]

V2  = 2.23416e-05
find V2 in [0.00001, 0.0001]

Km2 = 0.0980973
find Km2 in (0.01, 1)

~sdlx2 = 2 * SDLTSH

!! sdlx2

init : (SDLTSH = 7.69231E-05, HTA = 0.1357)

""")

    optimum = s_timate(mtransf, tc_dir=tcdir, timecourses=['tc_double.txt'],
                       names=['sdlx2', 'SDLTSH', 'HTA']) 

    print(optimum)
    optimum.plot()


    # --- an example with unknown initial values --------------------

    m2 = m1.copy()

    # Now, assume init.HTA is uncertain
    m2.set_bounds('init.HTA', (0.05, 0.25))
    # do not estimate Km1 and Km2, to help the analysis
    m2.reset_bounds('Km1')
    m2.reset_bounds('Km2')

    # VERY IMPORTANT:
    # only one time course can be used:
    # cannot fit one initial value using several timecourses!!!

    optimum = s_timate(m2, timecourses=['TSH2a.txt'], 
                       tc_dir=tcdir,
                       opt_settings={'pop_size': 60},
                       names=['SDLTSH', 'HTA'])

    print(optimum)
    optimum.plot(show=True)

if __name__ == "__main__":
    test()
