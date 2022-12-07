import math
import itertools

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from cycler import cycler

from stimator.utils import _is_string, _is_sequence


def _prepare_settigs(style, palette, font, fig_size):
    st_list = []

    more_custom_settings = {  # 'lines.markersize': 5,
                            'lines.markeredgewidth': 0.1}

    if style is not None:
        if not _is_sequence(style):
            style = [style]
        st_list.extend(style)
    else:
        st_list.extend(['seaborn-whitegrid', 'seaborn-notebook'])

    valid_styles = []
    for s in st_list:
        if not _is_string(s):
            valid_styles.append(s)
        else:
            if s in mpl.style.available:
                valid_styles.append(s)
    st_list = valid_styles

    if fig_size is not None:
        more_custom_settings['figure.figsize'] = fig_size
    else:
        more_custom_settings['figure.figsize'] = (8, 5.5)

    if palette is not None:
        more_custom_settings['axes.prop_cycle'] = cycler('color', list(palette))

    more_custom_settings['font.family'] = font

    st_list.append(more_custom_settings)

    return st_list

def _get_indexes_to_plot(timecourse, what):
    if what is None:
        return range(len(timecourse.names))
    if _is_string(what):
        what = [what]
    if not _is_sequence(what):
        raise ValueError(f"No data for {str(what)} in timecourse")
    indexes2plot = []
    for elem in what:
        if isinstance(elem, int):
            if (0 <= elem < len(timecourse.names)) and (elem not in indexes2plot):
                indexes2plot.append(elem)
            else:
                raise ValueError(f"Index {elem} out of range in timecourse")
        elif _is_string(elem):
            try:
                i = timecourse.names.index(elem)
            except ValueError:
                raise ValueError(f"No data for {str(elem)} in timecourse")
            if i not in indexes2plot:
                indexes2plot.append(i)
    return indexes2plot




def plot_timecourse(timecourse, what=None, ax=None, title=None, **kwargs):
    if ax is None:
        ax = plt.gca()

    indexes_to_plot = _get_indexes_to_plot(timecourse, what)
    nseries = len(indexes_to_plot)

    # cyl = [c['color'] for c in mpl.rcParams['axes.prop_cycle']]
    # cyclingcolors = itertools.cycle(cyl)
    for i in indexes_to_plot:
        name = timecourse.names[i]
        y = timecourse[i]
        data_loc = np.logical_not(np.isnan(y))
        x = timecourse.t[data_loc]
        y = y[data_loc]
        ax.plot(x, y, label=name, **kwargs)
    if title is not None:
        ax.set_title(title)
    elif timecourse.title:
        ax.set_title(timecourse.title)
    handles, lbls = ax.get_legend_handles_labels()
    ax.legend(handles, lbls, loc='best')
    return ax


#     for line in lines_desc:
#         sol = solutions[line['solution_index']]
#         y = sol[line['var_index']]
#         data_loc = np.logical_not(np.isnan(y))
#         x = sol.t[data_loc]
#         y = y[data_loc]
#         ax.plot(x, y, ls=ls, marker=marker,
#                 color=line['color'], label=line['name'],
#                 clip_on=False)
#     ax.set_title(title)


def _plotTC(lines_desc, solutions, title, ls, marker, ax):
    for line in lines_desc:
        sol = solutions[line['solution_index']]
        y = sol[line['var_index']]
        data_loc = np.logical_not(np.isnan(y))
        x = sol.t[data_loc]
        y = y[data_loc]
        ax.plot(x, y, ls=ls, marker=marker,
                color=line['color'], label=line['name'],
                clip_on=False)
    ax.set_title(title)


def plotTCs(solutions,
            show=False,
            figure=None,
            axis_set=None,
            fig_size=None,
            titles=None,
            ynormalize=False,
            yrange=None,
            group=False,
            suptitlegend=None,
            legend=True,
            force_dense=False,
            style=None,
            palette=None,
            font="sans-serif",
            save2file=None, **kwargs):

    """Generate a graph of the time course using matplotlib.

       Called by .plot() member function of class timecourse.Solutions"""

    settings = _prepare_settigs(style, palette, font, fig_size)

    with plt.style.context(settings):

        # handle names and titles
        nsolutions = len(solutions)
        pnames = ['time course %d' % (i+1) for i in range(nsolutions)]
        for i in range(nsolutions):
            if titles:
                pnames[i] = titles[i]
            else:
                if solutions[i].title:
                    pnames[i] = solutions[i].title

        # find how many plots
        nplts = len(group) if group else nsolutions

        # compute shape of grid
        ncols = int(math.ceil(math.sqrt(nplts)))
        nrows = int(math.ceil(float(nplts)/ncols))

        # handle axes
        if axis_set is None:
            figure, axis_set = plt.subplots(nrows, ncols, figsize=fig_size)
            if hasattr(axis_set, 'shape'):
                axis_set = list(axis_set.flat)
            else:
                axis_set = [axis_set]

        cyl = [c['color'] for c in mpl.rcParams['axes.prop_cycle']]
        cyclingcolors = itertools.cycle(cyl)

        # create "plot description" records
        plots_desc = []
        color_table = {}
        if not group:
            for k, sol in enumerate(solutions):
                line_desc = []
                for i in range(len(sol)):
                    name = sol.names[i]
                    line = {'name': name,
                            'solution_index': k,
                            'var_index': i}
                    if name in color_table:
                        c = color_table[name]
                    else:
                        c = next(cyclingcolors)
                        color_table[name] = c
                    line['color'] = c
                    line_desc.append(line)

                plots_desc.append({'title': pnames[k], 'lines': line_desc})
        else:
            for g in group:
                line_desc = []
                if _is_string(g):
                    pdesc = {'title': g}
                    for k, sol in enumerate(solutions):
                        if g in sol.names:
                            indx = sol.names.index(g)
                            line = {'name': pnames[k],
                                    'solution_index': k,
                                    'var_index': indx}
                            line_desc.append(line)
                else:
                    if not _is_sequence(g):
                        raise StimatorTCError('%s is not a str or seq' % str(g))

                    pdesc = {'title': ' '.join(g)}
                    for vvv in g:
                        for k, sol in enumerate(solutions):
                            if vvv in sol.names:
                                indx = sol.names.index(vvv)
                                line = {'solution_index': k, 'var_index': indx}
                                if len(solutions) > 1:
                                    line['name'] = "%s, %s" % (vvv, pnames[k])
                                else:
                                    line['name'] = "%s" % (vvv)
                                line_desc.append(line)
                for line in line_desc:
                    unique_id = (line['solution_index'], line['var_index'])
                    if unique_id in color_table:
                        c = color_table[unique_id]
                    else:
                        c = next(cyclingcolors)
                        color_table[unique_id] = c
                    line['color'] = c

                pdesc['lines'] = line_desc
                plots_desc.append(pdesc)

        # draw plots
        for i, p in enumerate(plots_desc):
            curraxis = axis_set[i]
            use_dots = not solutions[0].dense
            if force_dense:
                use_dots = False

            ls, marker = ('None', 'o') if use_dots else ('-', 'None')

            _plotTC(p['lines'], solutions, p['title'], ls, marker, curraxis)

            if yrange is not None:
                curraxis.set_ylim(yrange)
            if legend:
                handles, lbls = curraxis.get_legend_handles_labels()
                curraxis.legend(handles, lbls, loc='best')
            curraxis.set_xlabel('')
            curraxis.set_ylabel('')

        # draw suptitle (needs a figure object)
        fig_obj = plt.gcf()
        if suptitlegend is not None:
            fig_obj.suptitle(suptitlegend)
        elif hasattr(solutions, 'title'):
            fig_obj.suptitle(solutions.title)

        if ynormalize and not yrange:
            rs = [a.get_ylim() for a in axis_set]
            common_range = min([l for l, h in rs]), max([h for l, h in rs])
            for a in axis_set:
                a.set_ylim(common_range)

        # pl.tight_layout()

        if save2file is not None:
            figure.savefig(save2file)
        if show:
            if save2file is not None:
                if hasattr(save2file, 'read'):
                    save2file.close()
            plt.show()


def plot_estim_optimum(opt, figure=None,
                       axis_set=None,
                       fig_size=None,
                       style=None,
                       palette=None,
                       font="sans-serif",
                       save2file=None,
                       show=False):

    settings = _prepare_settigs(style, palette, font, fig_size)

    with plt.style.context(settings):
        if axis_set is None:
            if figure is None:
                figure = plt.figure()

        bestsols = opt.optimum_dense_tcs
        expsols = opt.optimizer.tc
        tcstats = opt.tcdata
        nplts = len(bestsols)
        ncols = int(math.ceil(math.sqrt(nplts)))
        nrows = int(math.ceil(float(nplts)/ncols))
        if axis_set is None:
            axis_set = [figure.add_subplot(nrows, ncols, i + 1) for i in range(nplts)]
        else:
            axis_set = axis_set

        for i in range(nplts):
            subplot = axis_set[i]
            # subplot.set_xlabel("time")
            subplot.set_title("%s (%d pt) %g" % tcstats[i], fontsize=12)
            expsol = expsols[i]
            symsol = bestsols[i]

            cyl = [c['color'] for c in mpl.rcParams['axes.prop_cycle']]
            cyclingcolors = itertools.cycle(cyl)

            for line in range(len(expsol)):
                # count NaN and do not plot if they are most of the timecourse
                yexp = expsol[line]
                nnan = len(yexp[np.isnan(yexp)])
                if nnan >= expsol.ntimes-1:
                    continue
                # otherwise plot lines
                xname = expsol.names[line]
                ysim = symsol[symsol.names.index(xname)]
                lsexp, mexp = 'None', 'o'
                lssim, msim = '-', 'None'

                color = next(cyclingcolors)

                subplot.plot(expsol.t, yexp,
                             marker=mexp, ls=lsexp, color=color, clip_on=False)
                subplot.plot(symsol.t, ysim,
                             marker=msim, ls=lssim, color=color,
                             label='%s' % xname, clip_on=False)
            subplot.legend(loc='best')

        if save2file is not None:
            figure.savefig(save2file)
        if show:
            if save2file is not None:
                if hasattr(save2file, 'read'):
                    save2file.close()
            plt.show()


def plot_generations(opt, generations=None,
                     pars=None, figure=None, show=False,
                     fig_size=None,
                     style=None, palette=None, font="sans-serif"):

    if not opt.generations_exist:
        raise IOError('file generations.txt was not generated')

    settings = _prepare_settigs(style, palette, font, fig_size)
    settings.append({'lines.markeredgewidth': 1.0})

    with plt.style.context(settings):

        if figure is None:
            figure = plt.figure()
        figure.clear()

        if generations is None:
            all_gens = list(range(opt.optimization_generations + 1))
            dump_generations = all_gens

        if pars is None:
            first2 = opt.parameters[:2]
            pars = [p[0] for p in first2]

        pnames = [p[0] for p in opt.parameters]

        colp0 = pnames.index(pars[0])
        colp1 = pnames.index(pars[1])

        scores_col = len(opt.parameters)

        # ax1 = pl.subplot(1,2,1)
        # ax2 = pl.subplot(1,2,2)
        ax2 = plt.subplot(1, 1, 1)

        # parse generations
        gen = -1
        f = open('generations.txt')
        solx = []
        soly = []
        objx = []
        objy = []
        reading = False
        for line in f:
            line = line.strip()
            if line == '' and reading:
                if len(solx) > 0:
                    # ax1.plot(solx, soly, marker='o', ls='None', label=gen)
                    # for px, py in zip(objx, objy):
                    #     ax2.axhline(py, xmin=px, xmax=px*0.01, color='black')
                    ax2.axvline(objx[0], lw=0.5, color='lightgray')
                    ax2.plot(objx, objy, '_', ls='None', label=gen)
                    solx = []
                    soly = []
                    objx = []
                    objy = []
                    reading = False
            elif line.startswith('generation'):
                gen = line.split()[1]
                igen = int(gen)
                if igen in dump_generations:
                    reading = True
                    # print 'generation', gen
            elif reading:
                line = [float(x) for x in line.split()]
                solx.append(line[colp0])
                soly.append(line[colp1])
                objx.append(igen)
                objy.append(line[scores_col])
            else:
                continue
        f.close()
        # ax1.legend(loc=0)
        # ax1.set_title('population')
        # ax1.set_xlabel(pars[0])
        # ax1.set_ylabel(pars[1])
        ax2.set_title('scores')
        ax2.set_yscale('log')
        ax2.set_xlabel('generation')
        if show:
            plt.show()

# ----------------------------------------------------------------------------
#         TESTING CODE
# ----------------------------------------------------------------------------


if __name__ == "__main__":
    from pathlib import Path
    from stimator.timecourse import Solution, Solutions, readTCs
    import stimator.examples.timecourses as expl_tcs

    demodata = """
#this is demo data with a header
t x y z
0       0.95 0         0
0.1                  0.1

  0.2 skip 0.2 skip this
nothing really usefull here
- 0.3 0.3 this line should be skipped
#0.4 0.4
0.3 0.4 0.5 0.55
0.4 0.5 0.6 0.7
0.5 0.6 0.8 0.9
0.55 0.7 0.85 0.95
0.6  - 0.5 - -

"""

    demodata2 = """
#this is demo data with a header
t x y z
0       0.95 0         0
0.1                  0.09

  0.2 skip 0.2 skip this
nothing really usefull here
- 0.3 0.3 this line should be skipped
#0.4 0.4
0.3 0.45 0.55 0.58
0.4 0.5 0.65 0.75
0.5 0.65 0.85 0.98
0.55 0.7 0.9 1.45
0.6  - 0.4 - -
"""
    sol = Solution().read_str(demodata)
    sol2 = Solution().read_str(demodata2)

    print(sol)

    example_tc = Path(expl_tcs.__path__[0]) / 'TSH2b.txt'

    stsh = Solution()
    stsh.read_from(example_tc)

    plot_timecourse(sol, title="simple TC plot")
    plt.show()
    plot_timecourse(sol, title="plot with different lw", linewidth=10)
    plt.show()
    plot_timecourse(stsh, title="plot with different marker", marker="+")
    plt.show()
    plot_timecourse(sol, what='z', title="plot with filtering")
    plt.show()
    plot_timecourse(sol, what=[0, 'z'], title="plot with filtering")


    # sols = Solutions([Solution(title='the first tc').read_str(demodata),
    #                   Solution().read_str(demodata2)],
    #                  title='all time courses')

    # sols.plot(suptitlegend="plotting the two time courses")
    # sols.plot(suptitlegend="with font=serif, palette='rgb'",
    #           font_scale=1.3, font='serif', palette='rgb')
    # p = ['crimson', 'mediumpurple', 'darkolivegreen']
    # sols.plot(suptitlegend=f"with font=serif, palette = {p}",
    #           font_scale=0.5, font='serif', palette=p)
    # sols.plot(suptitlegend="with style=default", style='default')
    # sols.plot(suptitlegend="with style=seaborn-darkgrid", style='seaborn-darkgrid')
    # sols.plot(suptitlegend="with style=bogus", style='bogus')

    # sols.plot(fig_size=(12, 6), suptitlegend="with fig_size=(12,6)")

    # sols.plot(suptitlegend="with force_dense=True", force_dense=True)
    # sols.plot(ynormalize=True, suptitlegend='with ynormalize=True')
    # sols.plot(yrange=(0, 2), suptitlegend='with yrange=(0,2)')

    # sols.plot(group=['z', 'x'], suptitlegend="with group=['z', 'x']")
    # sols.plot(group=['z', ('x', 'y')], suptitlegend="with group=['z', ('x','y')]")
    # sols.plot(group=['z', ('x', 'z')], suptitlegend="with group=['z', ('x','z')]")

    # f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    # sols.plot(suptitlegend="with given axis_set",
    #           force_dense=True,
    #           axis_set=[ax1, ax2])
    # ax1.set_ylabel('variables')
    # ax2.set_ylabel('variables')
    # ax2.set_xlabel('time')

    # sol = Solution().read_str(demodata)
    # sol.plot(group=['z', 'x'], suptitlegend="1 tc with group=['z', 'x']")
    # sol.plot(group=['z', ('x', 'y')],
    #          suptitlegend="1tc with group=['z', ('x','y')]")
    # sol.plot(group=['z', ('x', 'z')],
    #          suptitlegend="1tc with group=['z', ('x','z')]")
    
    # from pathlib import Path
    # # print(Path.cwd())
    # # print(Path(__file__).resolve())
    # examples_loc = str(Path(__file__).resolve().parent / 'examples' /'timecourses')
    # #print(locfile)

    # # sol.read_from('examples/timecourses/TSH2b.txt')
    # sol.read_from(examples_loc + '/TSH2b.txt')

    # f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    # sol.plot(suptitlegend="plotting on a given axes (1 TC)", axes=ax2)
    # ax2.set_ylabel('concentrations')
    # ax2.set_xlabel('time')

    # print('\n!! testing transformations ----------------')

    # sols = Solutions(title='all time courses')
    # s = Solution(title='original time course').read_str(demodata2)
    # sols += s

    # def average(x, t):
    #     # print ('applying transformation')
    #     return np.array([t/2.0, (x[0]+x[-1])/2.0])

    # s = s.transform(average,
    #                 newnames=['t/2', 'mid point'],
    #                 new_title='after transformation')
    # sols += s

    # sols.plot(suptitlegend="plotting original and transf", force_dense=True)

    # tcs = readTCs(['TSH2b.txt', 'TSH2a.txt'],
    #               examples_loc,
    #               names="SDLTSH HTA".split(),
    #               verbose=False)
    # tcs.plot(suptitlegend="read from file")
    # tcs.plot(group=['SDLTSH'], suptitlegend="read from file with group=['SDLTSH']")

    plt.show()
