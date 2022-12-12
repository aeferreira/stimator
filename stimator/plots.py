import math
import itertools

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from cycler import cycler, Cycler

from stimator.utils import _is_string, _is_sequence


def _get_cols_to_plot(timecourse, what):
    if what is None:
        return timecourse.names[:]
    if _is_string(what):
        what = [what]
    if not _is_sequence(what):
        raise ValueError(f"No data for {str(what)} in timecourse")
    cols2plot = []
    for elem in what:
        if isinstance(elem, int):
            if (0 <= elem < len(timecourse.names)):
                name = timecourse.names[elem]
                if name not in cols2plot:
                    cols2plot.append(name)
            else:
                raise ValueError(f"Index {elem} out of range in timecourse")
        elif _is_string(elem):
            if elem not in timecourse.names:
                raise ValueError(f"No data for {str(elem)} in timecourse")
            if elem not in cols2plot:
                cols2plot.append(elem)
    return cols2plot


def _get_overiding_prop_table(styling, prop_cycle):
    if styling is None:
        return None
    if not isinstance(styling, dict):
        raise TypeError("'styling' parameter must be a dict")

    props = list(prop_cycle)

    stylingdict = {}

    for k, v in styling.items():
        if isinstance(v, dict):
            stylingdict[k] = v
        elif isinstance(v, int):
            stylingdict[k] = props[v]
        elif _is_string(v):
            stylingdict[k] = {'color': v}
        else:
            raise ValueError(f"'{v}' is an invalid color or style ")
    return stylingdict


MPL_QUALIT = ['Pastel1', 'Pastel2', 'Paired', 'Accent',
              'Dark2', 'Set1', 'Set2', 'Set3', 'tab10',
              'tab20', 'tab20b', 'tab20c']

MPL_QUALITATIVE_TABLE = {cm: mpl.colormaps[cm].colors for cm in MPL_QUALIT}


def _get_list_from_palette(pname):
    """get a list of RGB values from a palette name.

    mpl qualitative names are accepted.
    TODO: accept seaborn and palettable generated palettes as pname"""
    if pname in MPL_QUALIT:
        return MPL_QUALITATIVE_TABLE[pname]
    raise ValueError(f'palette name {pname} not found')


def _get_new_prop_cycle(prop_cycle, palette):
    if prop_cycle is None and palette is None:
        return None
    if prop_cycle is not None:
        if isinstance(prop_cycle, dict):
            prop_cycle = cycler(**prop_cycle)
        if not isinstance(prop_cycle, Cycler):
            raise TypeError("'prop_cycle' must be a dict or a Cycler object")

    # include the palette in the prop_cycler
    if palette is not None:
        color_list = _get_list_from_palette(palette)
        color_cycle = cycler('color', color_list)
        if prop_cycle is None:
            return color_cycle

        # merge with specified cycler
        color_list = []
        # get color list as long as the prop_cycle
        for _, p in zip(range(len(prop_cycle)), color_cycle()):
            color_list.append(p['color'])

        # substitute colors
        trans = prop_cycle.by_key()
        trans['color'] = color_list
        prop_cycle = cycler(**trans)
    return prop_cycle


def plot_timecourse(timecourse, what=None, ax=None,
                    prop_cycle=None, palette=None, styling=None,
                    title=None, tight_t0=True, **kwargs):
    if ax is None:
        ax = plt.gca()

    # locally set the prop_cycler if prop_cycle or palette are given
    used_prop_cycle = _get_new_prop_cycle(prop_cycle, palette)
    if used_prop_cycle is not None:
        ax.set_prop_cycle(used_prop_cycle)
    else:
        used_prop_cycle = mpl.rcParams['axes.prop_cycle']

    # find what (which columns) to plot
    cols_to_plot = _get_cols_to_plot(timecourse, what)

    # plot (as lines)
    for name in cols_to_plot:
        x = timecourse.t
        y = timecourse[name]
        ax.plot(x, y, label=name, **kwargs)

    # overide prop_cycler styles if these are specified on the styling dict
    # for some variables
    if styling is not None:
        stylingdict = _get_overiding_prop_table(styling, used_prop_cycle)
        lines = ax.get_lines()
        for i, name in enumerate(cols_to_plot):
            if name in stylingdict:
                lines[i].set(label=name, **(stylingdict[name]))

    # draw legend
    handles, lbls = ax.get_legend_handles_labels()
    ax.legend(handles, lbls, loc='best')

    # draw plot title
    if title is not None:
        ax.set_title(title)
    elif timecourse.title:
        ax.set_title(timecourse.title)

    # remove the left margin for the x-axis
    # to start the plot at t0
    if tight_t0:
        if timecourse.t[0] < timecourse.t[-1]:
            ax.set_xlim(timecourse.t[0], None)
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
0.1   0.1  0.5        0.2

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

    # print(sol)

    example_tc = Path(expl_tcs.__path__[0]) / 'TSH2b.txt'

    stsh = Solution()
    stsh.read_from(example_tc)

    plot_timecourse(sol, title="simple TC plot")
    plt.show()
    plot_timecourse(sol,
                    title="plot with different lw, linewidth=10",
                    linewidth=10)
    plt.show()
    plot_timecourse(stsh,
                    title="plot with different marker, marker='+'",
                    tight_t0=False, marker='+')
    plt.show()
    plot_timecourse(sol, what='z',
                    title="plot with filtering, what='z'")
    plt.show()
    plot_timecourse(sol, what=[0, 'z'],
                    title="plot with filtering, what=[0, 'z']")
    plt.show()

    with plt.style.context('ggplot'):
        plot_timecourse(sol, title="simple TC plot with ggplot style")
    plt.show()

    custom_cycler = (cycler(color=['c', 'm', 'y', 'k']) +
                     cycler(lw=[1, 2, 3, 4]))

    plot_timecourse(sol,
                    title="plot with a custom cycler",
                    prop_cycle=custom_cycler)
    plt.show()

    st = {'x': 'royalblue'}
    plot_timecourse(sol, styling=st,
                    title="plot, styling {'x': 'royalblue'}")
    plt.show()

    st = {'x': {'c': 'royalblue', 'ls': '--'}}
    plot_timecourse(sol, styling=st,
                    title=f"plot, styling {st}")
    plt.show()

    st = {'x': 9}
    plot_timecourse(sol, styling=st,
                    title="plot, styling {'x': 9}")
    plt.show()

    plot_timecourse(sol,
                    title="plot with a custom palette (Dark2)",
                    palette='Dark2')
    plt.show()

    custom_cycler = (cycler(ls=['-', '--', ':', '-.']) +
                     cycler(lw=[1, 2, 4, 8]))

    plot_timecourse(sol,
                    title="plot with a custom cycler and palette",
                    prop_cycle=custom_cycler, palette='Dark2')
    plt.show()
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
