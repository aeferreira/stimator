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


def get_color_list(pname):
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
        color_list = get_color_list(palette)
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

    # keep a table of handles of lines
    line_handles = {}

    # plot (as lines)
    for name in cols_to_plot:
        x = timecourse.t
        y = timecourse[name]
        line, *_ = ax.plot(x, y, label=name, **kwargs)
        line_handles[name] = line

    # overide prop_cycler styles if these are specified on the styling dict
    # for some variables
    if styling is not None:
        stylingdict = _get_overiding_prop_table(styling, used_prop_cycle)
        for name in stylingdict:
            line_handles[name].set(label=name, **(stylingdict[name]))

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


def plot_estim_optimum_timecourse(opt, tc_index=0, ax=None, exp_style=None,
                                  opt_suffix='_pred', opt_prefix='',
                                  **kwargs):
    if ax is None:
        ax = plt.gca()

    tcstats = opt.tcdata[tc_index]
    expsol = opt.optimizer.tc[tc_index]
    symsol = opt.optimum_dense_tcs[tc_index]

    if exp_style is None:
        exp_style = dict(marker='o', linestyle='none')

    title = kwargs.pop('title', None)
    tight_t0 = kwargs.pop('tight_t0', True)

    # locally set the prop_cycler if prop_cycle or palette are given
    prop_cycle = kwargs.pop('prop_cycle', None)
    palette = kwargs.pop('palette', None)
    used_prop_cycle = _get_new_prop_cycle(prop_cycle, palette)
    if used_prop_cycle is not None:
        ax.set_prop_cycle(used_prop_cycle)
    else:
        used_prop_cycle = mpl.rcParams['axes.prop_cycle']

    # retrieve the 'styling' argument
    styling = kwargs.pop('styling', None)

    # ignore 'what' (for now...)
    _ = kwargs.pop('what', None)

    # keep a table of handles of lines
    line_handles = {}

    # find what to plot
    what2plot = []
    for line in range(len(expsol)):
        # count NaN and do not plot if they are most of the time series
        yexp = expsol[line]
        nnan = len(yexp[np.isnan(yexp)])
        if nnan >= expsol.ntimes-1:
            continue
        # otherwise keep the name
        what2plot.append(expsol.names[line])

    # plot experimental data
    ax.set_prop_cycle(used_prop_cycle)
    for name in what2plot:
        yexp = expsol[name]
        newstyles = dict(kwargs)
        newstyles.update(exp_style)
        line, *_ = ax.plot(expsol.t, yexp, label=name, **newstyles)
        line_handles[name] = line

    # plot predicted data
    ax.set_prop_cycle(used_prop_cycle)  # reuse props
    for name in what2plot:
        ysim = symsol[name]
        opt_name = opt_prefix + name + opt_suffix
        line, *_ = ax.plot(symsol.t, ysim, label=opt_name, **kwargs)
        line_handles[opt_name] = line

    # overide prop_cycler styles if these are specified on the styling dict
    # for some variables
    if styling is not None:
        stylingdict = _get_overiding_prop_table(styling, used_prop_cycle)
        for name in stylingdict:
            line_handles[name].set(label=name, **(stylingdict[name]))

    # draw legend
    handles, lbls = ax.get_legend_handles_labels()
    ax.legend(handles, lbls, loc='best')

    # draw plot title
    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title("%s (%d pt) %g" % tcstats)

    # remove the left margin for the x-axis
    # to start the plot at t0
    if tight_t0:
        if symsol.t[0] < symsol.t[-1]:
            ax.set_xlim(symsol.t[0], None)
    return ax


# ----------------------------------------------------------------------------
#         TESTING CODE
# ----------------------------------------------------------------------------


if __name__ == "__main__":
    from stimator import get_examples_path
    from stimator.timecourse import Solution

    demodata = """
#this is demo data with a header
t x y z
0       0.95 0         0
0.1   0.1  0.5        0.2

  0.2 0.2 0.2 skip this
nothing really usefull here
- 0.3 0.3 this line should be skipped
#0.4 0.4
0.3 0.4 0.5 0.55
0.4 0.5 0.6 0.7
0.5 0.6 0.8 0.9
0.55 0.7 0.85 0.95
0.6  0.65 0.5 - -

"""

    demodata2 = """
#this is demo data with a header
t x2 y2 z2
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

    example_tc = get_examples_path() / 'TSH2b.txt'

    stsh = Solution()
    stsh.read_from(example_tc)

    # plot_timecourse(sol, title="simple TC plot")
    sol.plot(title="simple TC plot")
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

    plot_timecourse(sol,
                    title="1st plot", palette='Dark2')
    plot_timecourse(sol2,
                    title="2nd plot, same axes as 1st plot", palette='Dark2')
    plt.show()

    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))
    plot_timecourse(sol,
                    title="1st plot", ax=ax1, palette='Pastel2')
    plot_timecourse(sol2,
                    title="2nd plot", ax=ax2, palette='tab10')
    plt.show()
