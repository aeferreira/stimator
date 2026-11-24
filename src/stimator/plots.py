import math
from itertools import cycle
from string import Template

import matplotlib as mpl
import numpy as np
from cycler import Cycler, cycler
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from packaging.version import Version

from . import style
from .utils import _is_sequence, _is_string


def _find_what_to_plot(timecourse, what):
    if what is None:
        return timecourse.names[:]
    if _is_string(what):
        what = [what]
    if not _is_sequence(what):
        raise ValueError(f"No data for {what:s} in timecourse")
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
                raise ValueError(f"No data for {elem:s} in timecourse")
            if elem not in cols2plot:
                cols2plot.append(elem)
    return cols2plot


_MPL_QUALIT = ['Pastel1', 'Pastel2', 'Paired', 'Accent',
               'Dark2', 'Set1', 'Set2', 'Set3', 'tab10',
               'tab20', 'tab20b', 'tab20c']


def _get_colormap_colortable():
    mpl_ver = Version(mpl.__version__)
    if mpl_ver >= Version('3.5'):
        return {cm: mpl.colormaps[cm].colors for cm in _MPL_QUALIT}
    cm_table = {}
    for cm in _MPL_QUALIT:
        cm_obj = mpl.cm.get_cmap(name=cm, lut=None)
        colors = [cm_obj(i) for i in range(cm_obj.N)]
        cm_table[cm] = colors
    return cm_table


MPL_QUALITATIVE_TABLE = _get_colormap_colortable()


def get_color_list(pname):
    """get a list of RGB values from a palette name.

    mpl qualitative names are accepted.
    TODO: accept seaborn and palettable generated palettes as pname"""
    if pname in _MPL_QUALIT:
        return MPL_QUALITATIVE_TABLE[pname]
    raise ValueError(f'palette name {pname} not found')


def _get_new_prop_cycle(prop_cycle, palette):
    if prop_cycle is None and palette is None:
        return mpl.rcParams['axes.prop_cycle']
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


def generate_line_handles(names, prop_cycle, palette, styling):
    """build table of Line2D handles with no data.

    Cycle through a prop_cycle obtained
    considering 'palette' and 'prop_cycle' arguments
    Overide styles using the styles in 'styling'
    """

    line_handles = {}

    used_prop_cycle = _get_new_prop_cycle(prop_cycle, palette)

    for name, stl in zip(names, cycle(used_prop_cycle)):
        line_handles[name] = Line2D([], [])
        line_handles[name].set(**stl)

    if styling is not None:

        if not isinstance(styling, dict):
            raise TypeError("'styling' parameter must be a dict")
        props = list(used_prop_cycle)

        for k, v in styling.items():
            # get matching handles
            if k.strip() == '*':
                set_handles = line_handles.values()
            else:
                set_handles = [line_handles[n] for n in line_handles if n == k]
                # TODO: allow regular expressions
            if isinstance(v, dict):
                for handle in set_handles:
                    handle.set(**v)
            elif isinstance(v, int):
                for handle in set_handles:
                    handle.set(**props[v])
            elif _is_string(v):
                for handle in set_handles:
                    handle.set(color=v)
            else:
                raise ValueError(f"'{v}' is an invalid color or style ")
    return line_handles


def draw_legend(ax, legend_argument):
    if legend_argument is False:
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()
        return
    if legend_argument is True:
        legend_argument = 'best'

    handles, lbls = ax.get_legend_handles_labels()

    if not legend_argument.startswith('out'):
        ax.legend(handles, lbls, loc=legend_argument)
    else:
        # must use tight layout, in this case
        f = ax.get_figure()
        if not (f.get_constrained_layout() or f.get_tight_layout()):
            f.set_tight_layout(True)
        ax.legend(handles, lbls, loc='upper left',
                  bbox_to_anchor=(1.01, 1),
                  borderaxespad=0)


def plot_timecourse(timecourse, what=None, ax=None,
                    prop_cycle=None, palette=None, styling=None,
                    legend='best', title=None, tight_t0=True,
                    **axes_settings):
    if ax is None:
        ax = plt.gca()

    # find what (which columns) to plot
    cols_to_plot = _find_what_to_plot(timecourse, what)

    # build table of Line2D handles with no data

    line_handles = generate_line_handles(cols_to_plot,
                                         prop_cycle, palette, styling)

    # plot the lines
    for name in cols_to_plot:
        x = timecourse.t
        y = timecourse[name]
        line_handles[name].set_data(x, y)
        line_handles[name].set(label=name)
        ax.add_line(line_handles[name])
    ax.autoscale(enable=None)

    # apply axes settings
    ax.set(**axes_settings)

    # draw legend
    draw_legend(ax, legend)

    # draw plot title
    if title is not None:
        if title is not False:
            ax.set_title(title)
    elif timecourse.title:
        ax.set_title(timecourse.title)

    # remove the left margin for the x-axis
    # to start the plot at t0, if 'tight_t0'
    if tight_t0 and (timecourse.t[0] < timecourse.t[-1]):
        ax.set_xlim(timecourse.t[0], None)
    return ax


def plot_estim_optimum_timecourse(opt, tc_index=0, ax=None, title=None,
                                  exp_style=None, tight_t0=True, legend='best',
                                  opt_label_fmt='Predicted ${name}',
                                  exp_label_fmt='$name',
                                  **axes_settings):
    if ax is None:
        ax = plt.gca()

    tcstats = opt.tcdata[tc_index]
    expsol = opt.optimizer.tc[tc_index]
    symsol = opt.optimum_dense_tcs[tc_index]

    if exp_style is None:
        exp_style = {"marker": 'o', "linestyle": 'none'}

    # retrieve 'prop_cycle' and 'palette' and 'styling' arguments
    prop_cycle = axes_settings.pop('prop_cycle', None)
    palette = axes_settings.pop('palette', None)
    styling = axes_settings.pop('styling', None)

    # ignore 'what' (for now...)
    _ = axes_settings.pop('what', None)
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

    # build tables of Line2D handles with no data
    line_handles = generate_line_handles(what2plot,
                                         prop_cycle, palette, styling)

    exp_line_handles = generate_line_handles(what2plot,
                                             prop_cycle, palette, styling)
    # overide the settings for experimental data
    for name in exp_line_handles:
        exp_line_handles[name].set(**exp_style)

    # defaults for line labels
    if exp_label_fmt is None:
        exp_label_fmt = '$name'
    if opt_label_fmt is None:
        opt_label_fmt = '${name}_pred'

    # plot experimental data
    for name in what2plot:
        x = expsol.t
        y = expsol[name]
        exp_line_handles[name].set_data(x, y)
        mapping = {'plot_indx': 0, 'name': name}
        explabel = Template(exp_label_fmt).safe_substitute(mapping)
        exp_line_handles[name].set(label=explabel)
        ax.add_line(exp_line_handles[name])

    # plot predicted data
    for name in what2plot:
        x = symsol.t
        y = symsol[name]
        line_handles[name].set_data(x, y)
        mapping = {'plot_indx': 1, 'name': name}
        opt_name = Template(opt_label_fmt).safe_substitute(mapping)
        line_handles[name].set(label=opt_name)
        ax.add_line(line_handles[name])

    ax.autoscale(enable=None)

    # apply axes settings
    ax.set(**axes_settings)

    # draw legend
    draw_legend(ax, legend)

    # draw plot title
    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title("%s (%d pt) %g" % tcstats)

    # remove the left margin for the x-axis
    # to start the plot at t0
    if tight_t0 and (symsol.t[0] < symsol.t[-1]):
        ax.set_xlim(symsol.t[0], None)

    return ax


def one_plot(time_courses, what=None, ax=None, label_fmt='${title}: ${name}',
             prop_cycle=None, palette=None, styling=None,
             legend='best', title=None, tight_t0=True,
             **axes_settings):

    if ax is None:
        ax = plt.gca()

    # find what (which columns) to plot
    cols_to_plot = []
    for tc in time_courses:
        cols_to_plot.append(_find_what_to_plot(tc, what))

    # generate names for the dict of handles with tc ref
    all_names = []
    for i, cols in enumerate(cols_to_plot):
        icols = [f'{i}:::{c}' for c in cols]
        all_names.extend(icols)

    # build table of Line2D handles with no data

    line_handles = generate_line_handles(all_names,
                                         prop_cycle, palette, styling)

    # plot the lines
    for fullname in all_names:
        # find tc index and name
        tc_index, name = fullname.split(':::', 1)
        tc_index = int(tc_index)
        timecourse = time_courses[tc_index]
        x = timecourse.t
        y = timecourse[name]
        line_handles[fullname].set_data(x, y)
        mapping = {'plot_indx': tc_index,
                   'title': timecourse.title,
                   'name': name}
        label = Template(label_fmt).safe_substitute(mapping)
        line_handles[fullname].set(label=label)
        ax.add_line(line_handles[fullname])
    ax.autoscale(enable=None)

    # apply axes settings
    ax.set(**axes_settings)

    # draw legend
    draw_legend(ax, legend)

    # draw plot title
    if title is not None:
        ax.set_title(title)

    # remove the left margin for the x-axis
    # to start the plot at t0, if 'tight_t0'
    if tight_t0 and (time_courses[0].t[0] < time_courses[0].t[-1]):
        ax.set_xlim(time_courses[0].t[0], None)
    return ax


def prepare_grid(sols, layout_style=None, **kwargs):

    is_more35 = Version(mpl.__version__) >= Version('3.5')
    subplots_kwargs = kwargs.copy()
    has_layoutkwds = 'tight_layout' in subplots_kwargs
    has_layoutkwds = has_layoutkwds or ('constrained_layout' in subplots_kwargs)
    has_layoutkwds = has_layoutkwds or ('layout' in subplots_kwargs)
    if not has_layoutkwds:
        if is_more35:
            subplots_kwargs['layout'] = 'tight'
        else:
            subplots_kwargs['tight_layout'] = True

    # find how many plots
    nplts = len(sols)

    # compute shape of grid
    if not ('nrows' in kwargs or 'ncols' in kwargs):
        ncols = math.ceil(math.sqrt(nplts))
        nrows = math.ceil(float(nplts)/ncols)
        subplots_kwargs['nrows'] = nrows
        subplots_kwargs['ncols'] = ncols

    if layout_style == 'compact':
        subplots_kwargs['sharey'] = 'all'
        subplots_kwargs['sharex'] = 'all'
        if is_more35:
            subplots_kwargs['layout'] = None
        else:
            subplots_kwargs['constrained_layout'] = False
            subplots_kwargs['tight_layout'] = False

    f, axs = plt.subplots(**subplots_kwargs)

    # hide empty axes
    if not isinstance(axs, np.ndarray):
        axs = np.array(axs)
    for i, ax in enumerate(axs.flat):
        if i >= len(sols):
            ax.set_visible(False)
    if layout_style == 'compact':
        f.subplots_adjust(hspace=0, wspace=0)
    return f, axs


def plotTCs(sols, axs=None, grid_kwds=None, titles=True, **kwargs):
    if axs is None:
        if grid_kwds is None:
            grid_kwds = {}
        f, axs = prepare_grid(sols, **grid_kwds)
    else:
        f = plt.gcf()

    # get legend argument
    # might be a sequence to cycle through
    legends = kwargs.pop('legend', 'best')
    if not _is_sequence(legends):
        legends = [legends]

    is_compact = False
    if grid_kwds is not None:
        layout = grid_kwds.get('layout_style', '')
        is_compact = (layout == 'compact')

    if not isinstance(axs, np.ndarray):
        usable_axs = np.array(axs)
    else:
        usable_axs = axs
    if titles is False or is_compact:
        kwargs['title'] = False
    for tc, ax, leg in zip(sols, usable_axs.flat, cycle(legends)):
        plot_timecourse(tc, ax=ax, legend=leg, **kwargs)
    return f, axs

