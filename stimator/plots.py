from __future__ import print_function, absolute_import, division
import math
import itertools

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as pl
from cycler import cycler

from stimator.utils import _is_string, _is_sequence, _is_integer

def _get_plot_grid_dict(g):
    if g is None:
        return {'nrows':1, 'ncols':1, 'names':[0]}
    if _is_integer(g) and g > 0:
        return {'nrows':g, 'ncols':1, 'names':[i for i in range(g)]}
    if _is_string(g):
        return {'nrows':1, 'ncols':len(g), 'names': list(g)}
    if _is_sequence(g):
        if len(g) == 2 and _is_integer(g[0]) and _is_integer(g[1]):
            if g[0] > 0 and g[1] > 0:
                nrows = g[0]
                ncols = g[1]
                names =[]
                for i in range(nrows):
                    for j in range(ncols):
                        names.append((i, j))
                return {'nrows':nrows, 'ncols':ncols, 'names':names}
        nrows = len(g)
        names = []
        for i in g:
            if not (_is_string(i) or _is_sequence(i)):
                raise ValueError('wrong grid specification')
        ncols = min([len(i) for i in g])
        for i in g:
            names.extend(list(i[:ncols]))
        return {'nrows': nrows, 'ncols':ncols, 'names':names}
    raise TypeError('wrong grid specification')


class PlotGrid(object):
    def __init__(self,
            solutions,
            grid=None,
            #axis_set=None,
            fig_size=None,
            ynormalize=False,
            yrange=None,
            suptitlegend=None,
            style=None, 
            palette=None,
            font="sans-serif",
            **kwargs):
        self.solutions = solutions
        self.fig_size = fig_size
        self.ynormalize = ynormalize
        self.yrange = yrange
        self.suptitlegend = suptitlegend
        self.style = style
        self.pallete = palette
        self.font = font
        #self.axis_set = axis_set
                
        grid_desc = _get_plot_grid_dict(grid)
        nrows = grid_desc['nrows']
        ncols = grid_desc['ncols']

        # create figure and axes
        fig, axes = pl.subplots(nrows, ncols, figsize=fig_size)
        if not hasattr(axes, 'shape'): # only one axes was created
            print(axes)
            self.axes_table = {grid_desc['names'][0]: axes}
        else:
            print(axes)
            print(axes.shape)
        
            self.axes_table = {}
            for a, name in zip(axes.ravel(), grid_desc['names']):
                self.axes_table[name] = a
        print(self.axes_table)
        
        self.fig = fig
        self.axes = axes
        print(self.fig.axes)

    def _get_axes_from_loc(self, loc):
        if loc is None:
            # the one and only or default to the first
            # fig.axes is a list
            return self.fig.axes[0]
        if loc in self.axes_table:
            return self.axes_table[loc]
        try:
            ax = self.axes[loc]
            if hasattr(ax, 'plot'):
                return ax
        except (IndexError, ValueError):
            pass
        raise ValueError('panel {} does not exist'.format(loc))
            
    
    def plot123(self, loc=None, show=False):
        ax = self._get_axes_from_loc(loc)
        ax.plot(range(3), range(1,4), 'or')
        if show:
            pl.show()
        return self
    
    def savefig(self, filename):
        self.fig.savefig(filename)
    
    def show(self):
        pl.show()

    def plot_solution(self, solution_index,
                loc=None, show=False,
                title=None,
                fig_size=None,
                yrange=None,
                legend=True,
                force_dense=False,
                style=None, 
                palette=None,
                font="sans-serif", **kwargs):

        """Generate a graph of the time course using matplotlib.
           
           Called by .plot() member function of class timecourse.Solutions"""

        ax = self._get_axes_from_loc(loc)
        sol = self.solutions[solution_index]
        
        pname = 'timecourse %d' % solution_index
        if title:
            pname = title
        else:
            if sol.title:
                pname = sol.title
        
        settings = _prepare_settigs(style, palette, font, fig_size)
        
        with pl.style.context(settings):
            
            cyl = [c['color'] for c in mpl.rcParams['axes.prop_cycle']]
            cyclingcolors = itertools.cycle(cyl)
            
            # create "plot description" record
            
            color_table = {}
            lines_desc = []
            for i in range(len(sol)):
                name = sol.names[i]
                line = {'name': name,
                        'solution_index': solution_index,
                        'var_index': i}
                if name in color_table:
                    c = color_table[name]
                else:
                    c = next(cyclingcolors)
                    color_table[name] = c
                line['color'] = c
                lines_desc.append(line)
            
            plot_desc = {'title': pname, 'lines': lines_desc}

            # draw plot
            use_dots = not self.solutions[0].dense
            if force_dense:
                use_dots = False
            
            ls, marker = ('None', 'o') if use_dots else ('-', 'None')
            
            _plotTC(plot_desc['lines'], self.solutions, plot_desc['title'], ls, marker, ax)

            if yrange is not None:
                ax.set_ylim(yrange)
            if legend:
                h, l = ax.get_legend_handles_labels()
                ax.legend(h, l, loc='best')
            ax.set_xlabel('')
            ax.set_ylabel('')
            
        if show:
            pl.show()
        return self

    def yscale_identical(self):
        rs = [a.get_ylim() for a in self.fig.axes]
        common_range = min([l for l,h in rs]), max([h for l,h in rs])
        for a in self.fig.axes:
            a.set_ylim(common_range)
        return self

def _prepare_settigs(style, palette, font, fig_size):
    st_list = []
    
    more_custom_settings = {#'lines.markersize': 5,
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
    
##     if fig_size is not None:
##         more_custom_settings['figure.figsize'] = fig_size
##     else:
##         more_custom_settings['figure.figsize'] = (8, 5.5)
    
    if palette is not None:
        more_custom_settings['axes.prop_cycle'] = cycler('color', list(palette))
    
    more_custom_settings['font.family'] = font
    
    st_list.append(more_custom_settings)
    
    return st_list
       

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
    
    settings = _prepare_settigs(style, palette, font, fig_size)
    
    with pl.style.context(settings):
        
        # handle axes
        if axis_set is None:
            if figure is None:
                figure = pl.figure()
            axis_set = [figure.add_subplot(nrows, ncols, i+1) for i in range(nplts)]

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
                            line={'name': pnames[k],
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
                h, l = curraxis.get_legend_handles_labels()
                curraxis.legend(h, l, loc='best')
            curraxis.set_xlabel('')
            curraxis.set_ylabel('')
        
        # draw suptitle (needs a figure object)
        fig_obj = pl.gcf()
        if suptitlegend is not None:
            fig_obj.suptitle(suptitlegend)
        elif hasattr(solutions, 'title'):
            fig_obj.suptitle(solutions.title)

        if ynormalize and not yrange:
            rs = [a.get_ylim() for a in axis_set]
            common_range = min([l for l,h in rs]), max([h for l,h in rs])
            for a in axis_set:
                a.set_ylim(common_range)

        # pl.tight_layout()

        if save2file is not None:
            figure.savefig(save2file)
        if show:
            if save2file is not None:
                if hasattr(save2file,'read'):
                    save2file.close()
            pl.show()


def plot_estim_optimum(opt, figure=None, 
                       axis_set=None,
                       fig_size=None,
                       style=None, 
                       palette=None,
                       font="sans-serif", 
                       save2file=None,
                       show=False):

    settings = _prepare_settigs(style, palette, font, fig_size)
    
    with pl.style.context(settings):
        if axis_set is None:
            if figure is None:
                figure = pl.figure()

        bestsols = opt.optimum_dense_tcs
        expsols = opt.optimizer.tc
        tcstats = opt.tcdata
        nplts = len(bestsols)
        ncols = int(math.ceil(math.sqrt(nplts)))
        nrows = int(math.ceil(float(nplts)/ncols))
        if axis_set is None:
            axis_set = [figure.add_subplot(nrows, ncols,i+1) for i in range(nplts)]
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
                             marker=mexp, ls=lsexp, color=color, clip_on = False)
                subplot.plot(symsol.t, ysim, 
                             marker=msim, ls=lssim, color= color,
                             label='%s' % xname, clip_on = False)
            subplot.legend(loc='best')

        if save2file is not None:
            figure.savefig(save2file)
        if show:
            if save2file is not None:
                if hasattr(save2file, 'read'):
                    save2file.close()
            pl.show()


def plot_generations(opt, generations=None,
                     pars=None, figure=None, show=False,
                     fig_size=None,
                     style=None, palette=None, font="sans-serif"):
    
    if not opt.generations_exist:
        raise IOError('file generations.txt was not generated')
    
    settings = _prepare_settigs(style, palette, font, fig_size)
    settings.append({'lines.markeredgewidth': 1.0})
    
    with pl.style.context(settings):

        if figure is None:
            figure = pl.figure()
        figure.clear()

        if generations is None:
            all_gens = list(range(opt.optimization_generations +1))
            dump_generations = all_gens

        n_gens = len(dump_generations)
        
        if pars is None:
            first2 = opt.parameters[:2]
            pars = [p[0] for p in first2]
        
        pnames = [p[0] for p in opt.parameters]
        
        colp0 = pnames.index(pars[0])
        colp1 = pnames.index(pars[1])
        
        scores_col = len(opt.parameters)
        
        #ax1 = pl.subplot(1,2,1)
        #ax2 = pl.subplot(1,2,2)
        ax2 = pl.subplot(1,1,1)
        
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
                    #ax1.plot(solx, soly, marker='o', ls='None', label=gen)
                    ## for px, py in zip(objx, objy):
                        ## ax2.axhline(py, xmin=px, xmax=px*0.01, color='black')
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
            pl.show()

# ----------------------------------------------------------------------------
#         TESTING CODE
# ----------------------------------------------------------------------------

if __name__ == "__main__":
    from stimator import read_model
    from stimator.timecourse import Solution, Solutions, readTCs

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

    t = np.linspace(0, 6, 11)
    y = np.sin(t)
    y2 = np.cos(t - 0.5)
    header = 't s1 s2'
    lines = ['{} {} {}'.format(i, j, k) for (i,j,k) in zip(t, y, y2)]
    demosin = '\n'.join([header]+lines)
    
    sine_sols = Solutions([Solution(title='sin and cos').read_str(demosin)])
    
    sols = Solutions([Solution(title='the first tc').read_str(demodata),
                      Solution().read_str(demodata2)],
                     title='all time courses')
    pg = PlotGrid(sols, grid=(2,2))
    pg.plot123((1,1))
    pg.plot123((0,0))
    pl.show()
    
    pg = PlotGrid(sols)
    pg.plot123()
    pl.show()
    
    pg = PlotGrid(sols, grid=('ABCD'))
    pg.plot123('B')
    pg.plot123(0)
    pl.show()

    pg = PlotGrid(sols, grid=('AB','CD'))
    pg.plot123('B')
    pg.plot123((1,0))
    #pg.plot123(0)
    pl.show()

    pg = PlotGrid(sine_sols)
    pg.plot_solution(0)
    pl.show()
    
    sols.plot(suptitlegend="plotting the two time courses")
    sols.plot(suptitlegend="with font=serif, palette='rgb'",
              font_scale=1.3, font='serif', palette='rgb')
    pl.show()
