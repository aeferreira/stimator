from __future__ import print_function, absolute_import, division
import math
import itertools

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as pl
from cycler import cycler

from stimator.utils import _is_string, _is_sequence


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
    
    if fig_size is not None:
        more_custom_settings['figure.figsize'] = fig_size
    else:
        more_custom_settings['figure.figsize'] = (8, 5.5)
    
    if palette is not None:
        more_custom_settings['axes.prop_cycle'] = cycler('color', list(palette))
    
    more_custom_settings['font.family'] = font
    
    st_list.append(more_custom_settings)
    
    return st_list


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
    
    with pl.style.context(settings):
        
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
                cdesc = []
                ldesc = [(sol.names[i], k, i) for i in range(len(sol))]
                for name, _, _ in ldesc:
                    if name in color_table:
                        c = color_table[name]
                    else:
                        c = next(cyclingcolors)
                        color_table[name] = c
                    cdesc.append(c)
                
                plots_desc.append({'title': pnames[k], 'lines': ldesc, 'colors': cdesc})
        else:
            for g in group:
                ldesc = []
                cdesc = []
                if _is_string(g):
                    pdesc = {'title': g}
                    for k, sol in enumerate(solutions):
                        if g in sol.names:
                            indx = sol.names.index(g)
                            ldesc.append((pnames[k], k, indx))
                else:
                    if not _is_sequence(g):
                        raise StimatorTCError('%s is not a str or seq' % str(g))

                    pdesc = {'title': ' '.join(g)}
                    for vvv in g:
                        for k, sol in enumerate(solutions):
                            if vvv in sol.names:
                                indx = sol.names.index(vvv)
                                if len(solutions) > 1:
                                    ldesc.append(("%s, %s" % (vvv, pnames[k]),
                                                  k,
                                                  indx))
                                else:
                                    ldesc.append(("%s" % (vvv), k, indx))
                for _, k, indx in ldesc:
                    if (k, indx) in color_table:
                        c = color_table[(k, indx)]
                    else:
                        c = next(cyclingcolors)
                        color_table[(k, indx)] = c
                    cdesc.append(c)

                pdesc.update({'lines': ldesc, 'colors': cdesc})
                plots_desc.append(pdesc)

        # draw plots
        for i, p in enumerate(plots_desc):
            curraxis = axis_set[i]
            nlines = len(p['lines'])
            use_dots = not solutions[0].dense
            if force_dense:
                use_dots = False
            
            ls, marker = ('None', 'o') if use_dots else ('-', 'None')

            for ((lname, ltc, li), c) in zip(p['lines'], p['colors']):
                y = solutions[ltc] [li]
                data_loc = np.logical_not(np.isnan(y))
                x = solutions[ltc].t[data_loc]
                y = y[data_loc]
                curraxis.plot(x, y, ls=ls, marker=marker, color=c,
                                    label=lname,
                                    clip_on = False)

            if yrange is not None:
                curraxis.set_ylim(yrange)
            curraxis.set_title(p['title'])
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
    from stimator.modelparser import read_model
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

    sols = Solutions([Solution(title='the first tc').read_str(demodata),
                      Solution().read_str(demodata2)],
                     title='all time courses') 
    
    sols.plot(suptitlegend="plotting the two time courses")
    sols.plot(suptitlegend="with font=serif, palette='rgb'",
              font_scale=1.3, font='serif', palette='rgb')
    sols.plot(suptitlegend="with style=default", style='default')
    sols.plot(suptitlegend="with style=seaborn-darkgrid", style='seaborn-darkgrid')
    sols.plot(suptitlegend="with style=bogus", style='bogus')
    
    sols.plot(fig_size=(12,6), suptitlegend="with fig_size=(12,6)")  
    
    sols.plot(suptitlegend="with force_dense=True", force_dense=True)
    sols.plot(ynormalize=True, suptitlegend='with ynormalize=True')    
    sols.plot(yrange=(0,2), suptitlegend='with yrange=(0,2)')
    
    sols.plot(group=['z', 'x'], suptitlegend="with group=['z', 'x']")
    sols.plot(group=['z', ('x','y')], suptitlegend="with group=['z', ('x','y')]")
    sols.plot(group=['z', ('x','z')], suptitlegend="with group=['z', ('x','z')]")
    
    f, (ax1, ax2) = pl.subplots(2, 1, sharex=True)
    
    sols.plot(suptitlegend="with given axis_set", 
              force_dense=True,
              axis_set=[ax1, ax2])
    ax1.set_ylabel('variables')
    ax2.set_ylabel('variables')
    ax2.set_xlabel('time')

    sol=Solution().read_str(demodata)
    sol.plot(group=['z', 'x'], suptitlegend="1 tc with group=['z', 'x']")
    sol.plot(group=['z', ('x','y')], 
             suptitlegend="1tc with group=['z', ('x','y')]")
    sol.plot(group=['z', ('x','z')], 
             suptitlegend="1tc with group=['z', ('x','z')]")

    sol.read_from('examples/timecourses/TSH2b.txt')
    
    f, (ax1, ax2) = pl.subplots(2, 1, sharex=True)
    
    sol.plot(suptitlegend="plotting on a given axes (1 TC)", axes=ax2)
    ax2.set_ylabel('concentrations')
    ax2.set_xlabel('time')

    print ('\n!! testing transformations ----------------')
       
    sols = Solutions(title='all time courses')
    s = Solution(title='original time course').read_str(demodata2)
    sols += s
    
    def average(x, t):
        # print ('applying transformation')
        return np.array([t/2.0, (x[0]+x[-1])/2.0])
    
    s = s.transform(average,
                    newnames=['t/2', 'mid point'], 
                    new_title='after transformation')
    sols += s 
    
    sols.plot(suptitlegend="plotting original and transf", force_dense=True)
    
    tcs = readTCs(['TSH2b.txt', 'TSH2a.txt'],
                  'examples/timecourses',
                  names="SDLTSH HTA".split(),
                  verbose=False)
    tcs.plot(suptitlegend="read from file")
    tcs.plot(group=['SDLTSH'], suptitlegend="read from file with group=['SDLTSH']")

    pl.show()
