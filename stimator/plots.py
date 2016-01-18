import math
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as pl
import seaborn as sns
sns.set(style='whitegrid')

def _is_string(a):
    return (isinstance(a, str) or
            isinstance(a, unicode))


def _is_sequence(arg):
    return (not hasattr(arg, "strip") and
            hasattr(arg, "__getitem__") or
            hasattr(arg, "__iter__"))

## try:
##     import seaborn as sns
## except ImportError:
##     import smallseaborn as sns

def plot_settings(*args, **kwargs):
    sns.set(*args, **kwargs)


def plotTCs(TCs,
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
            context=None, 
            style=None, 
            palette=None,
            font="sans-serif", 
            font_scale=1.0,
            save2file=None, **kwargs):

    """Generate a graph of the time course using matplotlib and seaborn."""

    # save seaborn data and figure size
    curr_axes_style = sns.axes_style()
    curr_plotting_context = sns.plotting_context()
    curr_color_palette = sns.color_palette()
    original_figsize = tuple(mpl.rcParams['figure.figsize'])
    if context is None:
        context = curr_plotting_context
    sns.set_context(context, font_scale=font_scale)
    if style is None:
        style = curr_axes_style
    sns.set_style(style, rc={"font.family": font})
    if palette is None:
        palette = curr_color_palette
    sns.set_palette(palette)
    mpl.rcParams['lines.markersize']=6
    mpl.rcParams['lines.markeredgewidth']=0.1
    
    if fig_size is not None:
        mpl.rcParams['figure.figsize'] = fig_size
    else:
        mpl.rcParams['figure.figsize'] = (8, 5.5)

    # handle names and titles
    ntc = len(TCs)
    pnames = ['time course %d' % (i+1) for i in range(ntc)]
    for i in range(ntc):
        if titles:
            pnames[i] = titles[i]
        else:
            if TCs[i].title:
                pnames[i] = TCs[i].title

    # find how many plots
    if group:
        nplots = len(group)
    else:
        nplots = ntc

    # compute rows and columns in grid of plots
    ncols = int(math.ceil(math.sqrt(nplots)))
    nrows = int(math.ceil(float(nplots)/ncols))

    # handle axes
    if axis_set is None:
        if figure is None:
            figure = pl.figure()
        axis_set = [figure.add_subplot(nrows, ncols, i+1) for i in range(nplots)]

    # create "plot description" records
    plots_desc = []
    if not group:
        for k, solution in enumerate(TCs):
            rsol = range(len(solution))
            pdesc = dict(title=pnames[k],
                         lines=[(solution.names[i], k, i) for i in rsol])
            plots_desc.append(pdesc)
    else:
        for g in group:
            if _is_string(g):
                pdesc = dict(title=g)
                plines = []
                for k, solution in enumerate(TCs):
                    if g in solution.names:
                        indx = solution.names.index(g)
                        plines.append((pnames[k], k, indx))
                pdesc['lines'] = plines
            else:
                if _is_sequence(g):
                    pdesc = dict(title=' '.join(g))
                    plines = []
                    for vvv in g:
                        for k, solution in enumerate(TCs):
                            if vvv in solution.names:
                                indx = solution.names.index(vvv)
                                if len(TCs) > 1:
                                    plines.append(("%s, %s" % (vvv, pnames[k]),
                                                  k,
                                                  indx))
                                else:
                                    plines.append(("%s" % (vvv), k, indx))
                    pdesc['lines'] = plines
                else:
                    raise StimatorTCError('%s is not a str or seq' % str(g))
            plots_desc.append(pdesc)
    
##         print ('---- plot descriptions |', suptitlegend)
##         for p in plots_desc:
##             print (p)
##         print ('---- end plot descriptions')

    # draw plots
    for i,p in enumerate(plots_desc):
        curraxis = axis_set[i]
        nlines = len(p['lines'])
        use_dots = not TCs[0].dense
        if force_dense:
            use_dots = False
        
        ls, marker = ('None', 'o') if use_dots else ('-', 'None')

        for lname, ltc, li in p['lines']:
            y = TCs[ltc] [li]
            data_loc = np.logical_not(np.isnan(y))
            x = TCs[ltc].t[data_loc]
            y = y[data_loc]
            curraxis.plot(x, y, ls=ls, marker=marker, label=lname)

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
    elif hasattr(TCs, 'title'):
        fig_obj.suptitle(TCs.title)

    if ynormalize and not yrange:
        rs = [a.get_ylim() for a in axis_set]
        common_range = min([l for l,h in rs]), max([h for l,h in rs])
        for a in axis_set:
            a.set_ylim(common_range)

    #pl.tight_layout()

    if save2file is not None:
        figure.savefig(save2file)
    if show:
        if save2file is not None:
            if hasattr(save2file,'read'):
                save2file.close()
        pl.show()

    # restore seaborn styles
    sns.set_context(curr_plotting_context)
    sns.set_style(curr_axes_style)
    sns.set_palette(curr_color_palette)
    mpl.rcParams['figure.figsize'] = original_figsize

#----------------------------------------------------------------------------
#         TESTING CODE
#----------------------------------------------------------------------------

if __name__ == "__main__":
    from modelparser import read_model
    from matplotlib import pyplot as pl
    import timecourse
    plot_settings(style='whitegrid')

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
    demodata_noheader = """
#this is demo data without a header
#t x y z
0       1 0         0
0.1                  0.1

  0.2 skip 0.2 skip this
nothing really usefull here
- 0.3 0.3 this line should be skipped
#0.4 0.4
0.5  - 0.5 - -
0.6 0.6 0.8 0.9

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

    sols = timecourse.Solutions(title='all time courses')
    sols += timecourse.Solution(title='the first tc').load_from_str(demodata)
    sols += timecourse.SolutionTimeCourse().load_from_str(demodata2)
    
    sols.plot(suptitlegend="plotting the two time courses", font_scale=1.5)
    sols.plot(fig_size=(12,6), suptitlegend="with fig_size=(12,6)")  
    sols.plot(group=['z', 'x'], suptitlegend="with group=['z', 'x']")
    sols.plot(group=['z', ('x','y')], suptitlegend="with group=['z', ('x','y')]")
    sols.plot(yrange=(0,2), suptitlegend='with yrange=(0,2)')
    sols.plot(ynormalize=True, suptitlegend='with ynormalize=True')    
    sols.plot(suptitlegend="with force_dense=True", force_dense=True)
    
    f, (ax1, ax2) = pl.subplots(2, 1, sharex=True)
    
    sols.plot(suptitlegend="with given axis_set", force_dense=True,
              axis_set=[ax1, ax2])
    ax1.set_ylabel('variables')
    ax2.set_ylabel('variables')
    ax2.set_xlabel('time')

    sol=timecourse.Solution().load_from_str(demodata)
    sol.plot(group=['z', 'x'], suptitlegend="1 tc with group=['z', 'x']")
    sol.plot(group=['z', ('x','y')], suptitlegend="1tc with group=['z', ('x','y')]")

    sol.load_from('examples/timecourses/TSH2b.txt')
    
    sol.plot(suptitlegend="plotting only one time course")

    f, (ax1, ax2) = pl.subplots(2, 1, sharex=True)
    sol.plot(suptitlegend="plotting on a given axes", axes=ax2)
    ax2.set_ylabel('concentrations')
    ax2.set_xlabel('time')

    print ('\n!! testing transformations ----------------')
    
    sols = timecourse.Solutions(title='all time courses')
    
    sols = timecourse.Solutions(title='all time courses')
    s = timecourse.SolutionTimeCourse(title='original time course').load_from_str(demodata2)
    sols += s
    
    def average(x, t):
        # print ('applying transformation')
        return np.array([t/2.0, (x[0]+x[-1])/2])
    
    s = s.transform(average,
                    newnames=['t/2', 'mid point'], 
                    new_title='after transformation')
    sols += s 
    
    sols.plot(suptitlegend="plotting the two time courses")
    sols.plot(suptitlegend="with force_dense=True", force_dense=True)
    
    tcs = timecourse.readTCs(['TSH2b.txt', 'TSH2a.txt'],
                               'examples/timecourses',
                               names="SDLTSH HTA".split(),
                               verbose=True)
    tcs.plot(suptitlegend="read from file")
    tcs.plot(group=['SDLTSH'], suptitlegend="read from file with group=['SDLTSH']")
    tcs.plot(force_dense=True, suptitlegend="read from file with force_dense=True")

    pl.show()
