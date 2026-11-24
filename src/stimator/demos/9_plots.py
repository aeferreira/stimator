# ----------------------------------------------------------------------------
#         TESTING CODE
# ----------------------------------------------------------------------------

from cycler import cycler
from matplotlib import pyplot as plt

from stimator import get_examples_path, style
from stimator.plots import one_plot, plot_timecourse, prepare_grid
from stimator.timecourse import Solution, Solutions

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
sol.title = 'TC 1'
sol2.title = 'TC 2'
sols = Solutions([sol, sol2, sol, sol, sol, sol2, sol2])
only_sol = Solutions([sol])

# print(sols)
# print(sol)

example_tc = get_examples_path() / 'TSH2b.txt'

stsh = Solution()
stsh.read_from(example_tc)
stsh.title = 'TSH2b.txt'

# plot_timecourse(sol, title="simple TC plot")
sol.plot(title="simple TC plot")
plt.show()

setts = {"xlabel": 't (s)', "ylabel": 'concentration (mM)',
                "box_aspect": 1}
sol.plot(title=f"TC plot, with settings\n{setts}", **setts)
plt.show()

ttt = "TC plot, modifying axes, using .set() on return value"
sol.plot(title=ttt).set(xlabel='t (s)', ylabel='concentration (mM)')
plt.show()

ttt = "plot with styling={'*': dict(marker='^', ls='none'), tight_t0=False"
plot_timecourse(stsh, styling={'*': {"marker": '^', "ls": 'none'}},
                title=ttt, tight_t0=False)
plt.show()

plot_timecourse(sol,
                title="plot with legend=False",
                tight_t0=False, legend=False)
plt.show()

plot_timecourse(sol,
                title="plot with legend='out'",
                tight_t0=False, legend='out')
plt.show()

plot_timecourse(sol, what='z',
                title="plot with filtering, what='z'")
plt.show()

plot_timecourse(sol, what=[0, 'z'],
                title="plot with filtering, what=[0, 'z']")
plt.show()

with plt.style.context('grayscale'):
    plot_timecourse(sol, title="simple TC plot with grayscale style")
plt.show()

styles = ['st-seaborn-whitegrid', 'seaborn-talk']
with style.context(styles):
    plot_timecourse(sol, title=f"plot with styles \n {styles}",
                    xlabel='time (s)')
plt.show()

custom_cycler = (cycler(color=['c', 'm', 'y', 'k']) +
                    cycler(lw=[1, 2, 3, 4]))

plot_timecourse(sol,
                title="plot with a custom cycler",
                prop_cycle=custom_cycler)
plt.show()

style = {'x': 'chocolate'}
plot_timecourse(sol, styling=style,
                title="plot, styling={'x': 'chocolate'}")
plt.show()

style = {'x_absent': 'chocolate'}
plot_timecourse(sol, styling=style,
                title="plot, styling={'x_absent': 'chocolate'}")
plt.show()

style = {'x': {'c': 'chocolate', 'ls': '--'}}
plot_timecourse(sol, styling=style,
                title=f"plot, styling {style}")
plt.show()

style = {'*': {'ls': '--'}}
plot_timecourse(sol, styling=style,
                title=f"plot, styling {style}")
plt.show()

style = {'x': 9}
plot_timecourse(sol, styling=style,
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

two_tcs = [sol, sol2]
one_plot(two_tcs, title='one_plot of 2 TCs')
plt.show()

f, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))
plot_timecourse(sol,
                title="1st plot", ax=ax1, palette='tab20')
plot_timecourse(sol2,
                title="2nd plot", ax=ax2, palette='tab10')
plt.show()

f, axs = prepare_grid(sols, figsize=(12, 10))
sols.plot(axs=axs, legend=False)
plt.gcf().suptitle('plot Solutions with prepare_grid(), no legend')
plt.show()

f, axs = prepare_grid(sols, figsize=(12, 10))
sols.plot(axs=axs, legend=[False, False, 'out', False, False, False])
plt.gcf().suptitle('plot Solutions with prepare_grid(), legend "out"')
plt.show()

f, axs = prepare_grid(only_sol)
style = {'*': {"marker": 'o'}}
only_sol.plot(axs=axs, styling=style, titles=False)
plt.gcf().suptitle('plot Solutions with only one solution')
plt.show()

grid_kwds = {"figsize": (12, 10), "layout_style": 'compact'}
sols.plot(styling=style, grid_kwds=grid_kwds)
title = 'Plot Solutions without prepare_grid(), "compact" style'
plt.gcf().suptitle(title)
plt.show()
