from stimator import read_model, get_examples_path
from matplotlib import pyplot as plt

# --- example: Glyoxalase system in L. Infantum --------------------

m1 = read_model("""
title example 2: Glyoxalase system in L. Infantum

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

print('=================================================')
print('Parameter estimation: glyoxalase system example')
print(m1)
print('-------- an example with two time courses --------------')

tcdir = get_examples_path()

optimum = m1.estimate(tc_dir=tcdir, names=['SDLTSH', 'HTA'],
                      dump_generations=True)
# convergence_noimprovement=40)
# ... intvarsorder=(0,2,1) ...

print(optimum)

plt.style.use('bmh')

_, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey='row')
optimum.plot(0, ax=ax1, xlabel='t (s)', ylabel='conc (microM)')
optimum.plot(1, ax=ax2, xlabel='t (s)', ylabel='conc (microM)')
plt.show()

optimum.print_generations()

# --- example 2, fitting a transformation --------------------

print('=================================================')
print('Parameter estimation: glyoxalase system example')
print('-------- an example with a transformation --------------')

mtransf = read_model("""
title example 2, fitting a transformation

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

~sdlx2 = 2 * SDLTSH # the transformation to fit

!! sdlx2

init : (SDLTSH = 7.69231E-05, HTA = 0.1357)

""")

optimum = mtransf.estimate(tc_dir=tcdir,
                           timecourses=['tc_double.txt'],
                           names=['sdlx2', 'SDLTSH', 'HTA'])

print(optimum)
optimum.plot()
plt.show()

# --- Same example with unknown initial values --------------------

print('-------- same example with unknown initial values --------------')

m2 = m1.copy()

# Now, assume init.HTA is uncertain
m2.set_bounds('init.HTA', (0.05, 0.25))
# do not estimate Km1 and Km2, to help the analysis
m2.reset_bounds('Km1')
m2.reset_bounds('Km2')

# VERY IMPORTANT:
# only one time course can be used:
# cannot fit one initial value using several timecourses!!!
# change the default pop_size(80) to 60

optimum = m2.estimate(timecourses=['TSH2a.txt'],
                      tc_dir=tcdir,
                      opt_settings={'pop_size': 60},
                      names=['SDLTSH', 'HTA'])

print(optimum)
optimum.plot()
plt.show()
