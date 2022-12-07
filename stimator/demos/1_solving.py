from stimator import read_model
from matplotlib import pyplot as plt

mdl = """# Example file for S-timator
title Example 1

#reactions (with stoichiometry and rate)
vin  : -> x1     , rate = k1
v2   : x1 ->  x2 , rate = k2 * x1
vout : x2 ->     , rate = k3 * x2

#parameters and initial state
k1 = 1
k2 = 2
k3 = 1
init: (x1=0, x2=0)

# filter what you want to plot
!! x1 x2"""

m = read_model(mdl)

print('========= model ========================================')
print(mdl)
print('--------------------------------------------------------')

ax = plt.gca()
print(plt.get_fignums())
m.solve(tf=5.0).plot(fig_size=(6.5, 6))
print(plt.get_fignums())
# plt.show()

print('========= check ========================================')
print('between plots')
print('--------------------------------------------------------')

m.solve(tf=20.0).plot(fig_size=(6.5, 6))
print(plt.get_fignums())
plt.show()
print(plt.get_fignums())
