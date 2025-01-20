import stimator as st
from matplotlib import pyplot as plt
import pandas as pd
from stimator.examples import models

print("---------------- EXAMPLE 1 ------------------")
mdl = """# Example file for S-timator
title Example 1

#reactions (with stoichiometry and rate)
vin  : -> x1     , rate = k1
v2   : x1 ->  x2 , rate = k2 * x1
vout : x2 ->     , rate = k3 * x2

# parameters and initial state
k1 = 1
k2 = 2
k3 = 1
init: (x1=0, x2=0)

# filter what you want to plot
!! x1 x2"""

m = st.read_model(mdl)

with st.style.context("st-bmh"):
    sol1 = m.solve(tf=5.0)
    sol1.plot(xlabel="time", ylabel="conc", palette="Dark2")

    plt.text(
        2.8,
        0.7,
        r" ⟶ x1 ⟶ x2 ⟶ ",
        fontsize="x-large",
        bbox=dict(fc="white", ec="steelblue", lw=2),
    )

    plt.show()

print("---------------- EXAMPLE 2 ------------------")
mtext = """
title a simple 2 enzyme system
v1 : A -> B, rate = Vin*A/(Km + A), V = 0.1, Km = 1
v2 : B -> C, rate = V*B/(Km + B), V = sqrt(4.0), Km = 20

init : A = 1
~ sum = A + B + C
~ sumAB = A + B
-> Vin = 0.1 * step(t, 10)
!! A B C ~
"""

print(mtext)

m1 = st.read_model(mtext)

solution1 = m1.solve(tf=50, title="two enzymes, use !! A C ~")
solution1a = m1.solve(
    tf=50, outputs="A B C sum".split(), title="explicit outputs=[A, B, C, sum]"
)
solution1v = m1.solve(tf=100, outputs=">>", title='outputs=">>"')

print("--- Last time point ----")
print("At t =", solution1.t[-1])
for x in solution1.last:
    print("%-8s= %f" % (x, solution1.last[x]))

print("---------------- EXAMPLE 3 ------------------")
m3 = st.read_model(models.ca.text)
print(models.ca.text)
ms = st.dynamics.ModelSolver(m3, tf=8.0, npoints=2000)
solution3 = ms.solve()

print("---------------- EXAMPLE 4 ------------------")
m4 = st.read_model(models.rossler.text)

print(m4)

solution4 = m4.solve(tf=100.0, npoints=2000, outputs="x1 x2 x3".split())
solution4b = m4.solve(
    tf=100.0, npoints=2000, outputs="~", title='Rossler, outputs="~"'
)


def transformation(vars, t):
    if t > 40.0:
        return (vars[0] - 5.0, vars[1], vars[2])
    else:
        return (-5.0, vars[1], vars[2])


solution4.apply_transf(
    transformation, new_title="Rossler, after a transformation"
)

sols = st.Solutions(
    [solution1, solution1a, solution1v, solution3, solution4b, solution4]
)

with st.style.context("st-seaborn-whitegrid"):
    f, axs = st.plots.prepare_grid(sols, figsize=(9, 6))
    sols.plot(axs=axs, palette="Set1")
    plt.show()

print("---------------- scanning example ------------------")
m3 = st.read_model(models.ca.text)
scans = 0.0, 0.1, 0.3, 0.5, 0.8, 1.0
# scans_k1 = 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9

sols2 = m3.scan({"B": scans}, tf=10.0)

with st.style.context("st-seaborn-whitegrid"):
    f, axs = st.plots.prepare_grid(sols2, figsize=(9, 6))
    sols2.plot(
        what="Ca", axs=axs, legend=False, ylim=(0, 1.5), xlabel="$t$ (min)"
    )
    suptitle = "Cytosolic $Ca^{2+}$ as a function of stimulus strength"
    f.suptitle(suptitle)
    plt.show()

print("---------------- stairway example ------------------")
mtext = """
title a simple 2 enzyme system
v1 : A -> B, rate = Vin*A/(Km + A), V = 0.1, Km = 1
v2 : B -> C, rate = V*B/(Km + B), V = 10, Km = 20
v3 : C ->, rate = kout * C, kout = 1
A = 1

init : B = 0, C = 0

-> Vin = stairway(t, [50, 100, 150, 200, 250], [1, 2, 3, 4, 5])
!! Vin B C
"""

mstair = st.read_model(mtext)

solstairs = mstair.solve(tf=300, title="stairway")

with st.style.context("st-seaborn"):
    f, ax = plt.subplots(figsize=(9, 6))
    solstairs.plot(ax=ax, legend="out")
    plt.show()

print("---- transformation of timecourses to Pandas dataframes -----------")

stairs = pd.DataFrame(sol1.to_dict()).set_index("t")

print(stairs)
