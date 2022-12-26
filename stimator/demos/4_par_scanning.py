"""S-timator : DEMO of  parameter scanning."""
from stimator import read_model
from stimator.plots import prepare_grid
from matplotlib import pyplot as plt


def run_normal():
    print(__doc__)

    mdl = """
    title Calcium Spikes
    v0         = -> Ca, 1
    v1         = -> Ca, k1*B*step(t, 1.0)
    k1         = 7.3
    B          = 0.4
    export     = Ca ->  , 10 ..
    leak       = CaComp -> Ca, 1 ..

    v2         = Ca -> CaComp, \
                      65 * Ca**2 / (1+Ca**2)
    v3         = CaComp -> Ca, \
                      500*CaComp**2/(CaComp**2+4) * Ca**4/(Ca**4 + 0.6561)
    init       : (Ca = 0.1, CaComp = 0.63655)
    """

    print(mdl)
    m = read_model(mdl)

    title = "CICR model: Effect of stimulus ($\\beta$) on citosolic calcium"

    Bvalues = 0.0, 0.1, 0.3, 0.5, 0.8, 1.0

    titles = [f'$\\beta$ = {b:g}' for b in Bvalues]

    s = m.scan({'B': Bvalues}, tf=10, npoints=1000, titles=titles)

    plt.style.use(['seaborn',
                  {'figure.figsize': (12, 8),
                   'xaxis.labellocation': 'right',
                   'legend.frameon': True,
                   'legend.facecolor': 'white'}])

    f, axs = prepare_grid(s, figsize=(14, 10))
    s.plot(axs=axs, what='Ca', legend=False,
           ylim=[0, 1.5], xlabel='$t$ (min)', box_aspect=1)
    f.suptitle(title, fontsize=16)
    plt.show()

    s.one_plot(what='Ca', title=title, xlabel='$t$ (min)',
               palette='tab20', ylim=(0, 2), xlim=(0,10), label_fmt='$title')
    plt.show()


if __name__ == "__main__":
    run_normal()
