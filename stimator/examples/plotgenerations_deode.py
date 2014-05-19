
from matplotlib import pylab as pl
import matplotlib.cm as cm

dump_generations = list(range(0,201,10))
n_gens = len(dump_generations)

if __name__ == "__main__":
    if dump_generations is not None:
        ax1 = pl.subplot(1,2,1)
        ax2 = pl.subplot(1,2,2)
        # parse generations
        gen = -1
        f = open('generations.txt')
        solx = []
        soly = []
        objx = []
        objy = []
        reading = False
        cindex = -1
        for line in f:
            line = line.strip()
            if line == '' and reading:
                if len(solx) > 0:
                    cindex += 1
                    c = cm.jet(cindex/float(n_gens), 1)
                    ax1.plot(solx, soly, color=c, marker='o', ls='None', label = gen)
                    ax2.plot(objx, objy, color=c, marker='o', ls='None', label = gen)
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
                    print 'generation', gen
            elif reading:
                line = [float(x) for x in line.split()]
                solx.append(line[0])
                soly.append(line[2])
                objx.append(igen)
                objy.append(line[4])
            else:
                continue
        f.close()
        pl.legend(loc=0)
        ax1.grid()
        ax1.set_title('population')
        ax1.set_xlabel('V1')
        ax1.set_ylabel('Km1')
        ax2.grid()
        ax2.set_title('scores')
        ax2.set_ylim(0, 0.01)
        ax2.set_xlabel('generation')
        pl.show()
