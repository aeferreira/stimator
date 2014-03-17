
from matplotlib import pylab as pl

dump_generations = [0, 1, 3, 5, 7, 8, 10, 13, 14]

if __name__ == "__main__":
    if dump_generations is not None:
        lcolors = ['white'] + ['0.8', '0.4']
        lcolors = lcolors + ['yellow', 'cyan', 'red', 'green', 'blue', 'black']
        lcolors = lcolors * 10
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
                    c = lcolors[cindex]
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
                print gen
            elif reading:
                line = [float(x) for x in line.split()]
                solx.append(line[0])
                soly.append(line[1])
                objx.append(-line[2])
                objy.append(-line[3])
            else:
                continue
        f.close()
        pl.legend(loc=0)
        ax1.grid()
        ax2.grid()
        pl.show()
