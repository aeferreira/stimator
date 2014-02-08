import numpy
import pylab as pl

def removeMostCrowded(x, knumber, pop_removed=False, distance_fn=None):

    """
    Removes a point with the smaller crowiding distance measured as
    the sum of the distances of the points to their k-nearest neighbours.
    x is a dictionary which values are the lists of objective values for each solution.
    distance_fn
    is an optional function that takes two points and returns the
    distance between them.  If distance_fn is None (the default), the
    Euclidean distance is used.
    This function was written by extensive modification of the Calculate function
    in the kNN module of the Biopython package.
    """
    if len(x) == 0 :
        return x
    if len(x) < 3:
        x.popitem()
        return x
    keys = x.keys()
    keys.sort()
    n_objs = len(x[keys[0]])
    
    # find keys for extremes
    extremes = []
    for i in range(n_objs):
        maximum = x[keys[0]][i]
        tempMax = [keys[0]]
        minimum = x[keys[0]][i]
        tempMin = [keys[0]]
        for j in keys[1:]:
            v = x[j][i]
            if v < minimum:
                tempMin = [j]
                minimum = v
            elif v == minimum:
                tempMin.append(j)                
            elif v > maximum:
                tempMax = [j]
                maximum = v
            elif v == maximum:
                tempMax.append(j)                
        tempMax.extend(tempMin)
        for kk in tempMax:
            if kk not in extremes:
                extremes.append(kk)
    extremes.sort()
    print 'extremes', extremes
    
    #compute distances
    
    distanceMatrix = []
    
    for i in range(len(keys)):
        distances = []
        for j in range(len(keys)):
            if j < i:
                distances.append(distanceMatrix[j][i])
            elif j ==i:
                distances.append(0)
            elif j > i:
                temp = numpy.array(x[keys[i]]) - numpy.array(x[keys[j]])
                d = numpy.sqrt(numpy.dot(temp,temp))
                distances.append(d)
        distanceMatrix.append(distances)
    
    #compute k shortest (position 0 after sorting iss always 0)
    lista =[]
    for i in range(len(keys)):
        distances = distanceMatrix[i]
        distances.sort()
        if knumber < len(distances):
            lista.append(distances[1:knumber+1])
        else:
            lista.append(distances[1:])
    
    #find and remove most crowded
    distancesAndKeys = []
    
    if keys != extremes:
        for i,k in enumerate(keys):
            if k in extremes:
                distancesAndKeys.append((10**300, k))
            else:
                distancesAndKeys.append(((sum(lista[i])), k))
    else:
        for i,k in enumerate(keys):
            distancesAndKeys.append(((sum(lista[i])), k))        
    mck = min(distancesAndKeys)[1]
    if pop_removed:
        mcv = x[mck]
        print 'Most crowded key:', mck, mcv
    del x[mck]
    return x


def pprint(x, showplot=False):
    if len(x) == 0:
        print 'empty dictionary'
    keys = x.keys()
    keys.sort()
    for k in keys:
        print k, x[k]
    if showplot:
        xx = [x[k][0] for k in keys]
        yy = [x[k][1] for k in keys]
        pl.figure()
        pl.plot(xx, yy, 'bo')

def test():
    x = dict(A=(0,1), B=(0,0), C=(1,1), D=(1.2,0.5), E=(1,0), F=(0,0.5),
    G=(0.25, 0.75), H=(0.25,0.5), I=(0.5,1), J=(0.5,0.75), K=(0.5,0.5), L=(0.75,0.5))
    print 'initial'
    pprint(x)
    for n in range(15):
        x = removeMostCrowded(x,3,pop_removed=True)
        print '\nafter removing %d points' % (n+1)
        pprint(x)

if __name__ == '__main__':
    test()
    pl.show()
