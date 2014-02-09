import numpy
import pylab as pl

def removeMostCrowded(x, knumber, remove_n = 1, print_removed = False):

    """
    Removes a point with the smaller crowiding distance measured as
    the sum of the distances of the points to their k-nearest neighbours.
    x is a dictionary which values are the lists of objective values for each solution.
    """
    
    n_points = len(x)
    if n_points == 0 :
        return x
    if n_points < 3:
        x.popitem()
        return x
    keys = x.keys()
    keys.sort() # TODO: not necessary??
    n_objs = len(x[keys[0]])
    
    points = range(n_points) #holds current indexes of points still present
    
    #compute matrix of objectives
    obj_matrix = numpy.empty((n_points, n_objs))
    for i,k in enumerate(keys):
        obj_matrix[i, :] = x[k]

    # find indexes of extreme points (max and min in each dimension)
    
    maxima = numpy.amax(obj_matrix, axis=0)
    minima = numpy.amin(obj_matrix, axis=0)
    
    extremes = []
    for i in range(n_objs):
        dimextremes = []
        for j in points:
            v = obj_matrix[j][i]
            if v == minima[i] or v == maxima[i]:
                dimextremes.append(j)                
        for j in dimextremes:
            if j not in extremes:
                extremes.append(j)
    extremes.sort()
    print 'extremes', [keys[i] for i in extremes]

    #compute distances
    # TODO: use numpy function
    distanceMatrix = []
    
    for i in range(n_points):
        distances = []
        for j in range(n_points):
            if j < i:
                distances.append(distanceMatrix[j][i])
            elif j ==i:
                distances.append(0.0)
            elif j > i:
                temp = numpy.array(obj_matrix[i]) - numpy.array(obj_matrix[j])
                d = numpy.sqrt(numpy.dot(temp,temp))
                distances.append(d)
        distanceMatrix.append(distances)
    
    #sort distances for each point
    distances_array = []
    for i in points:
        dd = zip(list(distanceMatrix[i]), range(n_points))
        dd.sort()
        distances_array.append(dd)
    
    #compute k shortest (note: position 0 after sorting is always 0.0)
    #TODO: repeat from here
    kdistances =[]
    ksums = []
    
    for i in points:
        kdistances.append(distances_array[i][1:knumber+1])
        ksums.append(sum([ d[0] for d in kdistances[-1]]))
    

    print '\nk-distances'
    for k,d in zip(points,kdistances):
        print keys[k],
        print [ ("%.3g"%t1, keys[t2]) for (t1,t2) in d]
    print
    
    #find and remove most crowded
    distancesAndKeys = []
    
    if points != extremes:
        for i,k in enumerate(points):
            if k in extremes:
                distancesAndKeys.append((10**300, k))
            else:
                distancesAndKeys.append((ksums[i], k))
    else:
        for i,k in enumerate(points):
            distancesAndKeys.append((ksums[i], k))      
    mcpoint = min(distancesAndKeys)[1]

    mck = keys[mcpoint]
    if print_removed:
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
    G=(0.25, 0.75), H=(0.25,0.5), I=(0.5,1), J=(0.5,0.75), K=(0.5,0.5), 
    L=(0.75,0.5))
    print 'initial'
    pprint(x)
    for n in range(15):
        x = removeMostCrowded(x, 3, print_removed=True)
        print '\n----------------------------\nafter removing %d points' % (n+1)
        pprint(x, showplot=False)

if __name__ == '__main__':
    test()
    pl.show()
