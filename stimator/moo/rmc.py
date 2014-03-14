import numpy
import pylab as pl

def removeMostCrowded(x, objs, knumber = 3, remove_n = 1, verbose = False):

    """
    Removes a point with the smaller crowiding distance measured as
    the sum of the distances of the points to their k-nearest neighbours.
    x is a dictionary which values are the point coordinates as a list.
    """
    n_points = len(x)
    if n_points == 0 :
        return x
    if remove_n < 1:
        return x
    keys = x[:]
    n_objs = len(objs[0])
    
    points = range(n_points) #holds current indexes of points still present
    
    #compute matrix of objectives
    obj_matrix = numpy.empty((n_points, n_objs))
    for i,k in enumerate(keys):
        obj_matrix[i, :] = objs[i]

    if verbose:
        print 'mapping\n'
        mapping = [(k,i,obj_matrix[i]) for (k,i) in zip(keys, points)]
        for k, i, o in mapping:
            print k,i,o

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
    
    if verbose:
        print 'extreme points', [keys[i] for i in extremes]

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
    distances = []
    for i in points:
        dd = zip(list(distanceMatrix[i]), range(n_points))
        dd.sort()
        distances.append(dd)
    
    #compute k shortest (note: position 0 after sorting is always 0.0)
    last_removed = None
    
    # loop to remove each point
    for n in range(remove_n):
        if len(points) == 0:
            if verbose:
                print '\nNo more points to remove'
            return []

        if verbose:
            print '---------------------------\nRemoving point #%d' % (n+1), ':'
            print '\nremaining points' #, points
            print [keys[k] for k in points]
        
        if last_removed is not None:
            #remove reference to last removed point in list of distances
            for i in points:
                indx_last_remove = -1
                d = distances[i]
                for j in range(len(d)):
                    if d[j][1] == last_removed:
                        indx_last_remove = j
                        break
                del(distances[i][indx_last_remove])
        
        if verbose:
            print '\nk-distances'
            for i in points:
                dd = distances[i][1:knumber+1]
                print keys[i],
                print ', '.join([ "(%-5.3g to %s)"% (t1, keys[t2]) for (t1,t2) in dd])
            print
        
        ksums = [sum([d[0] for d in distances[i][1:knumber+1]]) for i in points]
        
        #find and remove most crowded
        distancesAndKeys = []
        
        allextremes = True
        for i in points:
            if i not in extremes:
                allextremes = False
                break

        if not allextremes:
            for i,k in enumerate(points):
                if k in extremes:
                    distancesAndKeys.append((10**300, k))
                else:
                    distancesAndKeys.append((ksums[i], k))
        else:
            for i,k in enumerate(points):
                distancesAndKeys.append((ksums[i], k))      
        
        last_removed = min(distancesAndKeys)[1]
        points.remove(last_removed)
        mc_key = keys[last_removed]
        
        if verbose:
            mcv = objs[last_removed]
            print 'Point to remove:', last_removed, mc_key, mcv
        
    return [keys[p] for p in points]


def pprint(x, objs, showplot=False):
    if len(x) == 0:
        print 'empty data'
    keys = x[:]
    mapping = [(k, i) for (k,i) in zip(keys, range(len(keys)))]
    mapping.sort()
    keys = [i for (i,j) in mapping]
    points = [j for (i,j) in mapping]
    for k, p in zip(keys, points):
        print k, objs[p]
    if showplot:
        xx = [objs[p][0] for p in points]
        yy = [objs[p][1] for p in points]
        pl.figure()
        pl.plot(xx, yy, 'bo')
        for label, x, y in zip(keys, xx, yy):
            pl.annotate(label, xy = (x, y), xytext = (x, y),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            arrowprops = None)
        pl.show()

def test():
    x = dict(A=(0,1), B=(0,0), C=(1,1), D=(1.2,0.5), E=(1,0), F=(0,0.5),
    G=(0.25, 0.75), H=(0.25,0.5), I=(0.5,1), J=(0.5,0.75), K=(0.5,0.5), 
    L=(0.75,0.5))
    keys = list(x.keys())
    objs = list(x.values())
    print 'initial data\n'
    pprint(keys, objs, showplot=True)
    x = removeMostCrowded(keys, objs, knumber = 3, remove_n = 15, verbose=True)
    print 'remaining points:'
    print x

if __name__ == '__main__':
    test()

