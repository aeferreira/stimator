import numpy
import pylab as pl

def remove_most_crowded(coords, 
                        knumber = 3, remove_n = 1, 
                        indexes = None, labels = None, 
                        verbose = False):
    """Remove points with the smaller crowiding distance from a set.

    The crowding distance is measured as the sum of the distances 
    of the points to their k-nearest neighbours.
    
    coords -- a list or 2D numpy array of point coordinates
    remove_n -- the number of points to remove
    knumber -- the number of neighbours to consider in the crowding distance
    indexes -- indexes (list of integers) of points in coords to consider
    (Default(None): all points are considered)
    labels -- a list of strings (with the same size as indexes) used to label
    points in graphs and in verbose output. (Default: string representation of
    indexes)
    
    Inside the function, coords, indexes and labels are not changed,
    only indexed.
    
    Return is a list of integers (a subset of indexes) representing the
    remaining points after removal (can be used to index coords).
    """
    
    if indexes is None:
        indexes = range(len(coords))
    if labels is None:
        labels = [str(i) for i in indexes]
    n_points = len(indexes)
    if n_points == 0 or remove_n < 1:
        return indexes
    n_coords = len(coords[indexes[0]])
    
    points = range(len(indexes)) #holds current points still present

##     if verbose:
##         print 'mapping\n'
##         for k,i,o in [(labels[p], indexes[p], coords[indexes[p]]) for p in points]:
##             print k, i, o
    
    #unfortunatelly, a copy is necessary
    acoords = numpy.array([coords[i] for i in indexes])
    
    # find indexes of extreme points (max and min in each dimension)
    maxima = numpy.amax(acoords, axis=0)
    minima = numpy.amin(acoords, axis=0)
    
    extremes = []
    for i in range(n_coords):
        dimextremes = []
        for p in points:
            v = acoords[p][i]
            if v == minima[i] or v == maxima[i]:
                dimextremes.append(p)                
        for p in dimextremes:
            if p not in extremes:
                extremes.append(p)
    extremes.sort()
    
    if verbose:
        print '\nextreme points', [labels[x] for x in extremes]

    #compute distances
    # TODO: use scipy function
    distanceMatrix = []
    for i in points:
        distances = []
        for j in points:
            if j < i:
                distances.append(distanceMatrix[j][i])
            elif j ==i:
                distances.append(0.0)
            elif j > i:
                temp = acoords[i,:] - acoords[j,:]
                d = numpy.sqrt(numpy.dot(temp,temp))
                distances.append(d)
        distanceMatrix.append(distances)
    
    #sort distances for each point
    distances = []
    for p in points:
        dd = zip(list(distanceMatrix[p]), points)
        dd.sort()
        distances.append(dd)
    if verbose:
        print '\nDistances'
        for p in points:
            dd = distances[p]
            print labels[p],
            print ', '.join(["(%-5.3g to %s)"% (t1, labels[t2]) for (t1,t2) in dd])
        print
    
    #compute k shortest (note: position 0 after sorting is always 0.0)
    last_removed = None
    
    # loop to remove each point
    for n in range(remove_n):
        if len(points) == 0:
            if verbose:
                print '\nNo more points to remove'
            return points

        if verbose:
            print '---------------------------\nRemoving point #%d' % (n+1), ':'
            print '\nremaining points' #, points
            print [labels[p] for p in points]
        
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
                print labels[i],
                print ', '.join([ "(%-5.3g to %s)"% (t1, labels[t2]) for (t1,t2) in dd])
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
                if k not in extremes:
                    distancesAndKeys.append((ksums[i], k))
        else:
            for i,k in enumerate(points):
                distancesAndKeys.append((ksums[i], k))      
        
        last_removed = min(distancesAndKeys)[1]
        points.remove(last_removed)
        mc_key = labels[last_removed]
        mc_index = indexes[last_removed]
        
        if verbose:
            mcv = coords[mc_index]
            print 'Point to remove:', mc_key, mcv
        
    return [indexes[p] for p in points]


def pprint(coords, indexes = None, labels=None, showplot=False):
    if indexes is None:
        indexes = range(len(coords))
    if len(indexes) == 0:
        print 'empty data'
    if labels is None:
        labels = [str(i) for i in indexes]
    for k,p in zip(labels, indexes):
        print k, '(%d)'%p, coords[p]
    if showplot:
        plotlabels = labels
        xx = [coords[p][0] for p in indexes]
        yy = [coords[p][1] for p in indexes]
        pl.figure()
        pl.plot(xx, yy, 'bo')
        for label, x, y in zip(plotlabels, xx, yy):
            pl.text(x+0.01, y, label, 
                    ha = 'left', 
                    va = 'center',
                    bbox=dict(facecolor='white', alpha=0.5))
        pl.show()

def test():
    x = dict(A=(0,1), B=(0,0), Bogus=(-1,-1),
             C=(1,1), D=(1.2,0.5), E=(1,0), 
             F=(0,0.5), G=(0.25, 0.75), H=(0.25,0.5), 
             I=(0.5,1), J=(0.5,0.75), K=(0.5,0.5), L=(0.75,0.5))
    keys = list(x.keys())
    keys.sort()
    coords = [x[k] for k in keys]
    print 'all labels', keys
    points = range(len(coords))
    print 'all points', points
    
    #remove the bogus point
    bindx = keys.index('Bogus')
    del(points[bindx])
    del(keys[bindx])
    print
    print 'labels', keys
    print 'use only points'
    print points


    print
    print 'initial data\n'
    pprint(coords, indexes=points, labels=keys, showplot=True)
    
    x = remove_most_crowded(coords, 
                            indexes=points,
                            labels=keys,
                            knumber=3,
                            remove_n=15,
                            verbose=True)
    print 'remaining points:'
    print x

if __name__ == '__main__':
    test()

