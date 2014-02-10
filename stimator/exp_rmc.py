import numpy
import pylab as pl

def removeMostCrowded(x, knumber = 3, remove_n = 1, verbose = False):

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
    keys = list(x.keys())
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
            return x

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
            mcv = x[mc_key]
            print 'Point to remove:', mc_key, mcv
        
        del x[mc_key]
    
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
        pl.show()

def test():
    x = dict(A=(0,1), B=(0,0), C=(1,1), D=(1.2,0.5), E=(1,0), F=(0,0.5),
    G=(0.25, 0.75), H=(0.25,0.5), I=(0.5,1), J=(0.5,0.75), K=(0.5,0.5), 
    L=(0.75,0.5))
    print 'initial dictionary\n'
    pprint(x, showplot=True)
    x = removeMostCrowded(x, 3, 15, verbose=True)

if __name__ == '__main__':
    test()

