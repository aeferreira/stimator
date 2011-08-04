# -*- coding: ISO-8859-1 -*-

#------------------------------------------------------------------------------------------------------------------------------------
#Modified for the last time on May 15th, 2009


import numpy, random

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
    
    if len(x) < 3:
        x.popitem()
        return x
    else:
        keys = x.keys()
        keys.sort()
        extremes = []
        for i in range(len(random.choice(x.values()))):
            maximum = x[keys[0]][i]
            tempMax = [keys[0]]
            minimum = x[keys[0]][i]
            tempMin = [keys[0]]
            for j in keys:
                if x[j][i] < minimum:
                    tempMin = []
                    minimum = x[j][i]
                    tempMin.append(j)
                elif x[j][i] == minimum and j not in tempMin:
                    tempMin.append(j)                
                elif x[j][i] > maximum:
                    tempMax = []
                    maximum = x[j][i]
                    tempMax.append(j)
                elif x[j][i] == maximum and j not in tempMax:
                    tempMax.append(j)                
            for k in tempMin:
                if k not in extremes:
                    extremes.append(k)
            for k in tempMax:
                if k not in extremes:
                    extremes.append(k)
        for i in keys:
            x[i] = numpy.array(x[i])
        lista =[]
        #Use a fast implementation of the Euclidean distance
        temp = numpy.zeros(len(random.choice(x.values())))
        # Predefining temp allows reuse of this array, making this
        # function about twice as fast.
        distanceMatrix = []
        for i in range(len(keys)):
            distances = []  # list of (distance, index)
            for j in range(len(keys)):
                if j < i:
                    distances.append(distanceMatrix[j][i])
                elif j ==i:
                    distances.append(0)
                elif j > i:
                    temp[:] = x[keys[i]] - x[keys[j]]
                    dist = numpy.sqrt(numpy.dot(temp,temp))
                    distances.append(dist)
            distanceMatrix.append(numpy.copy(distances))
            distances.sort()
            if knumber < len(distances):
                lista.append(distances[1:knumber+1])
            if knumber >= len(distances):
                lista.append(distances[1:])
        distancesAndKeys = []
        tolook = []
        extremes.sort()
        if keys != extremes:
            for i in range(len(lista)):
                if keys[i] in extremes:
                    distancesAndKeys.append((10**300, keys[i]))
                    tolook.append(('inf', keys[i]))
                else:
                    distancesAndKeys.append(((sum(lista[i])), keys[i]))
                    tolook.append(((sum(lista[i])), keys[i]))
        else:
            for i in range(len(lista)):
                distancesAndKeys.append(((sum(lista[i])), keys[i]))        
        if pop_removed == True:
            print 'distancesAndKeys minimum', x[min(distancesAndKeys)[1]]
        del x[min(distancesAndKeys)[1]]
        for i in x.keys(): #Confirmar se este objecto deve ser retornado com listas ou arrays
            x[i] = list(x[i]) #de modo a ser usado pelas outras funções na geração seguinte.
        return x

#TEST

if __name__ == "__main__":
    xs = {1:[2,5],2:[3,4],3:[4,3],4:[4,5],5:[10,6],6:[11,4],7:[13,6],8:[13,3],9:[20,8],10:[22,5],11:[25,3],12:[25,8]}
    #xs = {1:[2,5],3:[4,3],8:[13,3],11:[25,3],12:[25,8]}
    while xs != {}:
        print removeMostCrowded(xs, 3, True)
    