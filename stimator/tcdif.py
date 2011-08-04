# -*- coding: ISO-8859-1 -*-

#Neste ficheiro serao escritas as funçoes objectivo alternativas, i.e. L2, ...

import analysis
from numpy import log, transpose, array, mean, std, float64, random, copy, diag, dot, linalg, trace, pi, append, abs, linspace, nansum, isnan, compress
from numpy.random import uniform
from numpy.linalg import inv, det, eigvals
from scipy import stats
from scipy.stats import norm
import sympy
from sympy import Symbol, diff
import re
from model import variable, parameters, variables, varnames
from dynamics import state2array

def generateBias(curveNumber, variableNumber, standardDeviation):
    bias = []
    for i in range(variableNumber):
        bias.append(norm.rvs(0,standardDeviation,size = curveNumber))
    return transpose(bias)

def generateUniformMeasurementError(absoluteError, npoints):
    errors = []
    for i in range(npoints):
        errors.append(uniform(-absoluteError,absoluteError))
    return tuple(errors)

def objDif(models, trial, t0, npoints, tf, objFunc, bias, measurementErrors):
    difs = []
    lenModels = len(models)
    i = 0
    while i < lenModels - 1:
        j = i + 1
        difs.append([])
        while j < lenModels:
            obj1, obj2 = Objective(models[i], t0, npoints, tf, objFunc, bias, measurementErrors), Objective(models[j], t0, npoints, tf, objFunc, bias, measurementErrors)
            difs[-1].append(obj1(trial) - obj2(trial))
            j += 1
        i += 1
    return difs

## def KLDiscrepancies(modelTCs, deltaT):
##     plusKLlist = []
##     minusKLlist = []
##     for i in range(len(modelTCs)-1):
##         for j in range(i+1, len(modelTCs)):
##             plusKL = 0.0
##             minusKL = 0.0
##             for m,n in zip(modelTCs[i],modelTCs[j]):
##                 plusTC = (m*(log(m/n)+n/m-1)*deltaT)
##                 plusKL -= nansum(plusTC)
##                 minusTC = (n*(log(n/m)+m/n-1)*deltaT)
##                 minusKL -= nansum(minusTC)
##             plusKLlist.append(plusKL)
##             minusKLlist.append(minusKL)
##     result = plusKLlist + minusKLlist
##     return result

## def KLs(modelTCs, deltaT):
##     plusKLlist = []
##     minusKLlist = []
##     for i in range(len(modelTCs)-1):
##         for j in range(i+1, len(modelTCs)):
##             plusKL = 0.0
##             minusKL = 0.0
##             for m,n in zip(modelTCs[i],modelTCs[j]):
##                 plusTC = float64(m*(log(m/n))*deltaT)
##                 plusKL -= nansum(plusTC)#[plusTC >=0.0])
##                 minusTC = float64(n*log(n/m)*deltaT)
##                 minusKL -= nansum(minusTC)#[minusTC >=0.0])
##             plusKLlist.append(plusKL)
##             minusKLlist.append(minusKL)
##     result = plusKLlist + minusKLlist
##     return result

def KLDiscrepancies(modelTCs, deltaT):
    plusKLlist = []
    minusKLlist = []
    for i in range(len(modelTCs)-1):
        for j in range(i+1, len(modelTCs)):
            m = modelTCs[i].data
            n = modelTCs[j].data
            plusKL = -deltaT * nansum(float64(m*(log(m/n)+n/m-1)))
            minusKL = -deltaT * nansum(float64(n*(log(n/m)+m/n-1)))
            plusKLlist.append(plusKL)
            minusKLlist.append(minusKL)
    result = plusKLlist + minusKLlist
    return result

def KLs(modelTCs, deltaT):
    plusKLlist = []
    minusKLlist = []
    for i in range(len(modelTCs)-1):
        for j in range(i+1, len(modelTCs)):
            m = modelTCs[i].data
            n = modelTCs[j].data
            plusKL = -deltaT * nansum(float64(m*log(m/n)))
            minusKL = -deltaT * nansum(float64(n*log(n/m)))
            plusKLlist.append(plusKL)
            minusKLlist.append(minusKL)
    result = plusKLlist + minusKLlist
    return result

def kremling(modelTCs, deltaT):
    #Maximizing this function is basically the same as maximizing a weighted L2 distance.
    result = []
    for i in range(len(modelTCs)-1):
        for j in range(i+1, len(modelTCs)):
            numResult = 0.0
            for tc1,tc2 in zip(modelTCs[i],modelTCs[j]):
                tempTC = float64((((tc1-tc2)**2)/(((tc1+tc2)/2)**2))*deltaT)
                numResult -= nansum(tempTC)
            result.append(numResult)
    return result

def L2(modelTCs, deltaT):
    #Maximizes this function is the same as maximizing the L2 distance.
    result = []
    for i in range(len(modelTCs)-1):
        for j in range(i+1, len(modelTCs)):
            numResult = 0.0
            for tc1,tc2 in zip(modelTCs[i],modelTCs[j]):
                tempTC = float64(((tc1-tc2)**2))*deltaT
                numResult -= nansum(tempTC)
            result.append(numResult)
    return result


class Objective:
    
    """Computes the objective function for a model and a trial vector 
    given the initial and final time points for the computation, the number of points for ODE integration,
    a user-defined bias for the initial conditions and measurement errors for each integration point.
    'model' is a model object;
    't0' and 'tf' are floats;
    'npoints' is a positive integer;
    'objFunc' is the name (a string) of the objective function to be used;
    'bias' is a list with the relative changes to be applied to the initial conditions;
    'measurementErrors' is a list of absolute changes to be added to the time course values."""
    
    def __init__(self, model, t0, npoints, tf, objFunc, optnames, observed, bias = None, measurementErrors = None):
        self.t0 = t0
        self.npoints = npoints
        self.tf = tf
        self.vector = state2array(model,"init")#model.vectorize("init")

        self.lenEstimatedPar = len(parameters(model)) #len(model.parameters)
        self.objFunc = objFunc
        
##         self.varNames = [i.name for i in model.variables]
        self.varNames = varnames(model)
        
        self.model = model
        self.optnames = optnames
        self.observed = observed
        
        self.optvars = [name for name in self.varNames if name in optnames]
        self.optvarsindexes = array([self.varNames.index(name) for name in self.optvars])
        
        self.obsIndex = []
        self.ordObsVarNames = []
##         for xname in self.varNames:
##             if xname in self.observed:
##                 self.obsIndex.append(self.varNames.index(xname))
##                 self.ordObsVarNames.append(xname)
        for (i,xname) in enumerate(self.varNames):
            if xname in self.observed:
                self.obsIndex.append(i)
                self.ordObsVarNames.append(xname)
        
        if self.objFunc in ('AIC', 'AICc', 'AICu', 'criterionA', 'modCriterionA', 'criterionD', 'criterionE', 'modCriterionE'):
            self.objective = getattr(self, objFunc)
            self.bias = bias
            self.measurementErrors = measurementErrors
    
    def __call__(self, trial):
        """Returns the result of the evaluation of the objective function for the particular trial candidate solution."""
        self.trial = trial # trial is a dicionary of optname:value items
                    
##         self.initTrial = []
        
##         for j in range(len(self.varNames)):
##             if self.varNames[j] in self.optvars:
##                 self.initTrial.append(self.trial[self.varNames[j]])
##             else:
##                 self.initTrial.append(self.vector[j])
##         self.initTrialArray = array(self.initTrial)

##         self.initTrialArray = copy(self.vector)
##         trialvec = array([self.trial[name] for name in self.optvars])
##         self.initTrialArray[self.optvarsindexes] = trialvec

        self.initTrialArray = copy(self.vector)
        for name,i in zip(self.optvars, self.optvarsindexes):
            self.initTrialArray[i] = value = self.trial[name]

        self.regularSimulation = analysis.solve(self.model, tf = self.tf, npoints = self.npoints, t0 = self.t0, initial = self.initTrialArray)
        
        if self.objFunc in ['KL','kremling','KLs','L2']:
            #Notice: the solution is returned and the actual objective is computed in the caller
            self.result = self.regularSimulation.copy(self.ordObsVarNames)
##             self.result = []
##             index = 0
##             for i in self.regularSimulation:
##                 if index in self.obsIndex:
##                     self.result.append(copy(i))
##                 index += 1
        elif self.objFunc in ['AIC', 'AICc', 'AICu']:
            self.result = self.objective(self.model, self.bias, self.measurementErrors)
        elif self.objFunc in ['criterionA', 'modCriterionA', 'criterionD', 'criterionE', 'modCriterionE']:
            self.lenVariables = len(self.varNames)
            self.parNames = []
            for i in range(self.lenEstimatedPar):
                self.parNames.append(self.model.parameters[i].name)
            self.senses = self.calculateSensitivities(self.model)
            self.sensModel = self.model.clone()
            self.addVars(self.sensModel, self.senses)
            self.TCandSenses = self.calculateSValues(self.sensModel)
            self.TCs, self.sensValues = self.separateTCFromS(self.sensModel, self.TCandSenses, self.lenVariables)
            self.sensMatrixNames, self.sensMatrixValues = self.calculateSMatrix(self.sensModel, self.sensValues, len(self.ordObsVarNames))
            self.vcMatrix = self.calculateVCMatrix(self.model, self.bias, self.measurementErrors)
            self.fim = self.calculateFIM(self.sensMatrixValues, self.vcMatrix)
            self.result = self.objective(self.fim)
        #TODO: implement the following objective functions in such a way that the sum of the derivates for each observed variable is maximized by Pareto MOO
        elif self.objFunc in ['dss','das']:
            transformation = self.model.getdXdt()
            self.regularSimulation.apply_transf(transformation)
            self.result = 0
            if self.objFunc == 'dss':
                for i in self.regularSimulation:
                    self.result += sum(num * num for num in i)
            elif self.objFunc == 'das':
                for i in self.regularSimulation:
                    self.result += sum(abs(num) for num in i)
        return self.result
    
    #Akaike's informaion criteria methods
    
    def generateRSS(self, model, bias, measurementErrors, returnList = False):
        """Generates the sum of the squared residuals from a user-defined bias (for the initial variable values) 
        and measurement errors"""
        alternativeSolutions = []
        for i in range(len(bias)):
            alternativeSolutions.append([])
            for j in range(len(bias[i])): #The length of bias[i] must be the same as the length of trialKeys, which is the number of variables to be optimized.
                alternativeSolutions[-1].append(self.trial[self.optvars[j]]*bias[i][j]+self.trial[self.optvars[j]])
        perturbedSimulations = []
        j = 0
        while j < len(bias):
            tempInit = []
            m = 0
            for i in (range(len(self.varNames))):
                if self.varNames[i] in self.optvars:
                    tempInit.append(alternativeSolutions[j][m])
                    m += 1
                else:
                    tempInit.append(self.vector[i])
            tempInit = array(tempInit)
            integrationResult = analysis.solve(model, tf = self.tf, npoints = self.npoints, t0 = self.t0, initial = tempInit)#[0]
            perturbedSimulations.append(integrationResult) 
            j = j + 1
        perturbedSimulations = transpose(perturbedSimulations)
        residuals=[]
        for i in range(len(perturbedSimulations)):
            residuals.append([])
            for j in range(len(perturbedSimulations[i])):
                if j in self.obsIndex:
                    residuals[-1].append(mean(perturbedSimulations[i][j])+measurementErrors[i]-self.regularSimulation[j][i])
        rss=0
        for i in residuals:
            for j in i:
                rss = rss + j**2
        if returnList == False:
            return rss
        elif returnList == True:
            return residuals

    def AIC(self, model, bias, measurementErrors):
        return 2*self.lenEstimatedPar+self.npoints*(log(2*pi*self.generateRSS(model, bias, measurementErrors)/self.npoints)+1)

    def AICc(self, model, bias, measurementErrors):
        return log(self.generateRSS(model, bias, measurementErrors)/self.npoints)+(self.npoints + self.lenEstimatedPar)/(self.npoints - self.lenEstimatedPar - 2)
    
    def AICu(self, model, bias, measurementErrors):
        return log(self.generateRSS(model, bias, measurementErrors)/(self.npoints - self.lenEstimatedPar))+(self.npoints + self.lenEstimatedPar)/(self.npoints - self.lenEstimatedPar - 2)

    #Fisher information matrix methods

    def calculateSensitivities(self, model):
        '''Returns a list of tuples with form 
        [(name of the sensitivity of a variable to a parameter, right-hand expression of the ODE defining the sensitivity as a function of time)].'''
        varSymb, parSymb = [], []
        varDer, parDer = [], []
        dxdtStrings = []
        sensitivityNames, sensitivityMatrix = [], []
        for i in range(self.lenVariables):
            if i in self.obsIndex:
                varSymb.append(Symbol(self.varNames[i]))
        for i in range(self.lenEstimatedPar):
            parSymb.append(Symbol(self.parNames[i]))
        strings = model.dXdt_strings()
        counter = 0
        for i,j in strings:
            if counter in self.obsIndex:
                tempString = j
                for k in self.varNames:
                    tempString = re.sub(r'\b'+k+r'\b',"Symbol('"+k+"')",tempString)
                for k in self.parNames:
                    tempString = re.sub(r'\b'+k+r'\b',"Symbol('"+k+"')",tempString)
                dxdtStrings.append(eval(tempString))
            counter += 1
        for i in range(len(dxdtStrings)):
            varDer.append('('+str(sympy.diff(dxdtStrings[i],varSymb[i]))+')')
            parDer.append([])
            for j in parSymb:
                parDer[-1].append('('+str(sympy.diff(dxdtStrings[i],j))+')')
        for i,j,k in zip(varDer, parDer, self.ordObsVarNames):
            for m in range(len(j)):
                sensitivityNames.append('S'+k+self.parNames[m])
                sensitivityMatrix.append(i+'*S'+k+self.parNames[m]+'+'+j[m])
        return zip(sensitivityNames, sensitivityMatrix) 

    def addVars(self, model, varList):
        '''Adds variables contained in a list of tuples of the form 
        [(name of a variable, right-hand expression of the ODE defining the variable as a function of time), ...] to a model .'''
        self.sensVector = self.vector
        self.sensVarNames = list(copy(self.varNames))
        for i,j in varList:
            setattr(model, i, variable(j))
            self.sensVector = append(self.sensVector,0)
            self.sensVarNames.append(i)

    def calculateSValues(self, model):
        """Calculates the time courses and the 
        time-dependent sensitivities of the variables to the parameters."""
        for j in range(len(self.ordObsVarNames)*self.lenEstimatedPar):
            self.initTrialArray = append(self.initTrialArray, 0)
        globalResult = analysis.solve(model, tf = self.tf, npoints = self.npoints, t0 = self.t0, initial = self.initTrialArray)
        return zip(self.sensVarNames, globalResult)

    def separateTCFromS(self, model, instantMatrix, TCNumber):
        """Separates the time courses from the sensitivities 
        and returns both in 2-element tuple."""
        return instantMatrix[:TCNumber], instantMatrix[TCNumber:]

    def calculateSMatrix(self, model, instantSMatrix, lenOrdObsVarNames):
        """Calculates the sensitivity matrix 
        in the form it must be to enter the calculation of the FIM."""
        globalSums = []
        for i,j in instantSMatrix:
            s = 0
            for k in j:
                s += abs(k)
            globalSums.append(s)
        k = 0
        globalMatrix = []
        nameMatrix = []
        for i in range(lenOrdObsVarNames):
            globalMatrix.append([])
            nameMatrix.append([])
            for j in model.parameters:
                globalMatrix[-1].append(globalSums[k])
                nameMatrix[-1].append(instantSMatrix[k][0])
                k += 1
        return nameMatrix, globalMatrix

    def calculateVCMatrix(self, model, bias, measurementErrors, calcCovariances = False):
        """Calculates the variance-covariance matrix that enters in the calculation of the FIM."""
        residuals = self.generateRSS(model, bias, measurementErrors, returnList = True)
        variances = []
        for i in transpose(residuals):
            s = 0
            for j in i:
                s += j**2
            variances.append(s/(self.npoints - self.lenEstimatedPar))
        if calcCovariances == False:
            return diag(variances)
        elif calcCovariances == True:
            print "Don't know how to calculate covariances yet..."

    def calculateFIM(self, sensMatrix, vcMatrix):
        """Computes the Fisher information matrix from the sensitivity matrix 
        and the variance-covariance matrix of the residuals."""
        return dot(dot(transpose(sensMatrix),vcMatrix),sensMatrix)

    def criterionA(self, fim):
        return trace(inv(fim))

    def modCriterionA(self, fim):
        return trace(fim)

    def criterionD(self, fim):
        print 'fim', fim
        return -det(fim)

    def criterionE(self, fim):
        result = abs(min(eigvals(fim)))
        return result

    def modCriterionE(self, fim):
        maxEV = max(eigvals(fim))
        minEV = min(eigvals(fim))
        return abs(maxEV/minEV)
    

if __name__ == "__main__":
    tcs = [[array([1,2,3,4,5])],[array([6.0,8.0,10.0,21.0,30.0])]]
    hr = 0.5
    print KLs(tcs,hr)
    print kremling(tcs,hr)
    print L2(tcs,hr)
