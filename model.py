#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator Model related classes
# Copyright António Ferreira 2006-2009
#----------------------------------------------------------------------------
import re
from numpy import *

class Model(object):
    def __init__(self):
        self.reset()
    
    def reset(self):
        self.problemname = ""   # the name of the problem

        self.rates       = []  # a list of {'name', 'reagents', 'products', 'rate', 'irreversible'}
        self.variables   = []  # a list of names of variables (order matters)
        self.constants   = {}  # a 'name':value dictionary
        self.parameters  = []  # a list of (name,min,max)
        self.atdefs      = []  # a list of (time,name,newvalue)

    def rateCalcString(self, rateString):
        #if self.error:
            #return ""
        nvars = len(self.variables)
        # replace varnames
        for i in range(nvars):
            rateString = re.sub(r"\b"+ self.variables[i]+r"\b", "variables[%d]"%i, rateString)
        # replace parameters
        for i in range(len(self.parameters)):
            rateString = re.sub(r"\b"+ self.parameters[i][0]+r"\b", "m_Parameters[%d]"%i, rateString)
        # replace constants
        for const in self.constants.keys():
            rateString = re.sub(r"\b"+ const + r"\b", "%e"% self.constants[const], rateString)
        return rateString
        
    def genStoichiometrymatrixOLD(self):
        self.stoichmatrixrows  = []   # sparse, using reactionname:coef dictionaries
        # build stoichiometry matrix row-wise
        for i in range(len(self.variables)):
            self.stoichmatrixrows.append({})
        for v in self.rates:
            for rORp, signedunit in [('reagents',-1.0),('products',1.0)]:
                for c in v[rORp]:
                    coef, var = (c[1]*signedunit, c[0])
                    if self.constants.has_key(var):
                        continue # there are no rows for constants in stoich. matrix
                    ivariable = self.variables.index(var) # error handling here
                    self.stoichmatrixrows[ivariable][v['name']] = coef

    def genStoichiometrymatrix(self):
        N = zeros((len(self.variables),len(self.rates)), dtype=float)
        for j,v in enumerate(self.rates):
            for rORp, signedunit in [('reagents',-1.0),('products',1.0)]:
                for c in v[rORp]:
                    coef, var = (c[1]*signedunit, c[0])
                    if self.constants.has_key(var):
                        continue # there are no rows for constants in stoich. matrix
                    ivariable = self.variables.index(var) # error handling here
                    N[ivariable, j] = coef
                    #self.stoichmatrixrows[ivariable][v['name']] = coef
        return N
