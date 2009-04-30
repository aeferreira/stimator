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

    def pprint(self):
        print
        print "the variables are" , self.variables
        print
        print "the constants are"
        for k in self.constants.keys():
               print "%s = %g" % (k, self.constants[k])
        print
        print "the parameters to find are"
        for k in self.parameters:
              print k[0],"from", k[1], "to", k[2]
        print
        print "the reactions are"
        for k in self.rates:
              irrstring = ""
              if k['irreversible']: irrstring = "(irreversible)"
              print k['name'], irrstring, ":"
              print " reagents:", k['reagents']
              print " products:", k['products']
              print " rate =", k['rate']
        print
        print "the @ definitions are"
        for k in self.atdefs:
              print "@", k[0], k[1], "=", k[2]
        print
        self.genStoichiometrymatrixOLD()
        print "the rows of the stoichiometry matrix are"
        for k,name in enumerate(self.variables):
              row = self.stoichmatrixrows[k]
              print "for", name, ":"
              for r in row.keys():
                  print "%g %s" % (row[r], r)
        print "the stoichiometry matrix as a numpy array is"
        N = self.genStoichiometrymatrix()
        print N
        print
        
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
