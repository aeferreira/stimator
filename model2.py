#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator Model related classes
# Copyright António Ferreira 2006-2009
#----------------------------------------------------------------------------
import os
import os.path
import re
import math
import timecourse
import utils
from numpy import *

#----------------------------------------------------------------------------
#         Functions to check the validity of math expressions
#----------------------------------------------------------------------------
def test_with_consts(valueexpr, consts={}):
    """Uses builtin eval function to check for the validity of a math expression.

       Constants previously defined can be used"""

    locs = consts.copy()
    try :
       value = float(eval(valueexpr, vars(math), locs))
    except Exception, e:
       return ("%s : %s"%(str(e.__class__.__name__),str(e)), 0.0)
    return ("", value)


def test_with_everything(valueexpr, parameters, unknown, varlist):
    locs = {}
    for p in parameters:
        locs[p.name] = p.value
    
    #part 1: nonpermissive, except for NameError
    try :
       value = float(eval(valueexpr, vars(math), locs))
    except NameError:
       pass
    except Exception, e:
       return ("%s : %s"%(str(e.__class__.__name__),str(e)), 0.0)
    #part 2: permissive, with dummy values (1.0) for vars and unknown parameters
    vardict = {}
    for i in varlist:
        vardict[i.name]=1.0
    for i in unknown:
        vardict[i.name]=1.0
    locs.update(vardict)
    try :
       value = float(eval(valueexpr, vars(math), locs))
    except (ArithmeticError, ValueError):
       pass # might fail but we don't know the values of vars
    except Exception, e:
       return ("%s : %s"%(str(e.__class__.__name__),str(e)), 0.0)
    return "", value

#----------------------------------------------------------------------------
#         Regular expressions for grammar elements and dispatchers
#----------------------------------------------------------------------------
identifierpattern = r"[_a-z]\w*"
fracnumberpattern = r"[-]?\d*[.]?\d+"
realnumberpattern = fracnumberpattern + r"(e[-]?\d+)?"
stoichiompattern   = r"^\s*(?P<reagents>.*)\s*(?P<irreversible>->|<=>)\s*(?P<products>.*)\s*$"

chemcomplexpattern = r"^\s*(?P<coef>\d*)\s*(?P<variable>[_a-z]\w*)\s*$"

nameErrorpattern = r"NameError : name '(?P<name>\S+)' is not defined"

identifier = re.compile(identifierpattern, re.IGNORECASE)
fracnumber = re.compile(fracnumberpattern, re.IGNORECASE)
realnumber = re.compile(realnumberpattern, re.IGNORECASE)
stoichiom  = re.compile(stoichiompattern,    re.IGNORECASE)

chemcomplex = re.compile(chemcomplexpattern, re.IGNORECASE)

nameErrormatch = re.compile(nameErrorpattern)

#----------------------------------------------------------------------------
#         Utility functions
#----------------------------------------------------------------------------
def processStoich(expr):
    match = stoichiom.match(expr)
    if not match:
        return None

    #process irreversible
    irrsign = match.group('irreversible')
    irreversible = irrsign == "->"
    reagents = []
    products = []

    #process stoichiometry
    fields = [(reagents,'reagents'),(products,'products')]
    for target,f in fields:
        complexesstring = match.group(f).strip()
        if len(complexesstring)==0:  #empty complexes allowed
            continue
        complexcomps = complexesstring.split("+")
        for c in complexcomps:
            m = chemcomplex.match(c)
            if m:
               coef = m.group('coef')
               var = m.group('variable')
               if coef == "":
                  coef = 1.0
               else:
                  coef = float(coef)
               if coef == 0.0: continue # a coef equal to zero means ignore
               target.append((var,coef))
            else:
               return None

    return reagents, products, irreversible

def massActionStr(k = 1.0, reagents = []):
    res = str(k)
    factors = []
    for var, coef in reagents:
        if coef == 0.0:
            factor = ''
        if coef == 1.0:
            factor = '%s'% var
        else:
            factor = '%s**%f' % (var,coef)
        factors.append(factor)
    factors = '*'.join(factors)
    if factors != '':
        res = res + '*' + factors
    return res

def findWithName(name, alist):
    res = None
    for i in alist:
        if i.name == name:
            return i
#----------------------------------------------------------------------------
#         Model and Model component classes
#----------------------------------------------------------------------------
class Reaction(object):
    def __init__(self, reagents, products, rate = 0.0, irreversible = False):
        self.reagents = reagents
        self.products = products
        self.rate = rate
        self.irrversible = irreversible
        self.name = '?'
    def __str__(self):
        return "\n%s:\n  reagents: %s\n  products: %s\n  rate = %s" % (self.name, str(self.reagents), str(self.products), str(self.rate)) 

class State(object):
    pass

class Transformation(object):
    def __init__(self, rate = 0.0):
        self.rate = rate
        self.name = '?'
    def __str__(self):
        return "\n%s:\n  rate = %s" % (self.name, str(self.rate))

class Constant(object):
    def __init__(self, value = 0.0):
        self.value = value
        self.name = '?'
    def __str__(self):
        return "\n%s = %f" % (self.name, self.value)

class PairConstants(object):
    def __init__(self, min = 0.0, max = 1.0):
        self.min = min
        self.max = max
        self.name = '?'
    def __str__(self):
        return "\n%s = ? (min=%f, max=%f)" % (self.name, self.min, self.max)

class Variable(object):
    def __init__(self, name):
        self.name = name
    def __str__(self):
        return self.name

class Model(object):
    def __init__(self, title = ""):
        self.__dict__['_Model__reactions']         = []
        self.__dict__['_Model__variables']         = []
        self.__dict__['_Model__extvariables']      = []
        self.__dict__['_Model__parameters']        = []
        self.__dict__['_Model__unknownparameters'] = []
        self.__dict__['_Model__transf']            = []
        self.__dict__['_Model__uncertain']         = []
        self.__title = title
    
    def __setattr__(self, name, value):
        if isinstance(value, Reaction):
            value.name = name
            self.__dict__['_Model__reactions'].append(value)
        elif isinstance(value, float) or isinstance(value, int):
            obj = Constant(float(value))
            obj.name = name
            self.__dict__['_Model__parameters'].append(obj)
        elif isinstance(value, Transformation):
            value.name = name
            self.__dict__['_Model__transf'].append(value)
        elif (isinstance(value, tuple) or isinstance(value, list)) and len(value)==2:
            obj = PairConstants(float(value[0]), float(value[1]))
            obj.name = name
            self.__dict__['_Model__unknownparameters'].append(obj)
        else:
            object.__setattr__(self, name, value)
        self.__refreshVars()

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        c = findWithName(name, self.__parameters)
        if c :
            return c.value
        c = findWithName(name, self.__reactions)
        if c :
            return c
        c = findWithName(name, self.__variables)
        if c :
            return c.name
        c = findWithName(name, self.__transf)
        if c :
            return c
        raise AttributeError, name + ' is not defined for this model'
    
    def varnames(self):
        return [i.name for i in self.__variables]
    
    def checkRates(self):
        for collection in (self.__reactions, self.__transf):
            for v in collection:
                resstring, value = test_with_everything(v.rate, self.__parameters, self.__unknownparameters, self.__variables)
                if resstring != "":
                    return False, resstring + '\nin rate of %s\n(%s)' % (v.name, v.rate)
        return True, 'OK'

    def rateCalcString(self, rateString):
        # replace varnames
        for i,v in enumerate(self.__variables):
            rateString = re.sub(r"\b"+ v.name+r"\b", "variables[%d]"%i, rateString)
        # replace parameters
        for i,u in enumerate(self.__unknownparameters):
            rateString = re.sub(r"\b"+ u.name+r"\b", "m_Parameters[%d]"%i, rateString)
        # replace constants
        for p in self.__parameters:
            rateString = re.sub(r"\b"+ p.name + r"\b", "%g"% p.value, rateString)
        return rateString

    def __str__(self):
        res = "%s\n"% self.__title
        for collection in (self.__reactions, self.__transf, self.__parameters,  self.__unknownparameters):
            for i in collection:
                res += str(i)
            res+='\n'
        res += "\nVariables: " + " ".join([i.name for i in self.__variables])
        res+='\n'
        res += "External variables: " + " ".join([i.name for i in self.__extvariables])
        res+='\n\n'
        check, msg = self.checkRates()
        if not check:
            res +="Problem in rates:\n"
            res += msg
        else:
            res +='rates are OK'

        return res
    
    def __refreshVars(self):
        del(self.__variables[:]) #can't use self.__variables= [] : Triggers __setattr__
        del(self.__extvariables[:])
        for v in self.__reactions:
            for rp in (v.reagents, v.products):
                for (name, coef) in rp:
                    if findWithName(name, self.__variables):
                        continue
                    else:
                        if findWithName(name, self.__parameters):
                            self.__extvariables.append(Variable(name))
                        else:
                            self.__variables.append(Variable(name))


    def genStoichiometryMatrix(self):
        varnames = self.varnames()
        N = zeros((len(self.__variables),len(self.__reactions)), dtype=float)
        for j,v in enumerate(self.__reactions):
            for rORp, signedunit in [(v.reagents,-1.0),(v.products,1.0)]:
                for c in rORp:
                    coef, var = (c[1]*signedunit, c[0])
                    if var in varnames:
                        ivariable = varnames.index(var) # error handling here
                        N[ivariable, j] = coef
                    else:
                        continue # there are no rows for extvariables in stoich. matrix
        return N

    def genRatesFunction(self, scale = 1.0):

        #compile rate laws
        self.__ratebytecode = [compile(self.rateCalcString(v.rate), 'bof.log','eval') for v in self.__reactions]
        
        # create array to hold v's
        self.__v = empty(len(self.reactions))
            
        def calcDerivs(variables, unknownparameters, t):
            m_Parameters = unknownparameters
            ratebytecode = self.__ratebytecode
            v = self.__v
            for i,r in enumerate(ratebytecode):
                v[i] = eval(r)
            return v
        return calcDerivs
        

    def __getReactions(self):
        return self.__reactions
    def __getVariables(self):
        return self.__variables
    def __getExtVariables(self):
        return self.__extvariables
    def __getParameters(self):
        return self.__parameters
    def __getUnknownParameters(self):
        return self.__unknownparameters
    def __getTransformations(self):
        return self.__transf
    
    variables    = property(__getVariables)
    extvariables = property(__getExtVariables)
    reactions    = property(__getReactions)
    parameters   = property(__getParameters)
    unknown      = property(__getUnknownParameters)
    transf       = property(__getTransformations)

class BadStoichError(Exception):
    """Used to flag a wrong stoichiometry expression"""

def react(stoichiometry, rate = 0.0):
    res = processStoich(stoichiometry)
    if not res:
        raise BadStoichError, "Bad stoichiometry definition:\n"+ stoichiometry
    if isinstance(rate, float) or isinstance(rate, int):
        rate = massActionStr(rate, res[0])
    return Reaction(res[0], res[1], rate, res[2])

def transf(rate = 0.0):
    if isinstance(rate, float) or isinstance(rate, int):
        rate = str(rate)
    return Transformation(rate)

def main():
    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   -> "   , "V3 * C / (Km3+ C)")
    m.t1 = transf("A*4 + C")
    m.B  = 2.2
    m.myconstant = 4
    m.V3 = (0.1, 0.9)
    m.Km3 = 4
    print m
    
    #~ print m.B * m.myconstant
    #~ print m.varnames()
    #~ print m.v1, '\nin the sky with', m.A
    #~ print m.v3.rate
    #~ print m.t1.rate
    print '\nStoichiometry matrix:'
    N = m.genStoichiometryMatrix()
    print '  ', '  '.join([v.name for v in m.reactions])
    for i,x in enumerate(m.variables):
        print x.name, N[i, :]
    
    print
    print 'calcstring of v3:'
    print m.rateCalcString(m.v3.rate)
    print
    vrates = m.genRatesFunction()
    vars = [0, 1]
    pars = [1.0]
    result = vrates(vars,pars,0)
    
    print 'with'
    for v,value in zip(m.variables, vars):
        print v.name, '=', value
    for v,value in zip(m.unknown, pars):
        print v.name, '=', value
    for p in m.parameters:
        print p.name, '=', p.value
    print 'then...'
    for v,r in zip(m.reactions, result):
        print v.name, '=', v.rate, '=', r
    
if __name__ == "__main__":
    main()
 