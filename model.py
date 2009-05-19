#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator Model related classes
# Copyright Ant�nio Ferreira 2006-2009
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

def test_with_everything(valueexpr, parameters, unknown, varlist):
    locs = {}
    for p in parameters:
        locs[p.name] = p #.value
    
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
stoichiompattern   = r"^\s*(?P<reagents>.*)\s*(?P<irreversible>->|<=>)\s*(?P<products>.*)\s*$"
chemcomplexpattern = r"^\s*(?P<coef>\d*)\s*(?P<variable>[_a-z]\w*)\s*$"

nameErrorpattern = r"NameError : name '(?P<name>\S+)' is not defined"

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
    for i in alist:
        if i.name == name:
            return i
    return None

def findWithNameIndex(name, alist):
    for i,elem in enumerate(alist):
        if elem.name == name:
            return i
    return -1
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
        return "%s:\n  reagents: %s\n  products: %s\n  rate    = %s\n" % (self.name, str(self.reagents), str(self.products), str(self.rate)) 

class StateArray(list):
    def __init__(self, data, varvalues, name):
        list.__init__(self, data)
        self.__dict__['varvalues'] = varvalues
        self.__dict__['name'] = name
    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        if name in self.varvalues:
            return self.varvalues[name]
        else:
            return object.__getattr__(self, name)
    def __setattr__(self, name, value):
        if name != 'name' and name != 'varvalues':
            self.__dict__['varvalues'][name]=value
        else:
            object.__setattr__(self, name, value)
    def __str__(self):
        res = "%s:\n"% self.name
        for x, value in self.varvalues.items():
            res += "   %s = %s\n" % (x, str(float(value)))
        return res

def state(**varvalues):
    return StateArray([], dict(**varvalues), '?')

class Transformation(object):
    def __init__(self, rate = 0.0):
        self.rate = rate
        self.name = '?'
    def __str__(self):
        return "%s:\n  rate = %s\n" % (self.name, str(self.rate))

class Constant(float):
    def __init__(self, value = 0.0):
        float.__init__(self,value)
        #self.value = value
        self.name = '?'
    #~ def __str__(self):
        #~ return "\n%s = %f" % (self.name, self) #.value)

class PairConstants(object):
    def __init__(self, min = 0.0, max = 1.0):
        self.min = min
        self.max = max
        self.name = '?'
    def __str__(self):
        return "%s = ? (min=%f, max=%f)\n" % (self.name, self.min, self.max)

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
        self.__dict__['_Model__states']            = []
        self.__dict__['title']                     = title
        self.__dict__['_Model__m_Parameters']      = None
    
    def __setattr__(self, name, value):
        if isinstance(value, float) or isinstance(value, int):
            value = Constant(float(value))
        if (isinstance(value, tuple) or isinstance(value, list)) and len(value)==2:
            value = PairConstants(float(value[0]), float(value[1]))
        assoc = ((Reaction,       '_Model__reactions'),
                 (Constant,       '_Model__parameters'),
                 (Transformation, '_Model__transf'),
                 (StateArray,     '_Model__states'),
                 (PairConstants,  '_Model__unknownparameters'))
        for t, dictname in assoc:
            c = findWithNameIndex(name, self.__dict__[dictname])
            if c > -1:
                if isinstance(value, t):
                    value.name = name
                    self.__dict__[dictname][c] = value
                    self.__refreshVars()
                    return
                else:
                    raise BadTypeComponent, name+ ' can not be assigned to ' + type(value).__name__
        for t, dictname in assoc:
            if isinstance(value, t):
                value.name = name
                self.__dict__[dictname].append(value)
                self.__refreshVars()
                return
        object.__setattr__(self, name, value)

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        c = findWithName(name, self.__parameters)
        if c :
            return c #.value
        c = findWithName(name, self.__reactions)
        if c :
            return c
        c = findWithName(name, self.__variables)
        if c :
            return c.name
        c = findWithName(name, self.__transf)
        if c :
            return c
        c = findWithName(name, self.__unknownparameters)
        if c :
            return c
        i = findWithNameIndex(name, self.__states)
        if i >-1 :
            #replace StateArray with new values as a list
            prov = self.__states[i]
            newlist = [prov.varvalues.get(var.name,0.0) for var in self.__variables]
            c = StateArray(newlist, prov.varvalues, prov.name)
            self.__states[i] = c
            return c
        raise AttributeError, name + ' is not defined for this model'
    
    def varnames(self):
        return [i.name for i in self.__variables]
    
    def checkRates(self):
        for collection in (self.__reactions, self.__transf):
            for v in collection:
                resstring, value = test_with_everything(v.rate, 
                                                        self.__parameters, 
                                                        self.__unknownparameters, 
                                                        self.__variables)
                if resstring != "":
                    return False, resstring + '\nin rate of %s\n(%s)' % (v.name, v.rate)
        return True, 'OK'

    def rateCalcString(self, rateString):
        # replace varnames
        for i,v in enumerate(self.__variables):
            rateString = re.sub(r"\b"+ v.name+r"\b", "variables[%d]"%i, rateString)
        # replace unknown parameters
        for i,u in enumerate(self.__unknownparameters):
            rateString = re.sub(r"\b"+ u.name+r"\b", "m_Parameters[%d]"%i, rateString)
        # replace parameters
        for p in self.__parameters:
            rateString = re.sub(r"\b"+ p.name + r"\b", "%g"% p, rateString)  #.value, rateString)
        return rateString

    def __str__(self):
        check, msg = self.checkRates()
        if not check:
            res ="Problem in rates:\n"
            res += msg
            return res
        res = "%s\n"% self.title
        res += "\nVariables: %s\n" % " ".join([i.name for i in self.__variables])
        if len(self.__extvariables) > 0:
            res += "External variables: %s\n" % " ".join([i.name for i in self.__extvariables])
        for collection in (self.__reactions, self.__transf, self.__unknownparameters, self.__states):
            for i in collection:
                res += str(i)
        for p in self.__parameters:
            res += p.name +' = '+ str(p) + '\n'
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
        #print "DEBUG, vars = ", [x.name for x in self.__variables]


    def genStoichiometryMatrix(self):
        check, msg = self.checkRates()
        if not check:
            raise BadRateError, msg

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

    def rates_func(self):
        """Generate function to compute rate vector for this model.
        
           Function has signature f(variables, t)"""

        check, msg = self.checkRates()
        if not check:
            raise BadRateError, msg

        #compile rate laws
        ratebytecode = [compile(self.rateCalcString(v.rate), 'bof.log','eval') for v in self.__reactions]
        # create array to hold v's
        v = empty(len(self.reactions))
            
        def f(variables, t):
            m_Parameters = self.__m_Parameters
            for i,r in enumerate(ratebytecode):
                v[i] = eval(r)
            return v
        return f

    def transf_func(self):
        """Generate function to compute transformations for this model.
        
           Function has signature f(variables, t)"""

        check, msg = self.checkRates()
        if not check:
            raise BadRateError, msg

        #compile rate laws
        transfbytecode = [compile(self.rateCalcString(v.rate), 'bof.log','eval') for v in self.__transf]
        # create array to hold v's
        m_Transformations = empty(len(self.transf))
            
        def f(variables, t):
            m_Parameters = self.__m_Parameters
            for i,r in enumerate(transfbytecode):
                m_Transformations[i] = eval(r)
            return m_Transformations
        return f

    def set_unknown(self, unknownparameters):
        self.__m_Parameters = unknownparameters
    
    def getdXdt(self):
        """Generate function to compute rhs of SODE for this model.
        
           Function has signature f(variables, t)
           This is compatible with scipy.integrate.odeint"""

        check, msg = self.checkRates()
        if not check:
            raise BadRateError, msg
        #compile rate laws
        ratebytecode = [compile(self.rateCalcString(v.rate), 'bof.log','eval') for v in self.__reactions]
        # compute stoichiometry matrix, scale and transpose
        N  = self.genStoichiometryMatrix()
        NT = N.transpose()

        # create array to hold v's
        v = empty(len(self.reactions))
        
        def f(variables, t):
            m_Parameters = self.__m_Parameters
            for i,r in enumerate(ratebytecode):
                v[i] = eval(r)
            return dot(v,NT)
        return f

    def scaled_dXdt(self, scale = 1.0):
        """Generate function to compute rhs of SODE for this model.
        
           Function has signature f(variables, t)
           This is compatible with scipy.integrate.odeint"""

        check, msg = self.checkRates()
        if not check:
            print "vars = "
            print [x.name for x in self.variables]
            raise BadRateError, msg
        #compile rate laws
        ratebytecode = [compile(self.rateCalcString(v.rate), 'bof.log','eval') for v in self.__reactions]
        # compute stoichiometry matrix, scale and transpose
        N  = self.genStoichiometryMatrix()
        N *= scale
        NT = N.transpose()

        # create array to hold v's
        v = empty(len(self.reactions))
        #self.__m_Parameters = unknownparameters

        def f(variables, t):
            m_Parameters = self.__m_Parameters
            for i,r in enumerate(ratebytecode):
                v[i] = eval(r)
            return dot(v,NT)
        return f

    def dXdt_with(self, unknownparameters, scale = 1.0):
        """Generate function to compute rhs of SODE for this model.
        
           Function has signature f(variables, t)
           This is compatible with scipy.integrate.odeint"""

        check, msg = self.checkRates()
        if not check:
            raise BadRateError, msg
        #compile rate laws
        ratebytecode = [compile(self.rateCalcString(v.rate), 'bof.log','eval') for v in self.__reactions]
        # compute stoichiometry matrix, scale and transpose
        N  = self.genStoichiometryMatrix()
        N *= scale
        NT = N.transpose()

        # create array to hold v's
        v = empty(len(self.reactions))

        def f(variables, t):
            m_Parameters = unknownparameters
            for i,r in enumerate(ratebytecode):
                v[i] = eval(r)
            return dot(v,NT)
        return f


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
    dXdt         = property(getdXdt)

class BadStoichError(Exception):
    """Used to flag a wrong stoichiometry expression"""

class BadRateError(Exception):
    """Used to flag a wrong stoichiometry expression"""

class BadTypeComponent(Exception):
    """Used to flag an assignement of a model component to a wrong type object"""


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
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.t1 = transf("A*4 + C")
    m.t2 = transf("sqrt(2*A)")
    m.B  = 2.2
    m.myconstant = 2 * m.B / 1.1 # should be 4.0
    m.V3 = [0.1, 0.9]
    m.Km3 = 4
    m.init = state(A = 1.0, C = 1.0)
    print '********** Testing model construction and printing **********'
    print m
    
    #~ print m.B * m.myconstant, m.v3.rate
    #~ print m.v1
    print '********** Testing component retrieval *********************'
    print 'm.K3 :',m.Km3
    print 'm.K3.name :',m.Km3.name, '(a float with a name attr)'
    print m.init
    print 'm.init[1] :',m.init[1]

    print '********** Testing component reassignement *****************'
    print 'm.myconstant :',m.myconstant
    print len(m.parameters), 'parameters total'
    print 'making m.myconstant = 5.0'
    m.myconstant = 5.0
    print 'm.myconstant :',m.myconstant
    print len(m.parameters), 'parameters total'

    print 'making m.myconstant = [5.0, 10.0]'
    try:
        m.myconstant = [5.0, 10.0]
    except BadTypeComponent:
        print 'Failed! BadTypeComponent was caught.'
    print 'm.myconstant :',m.myconstant, '(still!)'
    print len(m.parameters), 'parameters total'

    print m.V3
    print len(m.unknown), 'unknown parameters total'
    print 'making m.V3 = [0.1, 0.2]'
    m.V3 = [0.1, 0.2]
    print m.V3
    print len(m.unknown), 'unknown parameters total'
    print 

    print '********** Testing stoichiometry matrix ********************'
    print 'Stoichiometry matrix:'
    N = m.genStoichiometryMatrix()
    print '  ', '  '.join([v.name for v in m.reactions])
    for i,x in enumerate(m.variables):
        print x.name, N[i, :]
    print
    print '********** Testing rateCalcString ******'
    print 'calcstring for v3:\n', m.rateCalcString(m.v3.rate)
    print
    print '********** Testing rate and dXdt generating functions ******'
    print 'Operating point:'
    varvalues = [1.0, 1.0]
    pars      = [1.0]
    m.set_unknown(pars)
    print 'variables =', dict((v.name, value) for v,value in zip(m.variables, varvalues))
    print 'unknown   =', dict((v.name, value) for v,value in zip(m.unknown, pars))
    print 'parameters=', dict((p.name, p)     for p in m.parameters)
 
    print '---- rates using Model.rates_func() -------------------------'
    vratesfunc = m.rates_func()
    vrates = vratesfunc(varvalues,0)
    for v,r in zip(m.reactions, vrates):
        print "%s = %-20s = %s" % (v.name, v.rate, r)

    print '---- transformations using Model.transf_func() --------------'
    tratesfunc = m.transf_func()
    trates = tratesfunc(varvalues,0)
    for v,r in zip(m.transf, trates):
        print "%s = %-20s = %s" % (v.name, v.rate, r)

    print '---- dXdt using Model.dXdt() --------------------------------'
    #f = m.dXdt()
    m.set_unknown(pars)
    dXdt = m.dXdt(varvalues,0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)

    print '---- dXdt using Model.dXdt_with(pars) ------------------------'
    f = m.dXdt_with(pars)
    dXdt   = f(varvalues,0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)

    print '---- dXdt using Model.dXdt() with a state argument (m.init) --'
    print m.init
    f = m.dXdt
    m.set_unknown(pars)
    print 'dXdt = f(m.init,0)'
    dXdt = f(m.init,0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)
    print '---- same, changing state argument ---------------------------'
    m.init.A = 2.0
    print 'after m.init.A = 2.0'
    print m.init
    f = m.dXdt
    m.set_unknown(pars)
    print 'dXdt = f(m.init,0)'
    dXdt = f(m.init,0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)

    print '---------------- EXAMPLE 4 ------------------'
    m4 = Model("Rossler")
    m4.v1 = react(" -> X1", rate = "X2 - X3")
    m4.v2 = react(" -> X2", rate = "0.36 * X2 - X1")
    m4.v3 = react(" -> X3", rate = "X1 *X3 - 22.5 * X3 - 49.6 * X1 + 1117.8")
    m4.init = state(X1 = 19.0, X2 = 47, X3 = 50)

    print m4

if __name__ == "__main__":
    main()
 