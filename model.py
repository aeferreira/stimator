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

def test_with_everything(valueexpr, parameters, varlist): 
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
    #part 2: permissive, with dummy values (1.0) for vars
    vardict = {}
    for i in varlist:
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

stoichiom  = re.compile(stoichiompattern,    re.IGNORECASE)
chemcomplex = re.compile(chemcomplexpattern, re.IGNORECASE)

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
        self.irreversible = irreversible
        self.name = '?'
    def __str__(self):
        return "%s:\n  reagents: %s\n  products: %s\n  rate    = %s\n" % (self.name, str(self.reagents), str(self.products), str(self.rate)) 

class StateArray(object):
    def __init__(self, varvalues, name):
        self.__dict__['name'] = name
        for k,v in varvalues.items():
            varvalues[k] = constValue(value = v, name = k)
        self.__dict__['varvalues'] = varvalues
    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        if name in self.varvalues:
            return self.varvalues[name]
        else:
            raise AttributeError, name + ' is not a member of state'+ self.__dict__['name']
    def __setattr__(self, name, value):
        if name != 'name' and name != 'varvalues':
            value = constValue(value = value, name = name, into = self.__dict__['varvalues'].get(name, None))
            self.__dict__['varvalues'][name]=value
        else:
            object.__setattr__(self, name, value)
    def __str__(self):
        return '(%s)' % ", ".join(['%s = %s'% (x,str(float(value))) for (x,value) in  self.varvalues.items()])
    def __iter__(self):
        return iter(self.varvalues.items())

def state(**varvalues):
    return StateArray(varvalues, '?')

class Transformation(object):
    def __init__(self, rate = 0.0):
        self.rate = rate
        self.name = '?'
    def __str__(self):
        return "%s:\n  rate = %s\n" % (self.name, str(self.rate))

class ConstValue(float):
    def __init__(self, value):
        float.__init__(self,value)
        self.name = '?'
        self.bounds = None
    def uncertainty(self, *pars):
        if len(pars) == 0 or pars[0]==None:
            self.bounds = None
            return
        if len(pars) != 2:
            return #TODO raise exception
        self.bounds = Bounds(float(pars[0]), float(pars[1]))
        self.bounds.name = self.name
    def pprint(self):
        res = float.__str__(self)
        if self.bounds:
            res+= " ? (min = %f, max=%f)" % (self.bounds.min, self.bounds.max)
        return res

            
def constValue(value = None, name = None, into = None):
    if isinstance(value, float) or isinstance(value, int):
        v = float(value)
        res = ConstValue(v)
        if into:
            res.name = into.name
            res.bounds = into.bounds
    elif (isinstance(value, tuple) or isinstance(value, list)) and len(value)==2:
        bounds = Bounds(float(value[0]), float(value[1]))
        v = (bounds.min + bounds.max)/2.0
        if into:
            res = ConstValue(into)
            res.name = into.name
        else:
            res = ContsValue(v)
        res.bounds = bounds
    else:
        raise TypeError, value + ' is not a float or pair of floats'
    if name:
        res.name = name
    if res.bounds:
        res.bounds.name = res.name        
    return res
        

class Bounds(object):
    def __init__(self, min = 0.0, max = 1.0):
        self.min = min
        self.max = max
        self.name = '?'
    def __str__(self):
        return "(min=%f, max=%f)" % (self.min, self.max)

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
        self.__dict__['_Model__transf']            = []
        self.__dict__['_Model__states']            = []
        self.__dict__['title']                     = title
        self.__dict__['_Model__m_Parameters']      = None
    
    def __setattr__(self, name, value):
        if isinstance(value, float) or isinstance(value, int):
            value = ConstValue(float(value))
        if (isinstance(value, tuple) or isinstance(value, list)) and len(value)==2:
            value = Bounds(float(value[0]), float(value[1]))
        assoc = ((Reaction,       '_Model__reactions'),
                 (ConstValue,     '_Model__parameters'),
                 (Transformation, '_Model__transf'),
                 (StateArray,     '_Model__states'))
        for t, dictname in assoc:
            c = findWithNameIndex(name, self.__dict__[dictname])
            if c > -1:
                if isinstance(value, t) or (isinstance(value, Bounds) and dictname == '_Model__parameters'):
                    value.name = name
                    if isinstance(value, Bounds):
                        self.__dict__[dictname][c].bounds = value
                    else:
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
        if isinstance (value, Bounds):
            value.name = name
            newvalue = ConstValue((float(value.min)+float(value.max))/2.0)
            newvalue.name = name
            newvalue.bounds = value
            self.__dict__['_Model__parameters'].append(newvalue)
            self.__refreshVars()
            return
        object.__setattr__(self, name, value)

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        c = findWithName(name, self.__parameters)
        if c :
            return c
        c = findWithName(name, self.__reactions)
        if c :
            return c
        c = findWithName(name, self.__variables)
        if c :
            return c.name
        c = findWithName(name, self.__transf)
        if c :
            return c
        c = findWithName(name, self.__states)
        if c :
            return c
        raise AttributeError, name + ' is not defined for this model'
    
    def vectorize(self, state):
        if isinstance(state, str) or isinstance(state, unicode):
            if not hasattr(self, state):
                raise AttributeError, state + ' is not defined for this model'
            state = getattr(self, state)
        newlist = [state.varvalues.get(var.name,0.0) for var in self.__variables]
        return array(newlist)
        
    
    def varnames(self):
        return [i.name for i in self.__variables]
    
    def checkRates(self):
        for collection in (self.__reactions, self.__transf):
            for v in collection:
                resstring, value = test_with_everything(v.rate, 
                                                        self.__parameters, 
                                                        self.__variables)
                if resstring != "":
                    return False, resstring + '\nin rate of %s\n(%s)' % (v.name, v.rate)
        return True, 'OK'

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
        for collection in (self.__reactions, self.__transf):
            for i in collection:
                res += str(i)
        for p in self.__states:
            res += p.name +': '+ str(p) + '\n'
        for p in self.__parameters:
            res += p.name +' = '+ str(p) + '\n'
        for u in self.uncertain:
            res += u.name + ' = ? (' + str(u.min) + ', ' + str(u.max) + ')\n'
        return res
    
    def clone(self):
        m = Model(self.title)
        for r in self.reactions:
            setattr(m, r.name, Reaction(r.reagents, r.products, r.rate, r.irreversible))
        for p in self.parameters:
            setattr(m, p.name, ConstValue(p))
        for t in self.transf:
            setattr(m, t.name, Transformation(t.rate))
        for s in self.__states:
            newdict= {}
            for i in s:
                newdict[i[0]]=i[1]
            setattr(m, s.name, StateArray(newdict, s.name))
        #handle uncertainties
        for u in self.uncertain:
            loc = u.name.split('.')
            if len(loc) >1:
                s = getattr(m,loc[0])
                var = getattr(s,loc[1])
                var.uncertainty(u.min, u.max)
            else:
                getattr(m, loc[0]).uncertainty(u.min, u.max)
        return m
    
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

    def rateCalcString(self, rateString, with_uncertain = False):
        # replace varnames
        for i,v in enumerate(self.variables):
            rateString = re.sub(r"\b"+ v.name+r"\b", "variables[%d]"%i, rateString)
        # replace uncertain parameters
        if with_uncertain:
            for i,u in enumerate(self.uncertain):
                rateString = re.sub(r"\b"+ u.name+r"\b", "m_Parameters[%d]"%i, rateString)
        # replace parameters
        for p in self.parameters:
            if p.bounds and with_uncertain:
                continue
            rateString = re.sub(r"\b"+ p.name + r"\b", "%g"% p, rateString) 
        return rateString

    def rates_func(self, with_uncertain = False):
        """Generate function to compute rate vector for this model.
        
           Function has signature f(variables, t)"""

        check, msg = self.checkRates()
        if not check:
            raise BadRateError, msg

        #compile rate laws
        ratebytecode = [compile(self.rateCalcString(v.rate, with_uncertain=with_uncertain), 'bof.log','eval') for v in self.__reactions]
        # create array to hold v's
        v = empty(len(self.reactions))
            
        def f(variables, t):
            for i,r in enumerate(ratebytecode):
                v[i] = eval(r)
            return v
        def f2(variables, t):
            m_Parameters = self.__m_Parameters
            for i,r in enumerate(ratebytecode):
                v[i] = eval(r)
            return v

        if with_uncertain:
            return f2
        else:
            return f

    def transf_func(self, with_uncertain = False):
        """Generate function to compute transformations for this model.
        
           Function has signature f(variables, t)"""

        check, msg = self.checkRates()
        if not check:
            raise BadRateError, msg

        #compile rate laws
        transfbytecode = [compile(self.rateCalcString(v.rate, with_uncertain=with_uncertain), 'bof.log','eval') for v in self.__transf]
        # create array to hold v's
        m_Transformations = empty(len(self.transf))
            
        def f(variables, t):
            for i,r in enumerate(transfbytecode):
                m_Transformations[i] = eval(r)
            return m_Transformations
        def f2(variables, t):
            m_Parameters = self.__m_Parameters
            for i,r in enumerate(transfbytecode):
                m_Transformations[i] = eval(r)
            return m_Transformations

        if with_uncertain:
            return f2
        else:
            return f

    def set_uncertain(self, uncertainparameters):
        self.__m_Parameters = uncertainparameters
    
    def getdXdt(self, with_uncertain = False):
        """Generate function to compute rhs of SODE for this model.
        
           Function has signature f(variables, t)
           This is compatible with scipy.integrate.odeint"""

        check, msg = self.checkRates()
        if not check:
            raise BadRateError, msg
        #compile rate laws
        ratebytecode = [compile(self.rateCalcString(v.rate, with_uncertain = with_uncertain), 'bof.log','eval') for v in self.__reactions]
        # compute stoichiometry matrix, scale and transpose
        N  = self.genStoichiometryMatrix()
        NT = N.transpose()

        # create array to hold v's
        v = empty(len(self.reactions))
        
        def f(variables, t):
            for i,r in enumerate(ratebytecode):
                v[i] = eval(r)
            return dot(v,NT)

        def f2(variables, t):
            m_Parameters = self.__m_Parameters
            for i,r in enumerate(ratebytecode):
                v[i] = eval(r)
            return dot(v,NT)
        
        if with_uncertain:
            return f2
        else:
            return f

    def scaled_dXdt(self, scale = 1.0, with_uncertain = False):
        """Generate function to compute rhs of SODE for this model.
        
           Function has signature f(variables, t)
           This is compatible with scipy.integrate.odeint"""

        check, msg = self.checkRates()
        if not check:
            print "vars = "
            print [x.name for x in self.variables]
            raise BadRateError, msg
        #compile rate laws
        ratebytecode = [compile(self.rateCalcString(v.rate, with_uncertain = with_uncertain), 'bof.log','eval') for v in self.__reactions]
        # compute stoichiometry matrix, scale and transpose
        N  = self.genStoichiometryMatrix()
        N *= scale
        NT = N.transpose()

        # create array to hold v's
        v = empty(len(self.reactions))

        def f(variables, t):
            for i,r in enumerate(ratebytecode):
                v[i] = eval(r)
            return dot(v,NT)
        def f2(variables, t):
            m_Parameters = self.__m_Parameters
            for i,r in enumerate(ratebytecode):
                v[i] = eval(r)
            return dot(v,NT)
        if with_uncertain:
            return f2
        else:
            return f

    def dXdt_with(self, uncertainparameters, scale = 1.0):
        """Generate function to compute rhs of SODE for this model.
        
           Function has signature f(variables, t)
           This is compatible with scipy.integrate.odeint"""

        check, msg = self.checkRates()
        if not check:
            raise BadRateError, msg
        #compile rate laws
        ratebytecode = [compile(self.rateCalcString(v.rate, with_uncertain = True), 'bof.log','eval') for v in self.__reactions]
        # compute stoichiometry matrix, scale and transpose
        N  = self.genStoichiometryMatrix()
        N *= scale
        NT = N.transpose()

        # create array to hold v's
        v = empty(len(self.reactions))

        def f(variables, t):
            m_Parameters = uncertainparameters
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
    def __getTransformations(self):
        return self.__transf
    
    def __getUncertainValues(self):
        return list(getUncertainties(self)) #generator of uncertain values

    variables    = property(__getVariables)
    extvariables = property(__getExtVariables)
    reactions    = property(__getReactions)
    parameters   = property(__getParameters)
    uncertain    = property(__getUncertainValues)
    transf       = property(__getTransformations)
    dXdt         = property(getdXdt)

def getUncertainties(model):
    for p in model.parameters:
        if p.bounds:
            yield p.bounds
    for s in model._Model__states:
        for name, value in s:
            if value.bounds:
                ret = Bounds(value.bounds.min, value.bounds.max)
                ret.name = s.name + '.' + value.name
                yield ret

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

def test():
    
    m = Model('My first model')
    m.v1 = react("A+B -> C"  , 3)
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.t1 = transf("A*4 + C")
    m.t2 = transf("sqrt(2*A)")
    m.B  = 2.2
    m.myconstant = 2 * m.B / 1.1 # should be 4.0
    m.V3 = 0.5
    m.V3 = [0.1, 1.0]
    m.Km3 = 4
    m.init = state(A = 1.0, C = 1)
    m.afterwards = state(A = 1.0, C = 2)
    m.afterwards.C.uncertainty(1,3)
    
    print '********** Testing model construction and printing **********'
    print '------- result of model construction:\n'
    print m
    m2 = m.clone()
    print
    print '------- result of CLONING the model:\n'
    print m2
    
    print '********** Testing component retrieval *********************'
    print 'm.K3 :',m.Km3
    print 'm.K3.name :',m.Km3.name, '(a float with a name attr)'
    print 'm.init:',m.init
    print 'm.init.A :',m.init.A
    print 'iterating m.init'
    for name, x in m.init:
        print '\t', name, '=', x
    print

    print '********** Testing component reassignment *****************'
    print 'm.myconstant :',m.myconstant
    print len(m.parameters), 'parameters total'
    print 'making m.myconstant = 5.0'
    m.myconstant = 5.0
    print 'm.myconstant :',m.myconstant
    print len(m.parameters), 'parameters total'

    print 'making m.myconstant = react("A+B -> C"  , 3)'
    try:
        m.myconstant = react("A+B -> C"  , 3)
    except BadTypeComponent:
        print 'Failed! BadTypeComponent was caught.'
    print 'm.myconstant :',m.myconstant, '(still!)'
    print len(m.parameters), 'parameters total'
    print
    print 'm.V3 :', m.V3
    print 'm.V3.bounds:' , m.V3.bounds
    print 'iterating m.uncertain'
    for x in m.uncertain:
        print '\t', x.name, 'in (', x.min, ',', x.max, ')'
    print len(m.uncertain), 'uncertain parameters total'
    print 'making m.V3 = [0.1, 0.2]'
    m.V3 = [0.1, 0.2]
    print 'm.V3 :', m.V3
    print 'm.V3.bounds:' ,m.V3.bounds
    print len(m.uncertain), 'uncertain parameters total'
    print 'making m.V4 = [0.1, 0.6]'
    m.V4 = [0.1, 0.6]
    print 'm.V4 :', m.V4
    print 'm.V4.bounds:' ,m.V4.bounds
    print len(m.uncertain), 'uncertain parameters total'
    print 'iterating m.uncertain'
    for x in m.uncertain:
        print '\t', x.name, 'in (', x.min, ',', x.max, ')'
    print 'making m.init.A = 5.0'
    m.init.A = 5.0
    print 'iterating m.init'
    for name, x in m.init:
        print '\t', name, '=', x.pprint()
    print 'flagging init.A as uncertain with   m.init.A = (0.5, 2.5)'
    m.init.A = (0.5, 2.5)
    print 'iterating m.init'
    for name, x in m.init:
        print '\t', name, '=', x.pprint()
    print 'calling    m.init.A.uncertainy(0.5,3.0)'
    m.init.A.uncertainty(0.5,3.0)
    print 'iterating m.init'
    for name, x in m.init:
        print '\t', name, '=', x.pprint()
    print 'calling    m.init.A.uncertainy(None)'
    m.init.A.uncertainty(None)
    print 'iterating m.init'
    for name, x in m.init:
        print '\t', name, '=', x.pprint()
    print 'making m.init.A back to 1.0'
    m.init.A = 1.0
    print 'iterating m.init'
    for name, x in m.init:
        print '\t', name, '=', x.pprint()
    print 

    print '********** Testing stoichiometry matrix ********************'
    print 'Stoichiometry matrix:'
    N = m.genStoichiometryMatrix()
    print '  ', '  '.join([v.name for v in m.reactions])
    for i,x in enumerate(m.variables):
        print x.name, N[i, :]
    print
    print '********** Testing rateCalcString **************************'
    print 'calcstring for v3:\n', m.rateCalcString(m.v3.rate)
    print
    print 'calcstring for v3 with uncertain parameters:\n', m.rateCalcString(m.v3.rate, True)
    print

    print '********** Testing rate and dXdt generating functions ******'
    print 'Operating point:'
    varvalues = [1.0, 1.0]
    pars      = [1.0]

    print 'variables  =', dict((v.name, value) for v,value in zip(m.variables, varvalues))
    print 'parameters =', dict((p.name, p)     for p in m.parameters)
 
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
    #f = m.getdXdt()
    dXdt = m.dXdt(varvalues,0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)

    print '---- dXdt using Model.dXdt() setting uncertain parameters ---'
    print 'f = m.getdXdt(with_uncertain = True)'
    f = m.getdXdt(with_uncertain = True)
    print 'setting uncertain as', dict((v.name, value) for v,value in zip(m.uncertain, pars))
    print 'm.set_uncertain(pars)'
    m.set_uncertain(pars)
    dXdt = f(varvalues,0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)

    print '---- dXdt using Model.dXdt_with(pars) ------------------------'
    print 'f = m.dXdt_with(pars)'
    f = m.dXdt_with(pars)
    dXdt   = f(varvalues,0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)

    print '---- dXdt using Model.dXdt() with a state argument (m.init) --'
    print 'm.init:', m.init
    print 'making m.V3 = 1.0'
    m.V3 = 1.0
    print 'm.V3 :', m.V3
    print
    print 'f = m.dXdt'
    f = m.dXdt
    print 'dXdt = f(m.vectorize("init"),0)'
    dXdt = f(m.vectorize("init"),0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)
    print '---- same, changing state argument ---------------------------'
    m.init.A = 2.0
    print 'after m.init.A = 2.0'
    print 'm.init:', m.init
    print
    print 'f = m.dXdt'
    f = m.dXdt
    print 'dXdt = f(m.vectorize("init"),0)'
    dXdt = f(m.vectorize("init"),0)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)

if __name__ == "__main__":
    test()
 