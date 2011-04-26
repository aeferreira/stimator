#!/usr/bin/env python
# -*- coding: latin1-*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator Model related classes
# Copyright António Ferreira 2006-2011
#----------------------------------------------------------------------------
import re
import math
from kinetics import *
from numpy import *
import pprint

#----------------------------------------------------------------------------
#         Functions to check the validity of math expressions
#----------------------------------------------------------------------------
__globs = {}
__haskinetics = {}
for k, v in globals().items():
    if hasattr(v,"is_rate"):
        __haskinetics[k] = v
__globs.update(__haskinetics)
## pprint.pprint(__globs)
__globs.update(vars(math))
## pprint.pprint(__globs)

def register_kin_func(f):
    f.is_rate = True
    __globs[f.__name__] = f
    globals()[f.__name__] = f

def _test_with_everything(valueexpr, model): 
    locs = {}
    for p in parameters(model):
        locs[p.name] = p #.value
    
    #part 1: nonpermissive, except for NameError
    try :
       value = float(eval(valueexpr, __globs, locs))
    except NameError:
       pass
    except Exception, e:
       return ("%s : %s"%(str(e.__class__.__name__),str(e)), 0.0)
    #part 2: permissive, with dummy values (1.0) for vars
    vardict = {}
    for i in variables(model):
        vardict[i.name]=1.0
    vardict['t'] = 1.0
    locs.update(vardict)
    try :
       value = float(eval(valueexpr, __globs, locs))
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

identifier = re.compile(r"[_a-z]\w*", re.IGNORECASE)

def identifiersInExpr(_expr):
    iterator = identifier.finditer(_expr)
    return [_expr[m.span()[0]:m.span()[1]] for m in iterator]

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
    res = str(float(k))
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
    def __init__(self, reagents, products, rate, irreversible = False):
        self.reagents = reagents
        self.products = products
        self.rate = rate.strip()
        self.irreversible = irreversible
        self.name = '?'
    def __str__(self):
        return "%s:\n  reagents: %s\n  products: %s\n  rate    = %s\n" % (self.name, str(self.reagents), str(self.products), str(self.rate)) 

def react(stoichiometry, rate = 0.0):
    res = processStoich(stoichiometry)
    if not res:
        raise BadStoichError( "Bad stoichiometry definition:\n"+ stoichiometry)
    if isinstance(rate, float) or isinstance(rate, int):
        rate = massActionStr(rate, res[0])
    return Reaction(res[0], res[1], rate, res[2])


class Variable_dXdt(object):
    def __init__(self, rate = 0.0):
        self.rate = rate
        self.name = '?'
    def __str__(self):
        return "%s:\n  rate = %s\n" % (self.name, str(self.rate))

def variable(rate = 0.0):
    if isinstance(rate, float) or isinstance(rate, int):
        rate = str(rate)
    return Variable_dXdt(rate)

class Transformation(object):
    def __init__(self, rate):
        self.rate = rate.strip()
        self.name = '?'
    def __str__(self):
        return "%s:\n  rate = %s\n" % (self.name, str(self.rate))

def transf(rate = 0.0):
    if isinstance(rate, float) or isinstance(rate, int):
        rate = str(float(rate))
    return Transformation(rate)


class ConstValue(float):
    def __new__(cls, value):
        return float.__new__(cls,value)
##     def __init__(self):
##         self.name = '?'
##         self.bounds = None
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
        res.name = '?'
        res.bounds = None
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
            res.name = '?'
        res.bounds = bounds
    else:
        raise TypeError( value + ' is not a float or pair of floats')
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
            raise AttributeError( name + ' is not a member of state'+ self.__dict__['name'])
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


def isPairOfNums(value):
    if (isinstance(value, tuple) or isinstance(value, list)) and len(value)==2:
        for pos in (0,1):
            typeOK = False
            for numtype in (float,int,long):
                if isinstance(value[pos], numtype):
                    typeOK = True
                    break
            if not typeOK:
                return False
        return True
    return False

def ConvertPair2Reaction(value):
    if (isinstance(value, tuple) or isinstance(value, list)) and len(value)==2:
        if isinstance(value[0], str):
            good2nd = False
            for numtype in (str, float,int,long):
                if isinstance(value[1], numtype):
                    good2nd = True
                    break
            if not good2nd:
                return False
            res = processStoich(value[0])
            if not res:
                raise BadStoichError( "Bad stoichiometry definition:\n"+ value[0])
            rate = value[1]
            for numtype in (float,int,long):
                if isinstance(rate, numtype):
                    rate = massActionStr(rate, res[0])
            return Reaction(res[0], res[1], rate, res[2])
    return False

class Model(object):
    def __init__(self, title = ""):
        self.__dict__['_Model__reactions']         = []
        self.__dict__['_Model__variables']         = []
        self.__dict__['_Model__extvariables']      = []
        self.__dict__['_Model__parameters']        = []
        self.__dict__['_Model__transf']            = []
        self.__dict__['_Model__states']            = []
        self.__dict__['_Model__metadata']          = {}
        #self.__dict__['title']                     = title
        self.__dict__['_Model__m_Parameters']      = None
        self.setData('title', title)
    
    def __setattr__(self, name, value):
        for numtype in (float,int,long):
            if isinstance(value, numtype):
                value = constValue(float(value))
        if isPairOfNums(value):
            value = Bounds(float(value[0]), float(value[1]))
        if isinstance(value,Variable_dXdt):
            react_name = 'd_%s_dt'% name
            stoich = ' -> %s'% name
            name = react_name # hope this works...
            value = react(stoich, value.rate)
        r = ConvertPair2Reaction(value)
        if r:
            value = r
        assoc = ((Reaction,       '_Model__reactions'),
                 (ConstValue,     '_Model__parameters'),
                 (Transformation, '_Model__transf'),
                 (StateArray,     '_Model__states'))
        # find existing model object
        for t, listname in assoc:
            c = findWithNameIndex(name, self.__dict__[listname])
            if c > -1:
                if isinstance(value, t) or (isinstance(value, Bounds) and listname == '_Model__parameters'):
                    value.name = name
                    if isinstance(value, Bounds):
                        self.__dict__[listname][c].bounds = value
                    else:
                        self.__dict__[listname][c] = value
                    self.__refreshVars()
                    return
                else:
                    raise BadTypeComponent( name + ' can not be assigned to ' + type(value).__name__)
        # append new object to proper list
        for t, listname in assoc:
            if isinstance(value, t):
                value.name = name
                self.__dict__[listname].append(value)
                self.__refreshVars()
                return
        if isinstance (value, Bounds):
            value.name = name
            newvalue = constValue((float(value.min)+float(value.max))/2.0, name)
            newvalue.bounds = value
            self.__dict__['_Model__parameters'].append(newvalue)
            self.__refreshVars()
            return
        object.__setattr__(self, name, value)

    def __getattr__(self, name):
        if name == '__m_Parameters':
            return self.__dict__['_Model__m_Parameters']
        if name in self.__dict__:
            return self.__dict__[name]
        c = findWithNameIndex(name, self.__parameters)
        if c >=0:
            return self.__parameters[c]
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
        raise AttributeError( name + ' is not defined for this model')
    
    def setData(self, name, value):
        self.__metadata[name] = value
    
    def getData(self, name):
        if not name in self.__metadata:
            return None
        return self.__metadata[name]
        

    def __findComponent(self, name):
        c = findWithNameIndex(name, self.__parameters)
        if c>=0 :
            return c, 'parameters'
        c = findWithNameIndex(name, self.__reactions)
        if c>=0 :
            return c, 'reactions'
        c = findWithNameIndex(name, self.__variables)
        if c>=0 :
            return c, 'variables'
        c = findWithNameIndex(name, self.__transf)
        if c>=0 :
            return c, 'transf'
        raise AttributeError( name + ' is not a component in this model')

    def checkRates(self):
        for collection in (self.__reactions, self.__transf):
            for v in collection:
                resstring, value = _test_with_everything(v.rate,self)
                if resstring != "":
                    return False, '%s\nin rate of %s: %s' % (resstring, v.name, v.rate)
        return True, 'OK'

    def __str__(self):
        check, msg = self.checkRates()
        if not check:
            raise BadRateError(msg)
        res = "%s\n"% self.getData('title')
        #~ res = "%s\n"% self.title
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
        for u in uncertain(self):
            res += u.name + ' = ? (' + str(u.min) + ', ' + str(u.max) + ')\n'
        
        for k, v in self.__metadata.items():
            res += "%s: %s\n"%(str(k), str(v))
        return res
    
    def clone(self):
        m = Model(self.getData('title'))
        for r in reactions(self):
            setattr(m, r.name, Reaction(r.reagents, r.products, r.rate, r.irreversible))
        for p in parameters(self):
            setattr(m, p.name, constValue(p))
        for t in transformations(self):
            setattr(m, t.name, Transformation(t.rate))
        for s in self.__states:
            newdict= {}
            for i in s:
                newdict[i[0]]=i[1]
            setattr(m, s.name, StateArray(newdict, s.name))
        #handle uncertainties
        for u in uncertain(self):
            loc = u.name.split('.')
            if len(loc) >1:
                s = getattr(m,loc[0])
                var = getattr(s,loc[1])
                var.uncertainty(u.min, u.max)
            else:
                getattr(m, loc[0]).uncertainty(u.min, u.max)
        for k,v in self.__metadata.items():
            m.setData(k,v)
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
                            if not findWithName(name, self.__extvariables):
                                self.__extvariables.append(Variable(name))
                        else:
                            self.__variables.append(Variable(name))

    def set_uncertain(self, uncertainparameters):
        self.__m_Parameters = uncertainparameters
    


#----------------------------------------------------------------------------
#         Queries for Model network collections
#----------------------------------------------------------------------------

def variables(model):
    return model._Model__variables

def varnames(model):
    return [i.name for i in variables(model)]

def extvariables(model):
    return model._Model__extvariables

def reactions(model):
    return model._Model__reactions

def parameters(model):
    return model._Model__parameters

def transformations(model):
    return model._Model__transf

def uncertain(model):
    return list(iuncertain(model))

def iuncertain(model):
    for p in parameters(model):
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
    """Used to flag a wrong rate expression"""

class BadTypeComponent(Exception):
    """Used to flag an assignement of a model component to a wrong type object"""

def test():
    
    def force(A, t):
        return t*A
    register_kin_func(force)
##     pprint.pprint(__globs)

    m = Model('My first model')
    m.v1 = "A+B -> C", 3
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.v4 = "B   ->  "  , "2*4*step(t,at,top)"
    m.v5 = "C ->", "4.5*C*step(t,at,top)"
    m.t1 = transf("A*4 + C")
    m.t2 = transf("sqrt(2*A)")
    m.D  = variable("-2 * D")
    m.B  = 2.2
    m.myconstant = 2 * m.B / 1.1 # should be 4.0
    m.V3 = 0.5
    m.V3 = [0.1, 1.0]
    m.Km3 = 4
    m.Km3 = 1,6
    m.init = state(A = 1.0, C = 1, D = 1)
    m.afterwards = state(A = 1.0, C = 2, D = 1)
    m.afterwards.C.uncertainty(1,3)
    m.at = 1.0
    m.top = 2.0
    m.input2 = transf("4*step(t,at,top)")
    m.input3 = transf("force(top, t)")
    
    m.setData('where', 'in model')
    m.setData('title', 'My first model')
    
    print '********** Testing model construction and printing **********'
    print '------- result of model construction:\n'
    print m
    m2 = m.clone()
    print
    print '------- result of CLONING the model:\n'
    print m2
    
    print
    print '********** Testing iteration of components *****************'
    print 'iterating reactions(m)'
    for v in reactions(m):
        print v.name, ':', v.rate, '|', v.reagents, '->', v.products
    print '\niterating transformations(m)'
    for v in transformations(m):
        print v.name, ':', v.rate
    print '\niterating variables(m)'
    for x in variables(m):
        print x.name
    print '\niterating extvariables(m)'
    for x in extvariables(m):
        print x.name
    print '\niterating parameters(m)'
    for p in parameters(m):
        print p.name , '=',  p, 'bounds=', p.bounds
    print '\niterating uncertain(m)'
    for x in uncertain(m):
        print '\t', x.name, 'in (', x.min, ',', x.max, ')'
    
    print '\niterating m.init'
    for name, x in m.init:
        print '\t', name, '=', x
    print

    print '********** Testing component retrieval *********************'
    print 'm.K3 :',m.Km3
    print 'm.K3.name :',m.Km3.name, '(a float with a name attr)'
    print 'm.init:',m.init
    print 'm.init.A :',m.init.A
    print 'iterating m.init'
    for name, x in m.init:
        print '\t', name, '=', x


    print '********** Testing component reassignment *****************'
    print 'm.myconstant :',m.myconstant
    print len(parameters(m)), 'parameters total'
    print 'making m.myconstant = 5.0'
    m.myconstant = 5.0
    print 'm.myconstant :',m.myconstant
    print len(parameters(m)), 'parameters total'

    print 'making m.myconstant = react("A+B -> C"  , 3)'
    try:
        m.myconstant = react("A+B -> C"  , 3)
    except BadTypeComponent:
        print 'Failed! BadTypeComponent was caught.'
    print 'm.myconstant :',m.myconstant, '(still!)'
    print len(parameters(m)), 'parameters total'
    print
    print 'm.V3 :', m.V3
    print 'm.V3.bounds:' , m.V3.bounds
    print 'iterating m.uncertain'
    for x in uncertain(m):
        print '\t', x.name, 'in (', x.min, ',', x.max, ')'
    print len(uncertain(m)), 'uncertain parameters total'
    print 'making m.V3 = [0.1, 0.2]'
    m.V3 = [0.1, 0.2]
    print 'm.V3 :', m.V3
    print 'm.V3.bounds:' ,m.V3.bounds
    print len(uncertain(m)), 'uncertain parameters total'
    print 'making m.V4 = [0.1, 0.6]'
    m.V4 = [0.1, 0.6]
    print 'm.V4 :', m.V4
    print 'm.V4.bounds:' ,m.V4.bounds
    print len(uncertain(m)), 'uncertain parameters total'
    print 'iterating m.uncertain'
    for x in uncertain(m):
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
    m.Km3.uncertainty(None)

if __name__ == "__main__":
    test()
 