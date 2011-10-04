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

#----------------------------------------------------------------------------
#         Functions to check the validity of math expressions
#----------------------------------------------------------------------------
__globs = {}
__haskinetics = {}
for k, v in globals().items():
    if hasattr(v,"is_rate"):
        __haskinetics[k] = v
__globs.update(__haskinetics)
__globs.update(vars(math))

def register_kin_func(f):
    f.is_rate = True
    __globs[f.__name__] = f
    globals()[f.__name__] = f

def _test_with_everything(valueexpr, model): 
    locs = {}
    for p in parameters(model):
        locs[get_name(p)] = p #.value
    
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
        vardict[get_name(i)]=1.0
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
#         Regular expressions for stoichiometry patterns
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

def findWithName(aname, alist):
    for i in alist:
        if get_name(i) == aname:
            return i
    return None

def findWithNameIndex(aname, alist):
    for i,elem in enumerate(alist):
        if get_name(elem) == aname:
            return i
    return -1

#----------------------------------------------------------------------------
#         Model and Model component classes
#----------------------------------------------------------------------------

def get_name(obj):
    return obj._ModelObject__name

def set_name(obj, name):
    obj._ModelObject__name = name

class ModelObject(object):
    def __init__(self, name = '?'):
        self.__dict__['_ModelObject__metadata'] = {}
        self.__dict__['_ModelObject__name'] = name
    
    def __setitem__(self, key, value):
        if isinstance(key, str) or isinstance(key, unicode):
            self.__metadata[key] = value
        else:
            raise TypeError( "Keys must be strings.")

    def __delitem__(self, key):
        if isinstance(key, str) or isinstance(key, unicode):
            if self.__metadata.has_key(key):
                del(self.__metadata[key])
        else:
            raise TypeError( "Keys must be strings.")

    def __getitem__(self, key):
        """retrieves info by name"""
        if isinstance(key, str) or isinstance(key, unicode):
            if not key in self.__metadata:
                return None
            return self.__metadata[key]
        else:
            raise TypeError( "Keys must be strings.")

    
class _HasRate(ModelObject):
    def __init__(self, rate):
        ModelObject.__init__(self)
        self.__rate = rate.strip()
    def __str__(self):
        return "%s:\n  rate = %s\n" % (get_name(self), str(self()))
    def __call__(self):
        return self.__rate


class Reaction(_HasRate):
    def __init__(self, reagents, products, rate, irreversible = False):
        _HasRate.__init__(self, rate)
        self._reagents = reagents
        self._products = products
        self._irreversible = irreversible
    def __str__(self):
        return "%s:\n  reagents: %s\n  products: %s\n  rate    = %s\n" % (get_name(self), str(self._reagents), str(self._products), str(self())) 

def react(stoichiometry, rate = 0.0):
    res = processStoich(stoichiometry)
    if not res:
        raise BadStoichError( "Bad stoichiometry definition:\n"+ stoichiometry)
    if isinstance(rate, float) or isinstance(rate, int):
        rate = massActionStr(rate, res[0])
    return Reaction(res[0], res[1], rate, res[2])


class Variable_dXdt(_HasRate):
    def __init__(self, rate):
        _HasRate.__init__(self, rate)

class Transformation(_HasRate):
    def __init__(self, rate):
        _HasRate.__init__(self, rate)

def variable(rate = 0.0):
    if isinstance(rate, float) or isinstance(rate, int):
        rate = str(float(rate))
    return Variable_dXdt(rate)

def transf(rate = 0.0):
    if isinstance(rate, float) or isinstance(rate, int):
        rate = str(float(rate))
    return Transformation(rate)

class ConstValue(float,ModelObject):
    def __new__(cls, value):
        return float.__new__(cls,value)
    def initialize(self, aname = '?', into = None):
        ModelObject.__init__(self, aname)
        self.bounds = None
        if into:
            set_name(self, get_name(into))
            self.bounds = into.bounds
    
    def uncertainty(self, *pars):
        if len(pars) == 0 or pars[0]==None:
            self.bounds = None
            return
        if len(pars) != 2:
            return #TODO raise exception
        self.bounds = Bounds(float(pars[0]), float(pars[1]))
        set_name(self.bounds, get_name(self))
    def pprint(self):
        res = float.__str__(self)
        if self.bounds:
            res+= " ? (min = %f, max=%f)" % (self.bounds.min, self.bounds.max)
        return res

def constValue(value = None, name = '?', into = None):
    if isinstance(value, float) or isinstance(value, int):
        v = float(value)
        res = ConstValue(v)
        res.initialize(name)
    elif (isinstance(value, tuple) or isinstance(value, list)) and len(value)==2:
        bounds = Bounds(float(value[0]), float(value[1]))
        v = (bounds.min + bounds.max)/2.0
        if into:
            res = ConstValue(into)
            res.initialize(name, into = into)
        else:
            res = ContsValue(v)
            res.initialize(name)
        res.bounds = bounds
    else:
        raise TypeError( value + ' is not a float or pair of floats')
    if res.bounds:
        set_name(res.bounds, get_name(res))
    return res
        

class Bounds(ModelObject):
    def __init__(self, min = 0.0, max = 1.0):
        ModelObject.__init__(self)
        self.min = min
        self.max = max
    def __str__(self):
        return "(min=%f, max=%f)" % (self.min, self.max)

class Variable(ModelObject):
    def __init__(self, name):
        ModelObject.__init__(self,name)
    def __str__(self):
        return get_name(self)

class StateArray(ModelObject):
    def __init__(self, varvalues, name):
        ModelObject.__init__(self,name)
        for k,v in varvalues.items():
            varvalues[k] = constValue(value = v, name = k)
        self.__dict__['_varvalues'] = varvalues
    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        if name in self.__dict__['_varvalues']:
            return self.__dict__['_varvalues'][name]
        else:
            raise AttributeError( name + ' is not a member of state '+ get_name(self))
    def __setattr__(self, name, value):
        if name != '_varvalues' and name != '_ModelObject__name':
            value = constValue(value = value, name = name, into = self.__dict__['_varvalues'].get(name, None))
            self.__dict__['_varvalues'][name]=value
        else:
            object.__setattr__(self, name, value)
    def __str__(self):
        return '(%s)' % ", ".join(['%s = %s'% (x,str(float(value))) for (x,value) in  self._varvalues.items()])
    def __iter__(self):
        return iter(self._varvalues.items())

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

def _ConvertPair2Reaction(value):
    if (isinstance(value, tuple) or isinstance(value, list)) and len(value)==2:
        if isinstance(value[0], str):
            good2nd = False
            for numtype in (str, float,int,long):
                if isinstance(value[1], numtype):
                    good2nd = True
                    break
            if not good2nd:
                return value
            res = processStoich(value[0])
            if not res:
                raise BadStoichError( "Bad stoichiometry definition:\n"+ value[0])
            rate = value[1]
            for numtype in (float,int,long):
                if isinstance(rate, numtype):
                    rate = massActionStr(rate, res[0])
            return Reaction(res[0], res[1], rate, res[2])
    return value

class Model(ModelObject):
    def __init__(self, title = ""):
        self.__dict__['_Model__reactions']         = []
        self.__dict__['_Model__variables']         = []
        self.__dict__['_Model__extvariables']      = []
        self.__dict__['_Model__parameters']        = []
        self.__dict__['_Model__transf']            = []
        self.__dict__['_Model__states']            = []
        ModelObject.__init__(self)
        self.__dict__['_Model__m_Parameters']      = None
        self['title'] = title
        set_name(self, title)
    
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
            value = react(stoich, value())
        value = _ConvertPair2Reaction(value)

        # find if the model has an existing  object with that name
        assoc = ((Reaction,       '_Model__reactions'),
                 (ConstValue,     '_Model__parameters'),
                 (Transformation, '_Model__transf'),
                 (StateArray,     '_Model__states'))
        for t, listname in assoc:
            c = findWithNameIndex(name, self.__dict__[listname])
            if c > -1:
                if isinstance(value, t) or (isinstance(value, Bounds) and listname == '_Model__parameters'):
                    set_name(value,name)
                    if isinstance(value, Bounds):
                        self.__dict__[listname][c].bounds = value
                    else:
                        self.__dict__[listname][c] = value
                    self.__refreshVars()
                    return
                else:
                    raise BadTypeComponent( name + ' can not be assigned to ' + type(value).__name__)
        # else append new object to proper list
        for t, listname in assoc:
            if isinstance(value, t):
                set_name(value,name)
                self.__dict__[listname].append(value)
                self.__refreshVars()
                return
        if isinstance (value, Bounds):
            set_name(value,name)
            newvalue = constValue((float(value.min)+float(value.max))/2.0, name=name)
            newvalue.bounds = value
            self.__dict__['_Model__parameters'].append(newvalue)
            self.__refreshVars()
            return
        # Fail safe. This is a first level model atribute
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
            return get_name(c)
        c = findWithName(name, self.__transf)
        if c :
            return c
        c = findWithName(name, self.__states)
        if c :
            return c
        raise AttributeError( name + ' is not defined for this model')
    
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
                resstring, value = _test_with_everything(v(),self)
                if resstring != "":
                    return False, '%s\nin rate of %s: %s' % (resstring, get_name(v), v())
        return True, 'OK'

    def __str__(self):
        check, msg = self.checkRates()
        if not check:
            raise BadRateError(msg)
        res = "%s\n"% self['title']
        #~ res = "%s\n"% self.title
        res += "\nVariables: %s\n" % " ".join([get_name(i) for i in self.__variables])
        if len(self.__extvariables) > 0:
            res += "External variables: %s\n" % " ".join([get_name(i) for i in self.__extvariables])
        for collection in (self.__reactions, self.__transf):
            for i in collection:
                res += str(i)
        for p in self.__states:
            res += get_name(p) +': '+ str(p) + '\n'
        for p in self.__parameters:
            res += get_name(p) +' = '+ str(p) + '\n'
        for u in uncertain(self):
            res += get_name(u) + ' = ? (' + str(u.min) + ', ' + str(u.max) + ')\n'
        
        for k, v in self._ModelObject__metadata.items():
            res += "%s: %s\n"%(str(k), str(v))
        return res
    
    def clone(self):
        m = Model(self['title'])
        for r in reactions(self):
            setattr(m, get_name(r), Reaction(r._reagents, r._products, r(), r._irreversible))
        for p in parameters(self):
            setattr(m, get_name(p), constValue(p))
        for t in transformations(self):
            setattr(m, get_name(t), Transformation(t()))
        for s in self.__states:
            newdict= {}
            for i in s:
                newdict[i[0]]=i[1]
            setattr(m, get_name(s), StateArray(newdict, get_name(s)))
        #handle uncertainties
        for u in uncertain(self):
            loc = get_name(u).split('.')
            if len(loc) >1:
                s = getattr(m,loc[0])
                var = getattr(s,loc[1])
                var.uncertainty(u.min, u.max)
            else:
                getattr(m, loc[0]).uncertainty(u.min, u.max)
        for k,v in self._ModelObject__metadata.items():
            m[k] = v
        return m
    
    def __refreshVars(self):
        del(self.__variables[:]) #can't use self.__variables= [] : Triggers __setattr__
        del(self.__extvariables[:])
        for v in self.__reactions:
            for rp in (v._reagents, v._products):
                for (vname, coef) in rp:
                    if findWithName(vname, self.__variables):
                        continue
                    else:
                        if findWithName(vname, self.__parameters):
                            if not findWithName(vname, self.__extvariables):
                                self.__extvariables.append(Variable(vname))
                        else:
                            self.__variables.append(Variable(vname))

    def set_uncertain(self, uncertainparameters):
        self.__m_Parameters = uncertainparameters
    


#----------------------------------------------------------------------------
#         Queries for Model network collections
#----------------------------------------------------------------------------

def variables(model):
    return model._Model__variables

def varnames(model):
    return [get_name(i) for i in variables(model)]

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
        for iname, value in s:
##             print iname, value
            if value.bounds:
                ret = Bounds(value.bounds.min, value.bounds.max)
                set_name(ret, get_name(s) + '.' + iname)
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
    
    m['where'] = 'in model'
    m['for what'] = 'testing'
    
    print '********** Testing model construction and printing **********'
    print '------- result of model construction:\n'
    print m
    print "access to info as keys---------"
    print "m['for what'] =", m['for what']
    del(m['where'])
    print "\nafter del(m['where'])"
    print "m['where'] =", m['where']
    
    m2 = m.clone()
    print
    print '------- result of CLONING the model:\n'
    print m2
    
    print
    print '********** Testing iteration of components *****************'
    print 'iterating reactions(m)'
    for v in reactions(m):
        print get_name(v), ':', v(), '|', v._reagents, '->', v._products
    print '\niterating transformations(m)'
    for v in transformations(m):
        print get_name(v), ':', v()
    print '\niterating variables(m)'
    for x in variables(m):
        print get_name(x)
    print '\niterating extvariables(m)'
    for x in extvariables(m):
        print get_name(x)
    print '\niterating parameters(m)'
    for p in parameters(m):
        print get_name(p) , '=',  p, 'bounds=', p.bounds
    print '\niterating uncertain(m)'
    for x in uncertain(m):
        print '\t', get_name(x), 'in (', x.min, ',', x.max, ')'
    
    print '\niterating m.init'
    for xname, x in m.init:
        print '\t', xname, '=', x
    print

    print '********** Testing component retrieval *********************'
    print 'm.K3 :',m.Km3
    print 'get_name(m.K3) :',get_name(m.Km3)
    print 'm.init:',m.init
    print 'm.init.A :',m.init.A
    try:
        print 'm.init.B :',m.init.B
    except AttributeError:
        print 'm.init.B access raised exception AttributeError'
    print 'iterating m.init'
    for xname, x in m.init:
        print '\t', xname, '=', x


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
        print '\t', get_name(x), 'in (', x.min, ',', x.max, ')'
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
        print '\t', get_name(x), 'in (', x.min, ',', x.max, ')'
    print 'making m.init.A = 5.0'
    m.init.A = 5.0
    print 'iterating m.init'
    for xname, x in m.init:
        print '\t', xname, '=', x.pprint()
    print 'flagging init.A as uncertain with   m.init.A = (0.5, 2.5)'
    m.init.A = (0.5, 2.5)
    print 'iterating m.init'
    for xname, x in m.init:
        print '\t', xname, '=', x.pprint()
    print 'calling    m.init.A.uncertainy(0.5,3.0)'
    m.init.A.uncertainty(0.5,3.0)
    print 'iterating m.init'
    for xname, x in m.init:
        print '\t', xname, '=', x.pprint()
    print 'calling    m.init.A.uncertainy(None)'
    m.init.A.uncertainty(None)
    print 'iterating m.init'
    for xname, x in m.init:
        print '\t', xname, '=', x.pprint()
    print 'making m.init.A back to 1.0'
    m.init.A = 1.0
    print 'iterating m.init'
    for xname, x in m.init:
        print '\t', xname, '=', x.pprint()
    print 
    m.Km3.uncertainty(None)

if __name__ == "__main__":
    test()
 