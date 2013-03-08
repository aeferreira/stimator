#!/usr/bin/env python
# -*- coding: latin1-*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator Model related classes
# Copyright António Ferreira 2006-2013
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

#----------------------------------------------------------------------------
#         Regular expressions for stoichiometry patterns
#----------------------------------------------------------------------------
fracnumberpattern = r"[-]?\d*[.]?\d+"
realnumberpattern = fracnumberpattern + r"(e[-]?\d+)?"
fracnumber = re.compile(fracnumberpattern, re.IGNORECASE)
realnumber = re.compile(realnumberpattern, re.IGNORECASE)

stoichiompattern   = r"^\s*(?P<reagents>.*)\s*(?P<irreversible>->|<=>)\s*(?P<products>.*)\s*$"
chemcomplexpattern = r"^\s*(?P<coef>("+realnumberpattern+")?)\s*(?P<variable>[_a-z]\w*)\s*$"

stoichiom  = re.compile(stoichiompattern,    re.IGNORECASE)
chemcomplex = re.compile(chemcomplexpattern, re.IGNORECASE)

identifier = re.compile(r"[_a-z]\w*", re.IGNORECASE)

def identifiersInExpr(_expr):
    iterator = identifier.finditer(_expr)
    return [_expr[m.span()[0]:m.span()[1]] for m in iterator]

#----------------------------------------------------------------------------
#         Utility functions
#----------------------------------------------------------------------------
def isNumber(value):
    for numtype in (float,int,long):
        if isinstance(value, numtype):
            return True
    return False


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
    def __eq__(self, other):
        if get_name(self) != get_name(other):
            return False
        if len(self.__metadata) != len(other.__metadata):
            return False
        for k in self.__metadata:
            if repr(self.__metadata[k]) != repr(other.__metadata[k]):
                return False
        return True

def setPar(obj, name, value):
    # accept only ContsValue's or Bounds
    if isNumber(value):
        value = constValue(float(value))
    elif isPairOfNums(value):
        value = Bounds(float(value[0]), float(value[1]))
    else:
        raise BadTypeComponent( "%s.%s"%(get_name(obj),name) + ' can not be assigned to ' + type(value).__name__)
    
    collection = obj.__dict__['_ownparameters']
    already_exists = name in collection
    if not already_exists:
        if isinstance(value, ConstValue):
            newvalue = constValue(value, name=name)
        else: #Bounds object
            newvalue = constValue((float(value.min)+float(value.max))/2.0, name=name)
            newvalue.bounds = value
            set_name(newvalue.bounds, name)
    else: #aready exists
        if isinstance(value, ConstValue):
            newvalue = constValue(value, name=name)
            newvalue.bounds = collection[name].bounds
        else: #Bounds object
            newvalue = constValue(collection[name], name=name)
            newvalue.bounds = value
            set_name(newvalue.bounds, name)
    collection[name] = newvalue

class _HasOwnParameters(ModelObject):
    def __init__(self, name = '?', parvalues = {}):
        ModelObject.__init__(self)
        self._ownparameters = {}
        if not isinstance(parvalues, dict):
            parvalues = dict(parvalues)
        for k,v in parvalues.items():
            self._ownparameters[k] = constValue(value = v, name = k)
    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        if name in self.__dict__['_ownparameters']:
            return self.__dict__['_ownparameters'][name]
        else:
            raise AttributeError( name + ' is not a member of '+ get_name(self))
    def __setattr__(self, name, value):
        if not name.startswith('_'):
            setPar(self, name, value)
        else:
            object.__setattr__(self, name, value)
    def __iter__(self):
        return iter(self._ownparameters.items())
    def __eq__(self, other):
        if not ModelObject.__eq__(self, other):
            return False
        if len(self._ownparameters) != len(other._ownparameters):
            return False
        for k in self._ownparameters:
            if not (self._ownparameters[k]) == (other._ownparameters[k]):
                return False
        return True

class _HasRate(_HasOwnParameters):
    def __init__(self, rate, parvalues = {}):
        _HasOwnParameters.__init__(self, parvalues = parvalues)
        self.__rate = rate.strip()
    def __str__(self):
        res =  "%s:\n  rate = %s\n" % (get_name(self), str(self()))
        if len(self._ownparameters) > 0:
            res += "  Parameters:\n"
            for k, v in self._ownparameters.items():
                res += "    %s = %g\n"% (k, v)
        return res
    def __call__(self, fully_qualified = False):
        rate = self.__rate
        if fully_qualified:
            for localparname in self._ownparameters:
                fully = '%s.%s'%(get_name(self), localparname)
##                 rate = rate.replace(localparname, fully)
                rate = re.sub(r"(?<!\.)\b%s\b(?![.\[])"%localparname, fully, rate)
        return rate
        
    def __eq__(self, other):
        if not _HasOwnParameters.__eq__(self, other):
            return False
        if self.__rate != other.__rate:
            return False
        return True

class Reaction(_HasRate):
    def __init__(self, reagents, products, rate, parvalues = {}, irreversible = False):
        _HasRate.__init__(self, rate, parvalues = parvalues)
        self._reagents = reagents
        self._products = products
        self._irreversible = irreversible
    def __str__(self):
        res = "%s:\n  reagents: %s\n  products: %s\n  rate    = %s\n" % (get_name(self), str(self._reagents), str(self._products), str(self())) 
        if len(self._ownparameters) > 0:
            res += "  Parameters:\n"
            for k, v in self._ownparameters.items():
                res += "    %s = %g\n"% (k, v)
        return res
    def __eq__(self, other):
        if not _HasRate.__eq__(self, other):
            return False
        if (self._reagents != other._reagents) or (self._products != other._products) or (self._irreversible != other._irreversible):
            return False
        return True

def react(stoichiometry, rate = 0.0, pars = {}):
    res = processStoich(stoichiometry)
    if not res:
        raise BadStoichError( "Bad stoichiometry definition:\n"+ stoichiometry)
    if isinstance(rate, float) or isinstance(rate, int):
        rate = massActionStr(rate, res[0])
    return Reaction(res[0], res[1], rate, pars, res[2])

class Variable_dXdt(_HasRate):
    def __init__(self, rate, parvalues = {}):
        _HasRate.__init__(self, rate, parvalues = parvalues)

class Transformation(_HasRate):
    def __init__(self, rate, parvalues = {}):
        _HasRate.__init__(self, rate, parvalues = parvalues)

def variable(rate = 0.0, pars = {}):
    if isinstance(rate, float) or isinstance(rate, int):
        rate = str(float(rate))
    return Variable_dXdt(rate, parvalues = pars)

def transf(rate = 0.0, pars = {}):
    if isinstance(rate, float) or isinstance(rate, int):
        rate = str(float(rate))
    return Transformation(rate, parvalues = pars)

class ConstValue(float,ModelObject):
    def __new__(cls, value):
        return float.__new__(cls,value)
    def initialize(self, aname = '?'):
        ModelObject.__init__(self, aname)
        self.bounds = None
    
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
    def __eq__(self, other):
        if repr(self) != repr(other):
            return False
        if isinstance(other, ConstValue):
            sbounds = self.bounds is not None
            obounds = other.bounds is not None
            if sbounds != obounds:
                return False
            if self.bounds is not None:
    ##             print "====== checking bounds of ConstValue %s" % get_name(self)
                if (self.bounds.min != other.bounds.min) or (self.bounds.max != other.bounds.max):
                    return False
    ##         print "#### eq ConstValue %s passed" % get_name(self)
        return True

def constValue(value = None, name = '?'):
    if isinstance(value, float) or isinstance(value, int):
        v = float(value)
        res = ConstValue(v)
        res.initialize(name)
    else:
        raise TypeError( value + ' is not a float or int')
    return res
        
class Bounds(ModelObject):
    def __init__(self, min = 0.0, max = 1.0):
        ModelObject.__init__(self)
        self.min = min
        self.max = max
    def __str__(self):
        return "(min=%f, max=%f)" % (self.min, self.max)

## class Variable(ModelObject):
##     def __init__(self, name):
##         ModelObject.__init__(self,name)
##     def __str__(self):
##         return get_name(self)

class StateArray(_HasOwnParameters):
    def __init__(self, varvalues, name):
        _HasOwnParameters.__init__(self,name, varvalues)
    def __str__(self):
        return '(%s)' % ", ".join(['%s = %s'% (x,str(float(value))) for (x,value) in  self._ownparameters.items()])

def state(**varvalues):
    return StateArray(varvalues, '?')

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
            return Reaction(res[0], res[1], rate, {}, res[2])
    return value

class Model(ModelObject):
    def __init__(self, title = ""):
        self.__dict__['_Model__reactions']         = QueriableList()
        self.__dict__['_Model__variables']         = []
        self.__dict__['_Model__extvariables']      = []
        self.__dict__['_ownparameters']            = {}
        self.__dict__['_Model__transf']            = QueriableList()
        self.__dict__['_Model__states']            = QueriableList()
        ModelObject.__init__(self)
        self.__dict__['_Model__m_Parameters']      = None
        self['title'] = title
        set_name(self, title)
    
    def __setattr__(self, name, value):
        if name in self.__dict__:
            object.__setattr__(self, name, value)
            return
        if isinstance(value,Variable_dXdt):
            react_name = 'd_%s_dt'% name
            stoich = ' -> %s'% name
            name = react_name # hope this works...
            value = react(stoich, value())
        value = _ConvertPair2Reaction(value)

        # find if the model has an existing  object with that name
        # start with strict types
        assoc = ((Reaction,       '_Model__reactions'),
                 (Transformation, '_Model__transf'),
                 (StateArray,     '_Model__states'))
        for t, listname in assoc:
            c = findWithNameIndex(name, self.__dict__[listname])
            if c > -1:
                if isinstance(value, t):
                    set_name(value,name)
                    self.__dict__[listname][c] = value
                    self.__refreshVars()
                    return
                else:
                    raise BadTypeComponent( name + ' can not be assigned to ' + type(value).__name__)
        # move on to parameters, accepting ConstValue, numbers or pairs of numbers
        if name in self.__dict__['_ownparameters']:
            setPar(self, name, value)
            
        # else append new object to proper list
        # start with strict types
        for t, listname in assoc:
            if isinstance(value, t):
                set_name(value,name)
                self.__dict__[listname].append(value)
                self.__refreshVars()
                return
        # move on to parameters, accepting ConstValue, numbers or pairs of numbers
        setPar(self, name, value)
        self.__refreshVars()

    def __getattr__(self, name):
        if name == '__m_Parameters':
            return self.__dict__['_Model__m_Parameters']
        if name in self.__dict__:
            return self.__dict__[name]
        if name in self._ownparameters:
            return self._ownparameters[name]
        c = self.__reactions.get(name)
        if c :
            return c
        if name in self.__variables:
            return name
        c = self.__transf.get(name)
        if c :
            return c
        c = self.__states.get(name)
        if c :
            return c
        raise AttributeError( name + ' is not defined for this model')
    
    def __findComponent(self, name):
        if name in self._ownparameters:
            return -1, 'parameters'
        c = findWithNameIndex(name, self.__reactions)
        if c>=0 :
            return c, 'reactions'
        try:
            c = self.__variables.index(name)
            return c, 'variables'
        except:
            pass
        c = findWithNameIndex(name, self.__transf)
        if c>=0 :
            return c, 'transf'
        raise AttributeError( name + ' is not a component in this model')

    def checkRates(self):
        for collection in (self.__reactions, self.__transf):
            for v in collection:
                resstring, value = _test_with_everything(v(),self, v)
                if resstring != "":
                    return False, '%s\nin rate of %s: %s' % (resstring, get_name(v), v())
        return True, 'OK'

    def __str__(self):
        check, msg = self.checkRates()
        if not check:
            raise BadRateError(msg)
        res = "%s\n"% self['title']
        #~ res = "%s\n"% self.title
        res += "\nVariables: %s\n" % " ".join(self.__variables)
        if len(self.__extvariables) > 0:
            res += "External variables: %s\n" % " ".join(self.__extvariables)
        for collection in (self.__reactions, self.__transf):
            for i in collection:
                res += str(i)
        for p in self.__states:
            res += get_name(p) +': '+ str(p) + '\n'
        mq = self()
        for p in mq.parameters:
            res += get_name(p) +' = '+ str(p) + '\n'
        for u in mq.uncertain:
            res += get_name(u) + ' = ? (' + str(u.min) + ', ' + str(u.max) + ')\n'
        for k in self._ModelObject__metadata:
            res += "%s: %s\n"%(str(k), str(self._ModelObject__metadata[k]))
        return res
    
    def clone(self, new_title = None):
        m = Model(self['title'])
        if new_title is not None:
            m['title'] = new_title
        for r in self.__reactions:
            setattr(m, get_name(r), Reaction(r._reagents, r._products, r(), r._ownparameters, r._irreversible))
        for p in self._ownparameters.values():
            setattr(m, get_name(p), constValue(p))
        for t in self.__transf:
            setattr(m, get_name(t), Transformation(t(), t._ownparameters))
        for s in self.__states:
            newdict = {}
            for i in s:
                newdict[i[0]]=i[1]
            setattr(m, get_name(s), StateArray(newdict, get_name(s)))
        #handle uncertainties
        for u in self().uncertain:
            loc = get_name(u).split('.')
            currobj = m
            for attribute in loc:
                currobj = getattr(currobj, attribute)
            currobj.uncertainty(u.min, u.max)
        m._ModelObject__metadata.update(self._ModelObject__metadata)
        return m
    
    def __eq__(self, other):
        if not ModelObject.__eq__(self, other):
            return False
        cnames = ('reactions', 'transf', 'states', 'pars', 'vars', 'extvars')
        collections1 = [self.__reactions, self.__transf, self.__states, self._ownparameters, self.__variables, self.__extvariables]
        collections2 = [other.__reactions, other.__transf, other.__states, other._ownparameters, other.__variables, other.__extvariables]
        for cname, c1,c2 in zip(cnames, collections1, collections2):
            if len(c1) != len(c2):
                return False
            if isinstance(c1, dict):
                names = c1.keys()
            else:
                names = [v for v in c1]
            for vname in names:
                if isinstance(vname, ModelObject):
                    vname = get_name(vname)
                r = getattr(self, vname)
                ro = getattr(other, vname)
                if not ro == r:
                    return False
        return True        
    
    def set_uncertain(self, uncertainparameters):
        self.__m_Parameters = uncertainparameters

    def _genlocs4rate(self, obj = None):
        for p in self._ownparameters.items():
            yield p
        collections = [self.__reactions, self.__transf]
        for c in collections:
            for v in c:
                yield (get_name(v), v)
        
        if (obj is not None) and (len(obj._ownparameters) > 0):
            for p in obj._ownparameters.items():
                yield p

    def __refreshVars(self):
        del(self.__variables[:]) #can't use self.__variables= [] : Triggers __setattr__
        del(self.__extvariables[:])
        for v in self.__reactions:
            for rp in (v._reagents, v._products):
                for (vname, coef) in rp:
                    if vname in self.__variables:
                        continue
                    else:
                        if vname in self._ownparameters:
                            if not vname in self.__extvariables:
                                self.__extvariables.append(vname)
                        else:
                            self.__variables.append(vname)
    def __call__(self):
        mq = ModelQuerier(self)
        return mq

def _test_with_everything(valueexpr, model, obj): 
    locs = {}
    for (name, value) in model._genlocs4rate(obj):
        locs[name] = value
##     print "\nvalueexpr:", valueexpr
##     print "---locs"
##     for k in locs:
##         if isinstance(locs[k], _HasRate):
##             print k, '-->', type(locs[k])
##         else:
##             print k, '-->', locs[k]
    
    #part 1: nonpermissive, except for NameError
    try :
       value = float(eval(valueexpr, __globs, locs))
    except NameError:
       pass
    except Exception, e:
##        print 'returned on first pass'
       return ("%s : %s"%(str(e.__class__.__name__),str(e)), 0.0)
    #part 2: permissive, with dummy values (1.0) for vars
    vardict = {}
    for i in model._Model__variables:
        vardict[i]=1.0
    vardict['t'] = 1.0
    locs.update(vardict)
    try :
       value = float(eval(valueexpr, __globs, locs))
    except (ArithmeticError, ValueError):
       pass # might fail but we don't know the values of vars
    except Exception, e:
##        print 'returned on second pass'
       return ("%s : %s"%(str(e.__class__.__name__),str(e)), 0.0)
    return "", value

#----------------------------------------------------------------------------
#         Queries for Model network collections
#----------------------------------------------------------------------------
## def varnames(model):
##     return model._Model__variables

class QueriableList(list):
    def get(self, aname):
        for i in self:
            if get_name(i) == aname:
                return i
        return None
        

class ModelQuerier(object):
    """A class to query a model object as a whole."""
    
    def __init__(self, model):
        self.m = model
    
    def get_varnames(self):
        return self.m._Model__variables
    varnames = property(get_varnames)
    
    def get_extvariables(self):
        return self.m._Model__extvariables
    extvariables = property(get_extvariables)

    def get_reactions(self):
        return self.m._Model__reactions
    reactions = property(get_reactions)

    def get_parameters(self):
        return QueriableList(self.get_iparameters())
    parameters = property(get_parameters)

    def get_iparameters(self):
        for p in self.m._ownparameters.values():
            yield p
        collections = [self.m._Model__reactions, self.m._Model__transf]
        for c in collections:
            for v in c:
                for iname, value in v._ownparameters.items():
                    ret = constValue(value)
                    retname = get_name(v) + '.' + iname
                    ret.initialize(retname)
                    if value.bounds:
                        b = Bounds(value.bounds.min, value.bounds.max)
                        set_name(b, retname)
                        ret.bounds = b
                    yield ret
    iparameters = property(get_iparameters)

    def get_transformations(self):
        return self.m._Model__transf
    transformations = property(get_transformations)

    def get_uncertain(self):
        return QueriableList(self.get_iuncertain())
    uncertain = property(get_uncertain)

    def get_iuncertain(self):
        for p in self.get_iparameters():
            if p.bounds:
                yield p.bounds
        for s in self.m._Model__states:
            for iname, value in s:
                if value.bounds:
                    ret = Bounds(value.bounds.min, value.bounds.max)
                    set_name(ret, get_name(s) + '.' + iname)
                    yield ret
    iuncertain = property(get_iuncertain)

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
    m.v2 = react("    -> 4.5 A"  , rate = math.sqrt(4.0)/2)
    v3pars = (('V3',0.5),('Km', 4))
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)", pars = v3pars)
##     m.v3.V3 = 0.5
##     m.v3.Km3 = 4
    m.v3.V3 = [0.1, 1.0]
    m.v4 = "B   ->  "  , "2*4*step(t,at,top)"
    m.v5 = "C ->", "4.5*C*step(t,at,top)"
    m.t1 = transf("A*Vt + C", dict(Vt=4))
    m.t2 = transf("sqrt(2*A)")
    m.D  = variable("-2 * D")
    m.B  = 2.2
    m.myconstant = 2 * m.B / 1.1 # should be 4.0
    m.V3 = 0.6
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
    print '********** Testing equality of models *****************'
    print "m2 == m"
    print m2 == m

    print
    print '********** Testing iteration of components *****************'
    print 'iterating m().reactions'
    for v in m().reactions:
        print get_name(v), ':', v(), '|', v._reagents, '->', v._products
    print '\niterating m().reactions with fully qualified rates'
    for v in m().reactions:
        print get_name(v), ':', v(fully_qualified = True), '|', v._reagents, '->', v._products
    print '\niterating m().transformations'
    for v in m().transformations:
        print get_name(v), ':', v()
    print '\niterating m().transformations with fully qualified rates'
    for v in m().transformations:
        print get_name(v), ':', v(fully_qualified = True)
    print '\niterating m().varnames:'
    for x in m().varnames:
        print x
    print '\niterating m().extvariables'
    for x in m().extvariables:
        print x
    print '\niterating m().parameters'
    for p in m().parameters:
        print get_name(p) , '=',  p, '\n  bounds=', p.bounds
    print '\niterating m().uncertain'
    for x in m().uncertain:
        print '\t', get_name(x), 'in (', x.min, ',', x.max, ')'
    
    print '\niterating m.init (a state)'
    for xname, x in m.init:
        print '\t', xname, '=', x

    print '\niterating m.v3 (iterates parameters returning (name,value) tuples)'
    for xname, x in m.v3:
        print '\t', xname, '=', x
    print

    print '********** Testing component retrieval *********************'
    print 'm.K3 :',m.Km3
    print 'get_name(m.K3) :',get_name(m.Km3)
    print 'm.v3.K :',m.v3.Km
    print 'get_name(m.v3.K) :',get_name(m.v3.Km)
    print 'm.init:',m.init
    print 'm.init.A :',m.init.A
    print 'get_name(m.init.A) :',get_name(m.init.A)
    try:
        print 'm.init.B :',m.init.B
    except AttributeError:
        print 'm.init.B access raised exception AttributeError'
    print 'iterating m.init returning (name,value) tuples'
    for xname, x in m.init:
        print '\t', xname, '=', x

    print '********** Testing component reassignment *****************'
    print 'm.myconstant :',m.myconstant
    print len(m().parameters), 'parameters total'
    print '\nmaking m.myconstant = 5.0'
    m.myconstant = 5.0
    print 'm.myconstant :',m.myconstant
    print len(m().parameters), 'parameters total'

    print '\nmaking m.myconstant = react("A+B -> C"  , 3)'
    try:
        m.myconstant = react("A+B -> C"  , 3)
    except BadTypeComponent:
        print 'Failed! BadTypeComponent was caught.'
    print 'm.myconstant :',m.myconstant, '(still!)'
    print len(m().parameters), 'parameters total'
    print '\nmaking m.v2 = 3.14'
    try:
        m.v2 = 3.14
    except BadTypeComponent:
        print 'Failed! BadTypeComponent was caught.'
    print 'm.v2 :',type(m.v2), '\n(still a Reaction)'
    print len(m().parameters), 'parameters total'
    print
    print 'm.V3 :', m.V3
    print 'm.V3.bounds:' , m.V3.bounds
    print 'iterating m.uncertain'
    for x in m().uncertain:
        print '\t', get_name(x), 'in (', x.min, ',', x.max, ')'
    print len(m().uncertain), 'uncertain parameters total'
    print '\nmaking m.V3 = [0.1, 0.2]'
    m.V3 = [0.1, 0.2]
    print 'm.V3 :', m.V3
    print 'm.V3.bounds:' ,m.V3.bounds
    print len(m().uncertain), 'uncertain parameters total'
    print '\nmaking m.V4 = [0.1, 0.6]'
    m.V4 = [0.1, 0.6]
    print 'm.V4 :', m.V4
    print 'm.V4.bounds:' ,m.V4.bounds
    print len(m().uncertain), 'uncertain parameters total'
    print 'iterating m.uncertain'
    for x in m().uncertain:
        print '\t', get_name(x), 'in (', x.min, ',', x.max, ')'
    print '\nmaking m.V4 = 0.38'
    m.V4 = 0.38
    print 'm.V4 :', m.V4
    print 'm.V4.bounds:' ,m.V4.bounds
    print len(m().uncertain), 'uncertain parameters total'
    print 'iterating m.uncertain'
    for x in m().uncertain:
        print '\t', get_name(x), 'in (', x.min, ',', x.max, ')'
    print '\nmaking m.init.A = 5.0'
    m.init.A = 5.0
    print 'iterating m.init'
    for xname, x in m.init:
        print '\t', xname, '=', x.pprint()
    print '\nflagging init.A as uncertain with   m.init.A = (0.5, 2.5)'
    m.init.A = (0.5, 2.5)
    print 'iterating m.init'
    for xname, x in m.init:
        print '\t', xname, '=', x.pprint()
    print '\ncalling    m.init.A.uncertainy(0.5,3.0)'
    m.init.A.uncertainty(0.5,3.0)
    print 'iterating m.init'
    for xname, x in m.init:
        print '\t', xname, '=', x.pprint()
    print '\ncalling    m.init.A.uncertainy(None)'
    m.init.A.uncertainty(None)
    print 'iterating m.init'
    for xname, x in m.init:
        print '\t', xname, '=', x.pprint()
    print '\nmaking m.init.A back to 1.0'
    m.init.A = 1.0
    print 'iterating m.init'
    for xname, x in m.init:
        print '\t', xname, '=', x.pprint()
    print 


if __name__ == "__main__":
    test()
 