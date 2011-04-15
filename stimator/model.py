#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator Model related classes
# Copyright António Ferreira 2006-2009
#----------------------------------------------------------------------------
import sys
import os
import os.path
import re
import math
from kinetics import *
from numpy import *
from timeit import Timer
import pprint

#----------------------------------------------------------------------------
#         Functions to check the validity of math expressions
#----------------------------------------------------------------------------
__globs = vars(math)
__haskinetics = {}
for k, v in globals().items():
    if hasattr(v,"israte"):
        __haskinetics[k] = v
        
__globs.update(__haskinetics)
## pprint.pprint(__globs)

def _test_with_everything(valueexpr, model): 
    locs = {}
    for p in model.parameters:
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
    for i in model.variables:
        vardict[i.name]=1.0
    locs.update(vardict)
    for i in model.forcing:
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


class Forcing(object):
    def __init__(self, rate):
        self.rate = rate
        self.name = '?'
    def __str__(self):
        return "%s:\n  rate = %s\n" % (self.name, str(self.rate))

    def genForcingFunction(self, model):
        f = self.rate
        if not callable(f):
            if not isinstance(f,str):
                raise TypeError( 'rate must be a string or a callable.')
            ratebytecode = compile(model.rateCalcString(f, with_uncertain=False), '<string>','eval')
        
        else:
            cc = f.func_code
            nargs = cc.co_argcount
            argnames = cc.co_varnames[:nargs]
            codes = []
            for a in argnames:
                if a == 't':
                    codes.append(('t', 0))
                    continue
                i, collection = model.findComponent(a)
                if collection == 'parameters':
                    codes.append(('p', getattr(model, a)))
                elif collection == 'variables':
                    codes.append(('v', i))
                else:
                    raise AttributeError( a + ' is not a component in of the model')
            args = [0.0] * nargs
            en = [(i, c,d) for (i,(c,d)) in enumerate(codes)]
            for i,c,d in en:
                if c == 'p':
                    args[i]=d

        def reteval(variables, t):
            return eval(ratebytecode)

        def retf(variables, t):
            for i,c,d in en:
                if c == 't':
                    args[i] = t
                elif c == 'v':
                    args[i] = variables[d]
                else:
                    continue
            return f(*args)
                
        if callable(f):
            self.func = retf
        else:
            self.func = reteval

def forcing(rate):
    if isinstance(rate, str) or callable(rate):
        return Forcing(rate)
    else:
        raise ValueError('Bad forcing function')


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
        self.__dict__['_Model__forcing']           = []
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
                 (StateArray,     '_Model__states'),
                 (Forcing,        '_Model__forcing'))
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
                elif (isinstance(value, Forcing) and listname == '_Model__parameters'):
                    #delete the constant parameter. It is going to be substituted by a Forcing object
                    del (self.__dict__['_Model__parameters'][c])
                    break
                elif (isinstance(value, ConstValue) and listname == '_Model__forcing'):
                    #delete the forcing function. It is going to be substituted by constant parameter
                    del (self.__dict__['_Model__forcing'][c])
                    break
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
        c = findWithName(name, self.__forcing)
        if c :
            return c
        raise AttributeError( name + ' is not defined for this model')
    
    def setData(self, name, value):
        self.__metadata[name] = value
    
    def getData(self, name):
        if not name in self.__metadata:
            return None
        return self.__metadata[name]
        

    def findComponent(self, name):
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
        c = findWithNameIndex(name, self.__forcing)
        if c>=0 :
            return c, 'forcing'
        raise AttributeError( name + ' is not a component in this model')

    def vectorize(self, state):
        if isinstance(state, str) or isinstance(state, unicode):
            if not hasattr(self, state):
                raise AttributeError( state + ' is not defined for this model')
            state = getattr(self, state)
        newlist = [state.varvalues.get(var.name,0.0) for var in self.__variables]
        return array(newlist)

    def varnames(self):
        return [i.name for i in self.__variables]
    
    def checkRates(self):
        for collection in (self.__reactions, self.__transf):
            for v in collection:
                resstring, value = _test_with_everything(v.rate,self)
                if resstring != "":
                    return False, '%s\nin rate of %s: %s' % (resstring, v.name, v.rate)
        for t in self.__forcing:
            t.genForcingFunction(self)
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
        for collection in (self.__reactions, self.__transf, self.__forcing):
            for i in collection:
                res += str(i)
        for p in self.__states:
            res += p.name +': '+ str(p) + '\n'
        for p in self.__parameters:
            res += p.name +' = '+ str(p) + '\n'
        for u in self.uncertain:
            res += u.name + ' = ? (' + str(u.min) + ', ' + str(u.max) + ')\n'
        
        for k, v in self.__metadata.items():
            res += "%s: %s\n"%(str(k), str(v))
        return res
    
    def clone(self):
        m = Model(self.getData('title'))
        for r in self.reactions:
            setattr(m, r.name, Reaction(r.reagents, r.products, r.rate, r.irreversible))
        for p in self.parameters:
            setattr(m, p.name, constValue(p))
        for t in self.transf:
            setattr(m, t.name, Transformation(t.rate))
        for f in self.forcing:
            setattr(m, f.name, Forcing(f.rate))
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


    def genStoichiometryMatrix(self):
        check, msg = self.checkRates()
        if not check:
            raise BadRateError(msg)

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
        # replace forcing functions
        for i,v in enumerate(self.forcing):
            rateString = re.sub(r"\b"+ v.name+r"\b", "forcing_function[%d]"%i, rateString)
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
    
    def rates_strings(self):
        """Generate a tuple of tuples of
           (name, rate) where
           'name' is the name of a reaction
           'rhs' is the string of the rate of the reaction.
        """
        check, msg = self.checkRates()
        if not check:
            raise BadRateError(msg)

        return tuple([(v.name, v.rate) for v in self.__reactions])

    def dXdt_strings(self):
        """Generate a tuple of tuples of
           (name, rhs) where
           'name' is the name of a variable
           'rhs' is the string of the rhs of that variable
           in the SODE for this model.
        """
        check, msg = self.checkRates()
        if not check:
            raise BadRateError(msg)
        N = self.genStoichiometryMatrix()
        res = []
        for i,x in enumerate(self.variables):
            name = x.name
            dXdtstring = ''
            for j,v in enumerate(self.reactions):
                coef = N[i,j]
                if coef == 0.0: continue
                ratestring = '(%s)'% v.rate
                if coef == 1.0:
                    ratestring = '+'+ratestring
                else:
                    ratestring = "%g*%s" % (coef,ratestring)
                    if coef > 0.0:
                        ratestring = '%s%s'%('+', ratestring)
                dXdtstring += ratestring
            res.append((name, dXdtstring))
        return tuple(res)

    def Jacobian_strings(self, _scale = 1.0):
        """Generate a matrix (list of lists) of strings
           to compute the jacobian for this model.
        
           IMPORTANT: sympy module must be installed!"""

        check, msg = self.checkRates()
        if not check:
            raise BadRateError(msg)
        try:
            import sympy
        except:
            print 'ERROR: sympy module must be installed to generate Jacobian strings'
            raise
        _dxdtstrings = self.dXdt_strings()
        _symbs = {}
        for x in self.__variables:
            _symbs[x.name] = sympy.Symbol(str(x.name))
        for p in self.__parameters:
            _symbs[p.name] = sympy.Symbol(str(p.name))
        for f in self.__forcing:
            _symbs[f.name] = sympy.Symbol(str(f.name))
        _nvars = len(self.__variables)
        _nvars = len(self.__variables)
        _jfuncs = []
        for _i in range(_nvars):
            _jfuncs.append([])
            _ids = identifiersInExpr(_dxdtstrings[_i][1])
            if len(_ids) == 0:
                for _j in range(_nvars):
                    _jfuncs[_i].append('0.0')
            else:
            
                for _j in range(_nvars):
                    _res = eval(_dxdtstrings[_i][1], None, _symbs)
                    _res = _res * _scale
                    _dres = str(sympy.diff(_res, _symbs[self.__variables[_j].name]))
                    if _dres == '0':
                        _dres = '0.0'
                    _jfuncs[_i].append(_dres)
                    #~ print _jfuncs
        return _jfuncs
        
    def dfdp_strings(self, _parnames, _scale = 1.0):
        """Generate a matrix (list of lists) of strings
           to compute the partial derivatives of rhs o SODE
           with respect to a list of parameters.
           _parnames is a list of parameter names.
        
           IMPORTANT: sympy module must be installed!"""

        check, msg = self.checkRates()
        if not check:
            raise BadRateError(msg)
        try:
            import sympy
        except:
            print 'ERROR: sympy module must be installed to generate partial derivative strings'
            raise
        _dxdtstrings = self.dXdt_strings()
        _symbs = {}
        for x in self.__variables:
            _symbs[x.name] = sympy.Symbol(str(x.name))
        for p in self.__parameters:
            _symbs[p.name] = sympy.Symbol(str(p.name))
        for f in self.__forcing:
            _symbs[f.name] = sympy.Symbol(str(f.name))
        _nvars = len(self.__variables)
        _npars = len(_parnames)
        _jfuncs = []
        for _i in range(_nvars):
            _jfuncs.append([])
            _ids = identifiersInExpr(_dxdtstrings[_i][1])
            if len(_ids) == 0:
                for _j in range(_npars):
                    _jfuncs[_i].append('0.0')
            else:
            
                for _j in range(_npars):
                    _res = eval(_dxdtstrings[_i][1], None, _symbs)
                    _res = _res * _scale
                    if not _symbs.has_key(_parnames[_j]):
                        _dres = '0.0'
                    else:
                        _dres = str(sympy.diff(_res, _symbs[_parnames[_j]]))
                    if _dres == '0':
                        _dres = '0.0'
                    _jfuncs[_i].append(_dres)
                    #~ print _jfuncs
        return _jfuncs
        
    def rates_func(self, with_uncertain = False, scale = 1.0, t0=0.0):
        """Generate function to compute rate vector for this model.
        
           Function has signature f(variables, t)"""

        check, msg = self.checkRates()
        if not check:
            raise BadRateError(msg)

        #compile rate laws
        ratebytecode = [compile(self.rateCalcString(v.rate, with_uncertain=with_uncertain), '<string>','eval') for v in self.__reactions]
        # create array to hold v's
        v = empty(len(self.reactions))
        en = list(enumerate(ratebytecode))
        # create array to hold forcing functions
        forcing_function = empty(len(self.forcing))
        enf = list(enumerate(self.forcing))
            
        def f(variables, t):
            t = t*scale + t0
            for i,fi in enf:
                forcing_function[i] = fi.func(variables,t)
            for i,r in en:
                v[i] = eval(r)
            return v
        def f2(variables, t):
            m_Parameters = self.__m_Parameters
            t = t*scale + t0
            for i,fi in enf:
                forcing_function[i] = fi.func(variables,t)
            for i,r in en:
                v[i] = eval(r)
            return v

        if with_uncertain:
            return f2
        else:
            return f

    def transf_func(self, with_uncertain = False, scale = 1.0, t0=0.0):
        """Generate function to compute transformations for this model.
        
           Function has signature f(variables, t)"""

        check, msg = self.checkRates()
        if not check:
            raise BadRateError(msg)

        #compile rate laws
        transfbytecode = [compile(self.rateCalcString(v.rate, with_uncertain=with_uncertain), '<string>','eval') for v in self.__transf]
        # create array to hold v's
        m_Transformations = empty(len(self.transf))
        # create array to hold forcing functions
        forcing_function = empty(len(self.forcing))
        enf = list(enumerate(self.forcing))
            
        def f(variables, t):
            t = t*scale + t0
            for i,fi in enf:
                forcing_function[i] = fi.func(variables,t)
            for i,r in enumerate(transfbytecode):
                m_Transformations[i] = eval(r)
            return m_Transformations
        def f2(variables, t):
            m_Parameters = self.__m_Parameters
            t = t*scale + t0
            for i,fi in enf:
                forcing_function[i] = fi.func(variables,t)
            for i,r in enumerate(transfbytecode):
                m_Transformations[i] = eval(r)
            return m_Transformations

        if with_uncertain:
            return f2
        else:
            return f

    def genTransformationFunction(self, f):
        if not callable(f):
            if not isinstance(f,str):
                raise TypeError('argument must be a string or a callable.')
            argnames = f.split()
            names = argnames
            nargs = len(argnames)
        else:
            cc = f.func_code
            nargs = cc.co_argcount
            argnames = cc.co_varnames[:nargs]
            names = list(argnames[:])
            if hasattr(f, 'names'):
                names[:len(f.names)] = f.names
        data = []
        for a in argnames:
            i, collection = self.findComponent(a)
            if collection == 'parameters':
                data.append(('p', getattr(self, a)))
            elif collection == 'variables':
                data.append(('v', i))
            elif collection == 'transf':
                data.append(('t', compile(self.rateCalcString(self.transf[i].rate), '<string>','eval')))
            elif collection == 'reactions':
                data.append(('r', compile(self.rateCalcString(self.reactions[i].rate), '<string>','eval')))
            elif collection == 'forcing':
                data.append(('f', self.forcing[i].func))
            else:
                raise AttributeError(a + ' is not a component in this model')
        args = [0.0]*nargs
        en = [(i,c,d) for (i, (c,d)) in enumerate(data)]
        for i,c,d in en:
            if c == 'p':
                args[i] =d

        def retf(variables, t):
            for i,c,d in en:
                if c == 'p':
                    continue
                elif c == 'v':
                    args[i] = variables[d]
                elif c == 'f':
                    args[i] = d(variables,t)
                else:
                    args[i] = eval(d)
            return f(*args)
        def retargs(variables, t):
            for i,c,d in en:
                if c == 'p':
                    continue
                elif c == 'v':
                    args[i] = variables[d]
                elif c == 'f':
                    args[i] = d(variables,t)
                else:
                    args[i] = eval(d)
            return args
                
        if callable(f):
            result = retf
            result.names = names
            return result
        else:
            result = retargs
            result.names = names
            return result

    def set_uncertain(self, uncertainparameters):
        self.__m_Parameters = uncertainparameters
    
    def getdXdt(self, with_uncertain = False, scale = 1.0, t0=0.0):
        """Generate function to compute rhs of SODE for this model.
        
           Function has signature f(variables, t)
           This is compatible with scipy.integrate.odeint"""

        check, msg = self.checkRates()
        if not check:
            print "vars = "
            print [x.name for x in self.variables]
            raise BadRateError(msg)
        #compile rate laws
        #m_Parameters = self.__m_Parameters
        ratebytecode = [compile(self.rateCalcString(v.rate, with_uncertain = with_uncertain), '<string>','eval') for v in self.__reactions]
        # compute stoichiometry matrix, scale and transpose
        N  = self.genStoichiometryMatrix()
        N *= scale
        NT = N.transpose()
        # create array to hold v's
        v = empty(len(self.reactions))
        en = list(enumerate(ratebytecode))
        # create array to hold forcing functions
        nforce = len(self.forcing)
        forcing_function = empty(nforce)
        enf = list(enumerate(self.forcing))
        
        
        def f2(variables, t):
            m_Parameters = self.__m_Parameters
            t = t*scale + t0
            for i,fi in enf:
                forcing_function[i] = fi.func(variables,t)
            for i,r in en:
                v[i] = eval(r)
            return dot(v,NT)
        def f3(variables, t):
            m_Parameters = self.__m_Parameters
            t = t*scale + t0
            for i,r in en:
                v[i] = eval(r)
            return dot(v,NT)

        if nforce >0:
            return f2
        else:
            return f3
    
    def getdXdt2(self, with_uncertain = False, scale = 1.0, t0=0.0):
        """Generate function to compute rhs of SODE for this model.
        
           Function has signature f(variables, t)
           This is compatible with scipy.integrate.odeint"""

        check, msg = self.checkRates()
        if not check:
            print "vars = "
            print [x.name for x in self.variables]
            raise BadRateError(msg)
        
        x = empty(len(self.variables))
        nforce = len(self.forcing)
        forcing_function = empty(nforce)
        m_Parameters = self.__m_Parameters
       
        ODEstr = "def f4(variables,t):\n"
        #ODEstr += '\tm_Parameters = self.__m_Parameters\n'
        if nforce >0:
            ODEstr += '\tt = t*scale + t0\n'
            for i in range(nforce):
                ODEstr += '\tforcing_function[%d] = self.forcing[%d].func(variables,t)\n' % (i,i)

        for i, xstr in enumerate(self.dXdt_strings()):
            ODEstr += "\tx[%d] = "% i + self.rateCalcString(xstr[1], with_uncertain = with_uncertain)+'\n'
        #ODEstr += "\tself.fcount+=1\n"            
        #ODEstr += "\tprint self.fcount\n"            
        ODEstr += "\treturn x"            
        
        #print ODEstr
        dct = locals()
        exec ODEstr in dct#{'self':self, 'nforce':nforce, 'x':x, 'forcing_function': forcing_function}
        return f4
        
    def dXdt_with(self, uncertainparameters, scale = 1.0, t0=0.0):
        """Generate function to compute rhs of SODE for this model.
        
           Function has signature f(variables, t)
           This is compatible with scipy.integrate.odeint"""

        check, msg = self.checkRates()
        if not check:
            raise BadRateError(msg)
        #compile rate laws
        ratebytecode = [compile(self.rateCalcString(v.rate, with_uncertain = True), '<string>','eval') for v in self.__reactions]
        # compute stoichiometry matrix, scale and transpose
        N  = self.genStoichiometryMatrix()
        N *= scale
        NT = N.transpose()

        # create array to hold v's
        v = empty(len(self.reactions))
        en = list(enumerate(ratebytecode))
        # create array to hold forcing functions
        forcing_function = empty(len(self.forcing))
        enf = list(enumerate(self.forcing))

        def f(variables, t):
            m_Parameters = uncertainparameters
            t = t*scale + t0
            for i,fi in enf:
                forcing_function[i] = fi.func(variables,t)
            for i,r in en:
                v[i] = eval(r)
            return dot(v,NT)
        return f

    def getJacobian(self, with_uncertain = False, scale = 1.0, t0=0.0):
        """Generate function to compute the jacobian for this model.
        
           Function has signature J(variables, t)
           and returns an nvars x nvars numpy array
           IMPORTANT: sympy module must be installed!"""

        Jstrings = self.Jacobian_strings(_scale = scale)
        nvars = len(Jstrings)
        
        #compile rate laws
        ratebytecode = [[compile(self.rateCalcString(col, with_uncertain = with_uncertain), '<string>','eval') for col in line] for line in Jstrings]

        forcing_function = empty(len(self.forcing))
        enf = list(enumerate(self.forcing))
        
        def J(variables, t):
            Jarray = empty((nvars,nvars), float)
            t = t*scale + t0
            for i,fi in enf:
                forcing_function[i] = fi.func(variables,t)
            for i in range(nvars):
                for j in range(nvars):
                    Jarray[i,j] = eval(ratebytecode[i][j])
            return Jarray
        def J2(variables, t):
            m_Parameters = self.__m_Parameters
            Jarray = empty((nvars,nvars), float)
            t = t*scale + t0
            for i,fi in enf:
                forcing_function[i] = fi.func(variables,t)
            for i in range(nvars):
                for j in range(nvars):
                    Jarray[i,j] = eval(ratebytecode[i][j])
            return Jarray
        if with_uncertain:
            return J2
        else:
            return J


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
    def __getForcing(self):
        return self.__forcing
    
    def __getUncertainValues(self):
        return list(getUncertainties(self)) #generator of uncertain values

    variables    = property(__getVariables)
    extvariables = property(__getExtVariables)
    reactions    = property(__getReactions)
    parameters   = property(__getParameters)
    uncertain    = property(__getUncertainValues)
    transf       = property(__getTransformations)
    forcing      = property(__getForcing)
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

def test():
    
    def force(A, t):
        return t*A
    m = Model('My first model')
    m.v1 = "A+B -> C", 3
    m.v2 = react("    -> A"  , rate = math.sqrt(4.0)/2)
    m.v3 = react("C   ->  "  , "V3 * C / (Km3 + C)")
    m.v4 = "B   ->  "  , "2*input1"
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
    m.input1 = forcing(force)
    m.at = 1.0
    m.top = 2.0
    m.input2 = transf("4*step(t,at,top)")
    
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
    print 'iterating m.reactions'
    for v in m.reactions:
        print v.name, ':', v.rate, '|', v.reagents, '->', v.products
    print '\niterating m.transf'
    for v in m.transf:
        print v.name, ':', v.rate
    print '\niterating m.variables'
    for x in m.variables:
        print x.name
    print '\niterating m.extvariables'
    for x in m.extvariables:
        print x.name
    print '\niterating m.parameters'
    for p in m.parameters:
        print p.name , '=',  p, 'bounds=', p.bounds
    print '\niterating m.uncertain'
    for x in m.uncertain:
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
    m.Km3.uncertainty(None)

    print '********** Testing rate and dXdt strings *******************'
    print 'rates_strings(): -------------------------'
    for v in m.rates_strings():
        print v
    print '\ndXdt_strings(): -------------------------'
    for dxdt in m.dXdt_strings():
        print dxdt
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
    print 'calcstring for v5:\n', m.rateCalcString(m.v5.rate)
    print
    print 'calcstring for v3 with uncertain parameters:\n', m.rateCalcString(m.v3.rate, True)
    print
    print 'calcstring for t1:\n', m.rateCalcString(m.t1.rate)
    print
    print 'calcstring for t2:\n', m.rateCalcString(m.t2.rate)
    print
    print 'calcstring for input2:\n', m.rateCalcString(m.input2.rate)
    print

    print '********** Testing rate and dXdt generating functions ******'
    print 'Operating point:'
    varvalues = [1.0, 1.0, 1.0]
    pars      = [1.0]
    t         = 2.0
    
    dxdtstrs = [b for (a,b) in m.dXdt_strings()]

    print "t =", t
    print 'variables  ='
    pprint.pprint(dict((v.name, value) for v,value in zip(m.variables, varvalues)))
    print 'parameters ='
    pprint.pprint(dict((p.name, p)     for p in m.parameters))
 
    print '---- rates using Model.rates_func() -------------------------'
    vratesfunc = m.rates_func()
    vrates = vratesfunc(varvalues,t)
    for v,r in zip(m.reactions, vrates):
        print "%s = %-20s = %s" % (v.name, v.rate, r)

    print '---- transformations using Model.transf_func() --------------'
    tratesfunc = m.transf_func()
    trates = tratesfunc(varvalues,t)
    for v,r in zip(m.transf, trates):
        print "%s = %-20s = %s" % (v.name, v.rate, r)

    print '---- dXdt using Model.dXdt() --------------------------------'
    f = m.getdXdt()
    #~ dXdt = m.dXdt(varvalues,0)
    dXdt = f(varvalues,t)
    for x,s,r in zip(m.variables, dxdtstrs, dXdt):
        print "d%s/dt = %s = %s" % (x.name, s,r)

##     print '---- dXdt using Model.dXdt2() --------------------------------'
##     print 'f = m.getdXdt2()'
##     f = m.getdXdt2()
##     dXdt = f(varvalues,t)
##     for x,s,r in zip(m.variables, dxdtstrs, dXdt):
##         print "d%s/dt = %s = %s" % (x.name, s,r)
    
##     def t1(ntimes):
##         f = m.getdXdt()
##         for i in range(ntimes):
##             dXdt = f(varvalues,t)
##         return dXdt
##     def t2(ntimes):
##         f = m.getdXdt2()
##         for i in range(ntimes):
##             dXdt = f(varvalues,t)
##         return dXdt

    print '---- dXdt using Model.dXdt() setting uncertain parameters ---'
    print 'f = m.getdXdt(with_uncertain = True)'
    f = m.getdXdt(with_uncertain = True)
    print 'setting uncertain as', dict((v.name, value) for v,value in zip(m.uncertain, pars))
    print 'm.set_uncertain(pars)'
    m.set_uncertain(pars)
    dXdt = f(varvalues,t)
    for x,s,r in zip(m.variables, dxdtstrs, dXdt):
        print "d%s/dt = %s = %s" % (x.name, s,r)

    print '---- dXdt using Model.dXdt_with(pars) ------------------------'
    print 'f = m.dXdt_with(pars)'
    f = m.dXdt_with(pars)
    dXdt   = f(varvalues,t)
    for x,s,r in zip(m.variables, dxdtstrs, dXdt):
        print "d%s/dt = %s = %s" % (x.name, s,r)

    print '---- dXdt using Model.dXdt() with a state argument (m.init) --'
    print 'm.init:', m.init
    print 'making m.V3 = 1.0'
    m.V3 = 1.0
    print 'm.V3 :', m.V3
    print
    print 'f = m.dXdt'
    f = m.dXdt
    print 'dXdt = f(m.vectorize("init"),2.0)'
    dXdt = f(m.vectorize("init"),t)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)
    print '---- same, changing state argument ---------------------------'
    m.init.A = 2.0
    print 'after m.init.A = 2.0'
    print 'm.init:', m.init
    print
    print 'f = m.dXdt'
    f = m.dXdt
    print 'dXdt = f(m.vectorize("init"),2.0)'
    dXdt = f(m.vectorize("init"),t)
    for x,r in zip(m.variables, dXdt):
        print "d%s/dt = %s" % (x.name, r)

    #~ print '\n\n\nProfiling dXdt...'
    #~ ntimes = 100000
    #~ import cProfile, pstats, StringIO
    #~ prof = cProfile.Profile()
    #~ prof = prof.runctx("t1(ntimes)", globals(), locals())
    #~ print '------------- for dXdt --------------------'
    #~ stream = StringIO.StringIO()
    #~ stats = pstats.Stats(prof, stream=stream)
    #~ stats.sort_stats("cumulative")  # Or cumulative
    #~ stats.print_stats(40)  # 40 = how many to print
    #~ # The rest is optional.
    #~ #stats.print_callees()
    #~ #stats.print_callers()
    #~ print stream.getvalue()
    #~ #logging.info("Profile data:\n%s", stream.getvalue())
    #~ prof = cProfile.Profile()
    #~ prof = prof.runctx("t2(ntimes)", globals(), locals())
    #~ print '\n-------------- for dXdt2 -------------------'
    #~ stream = StringIO.StringIO()
    #~ stats = pstats.Stats(prof, stream=stream)
    #~ stats.sort_stats("cumulative")  # Or cumulative
    #~ stats.print_stats(40)  # 40 = how many to print
    #~ # The rest is optional.
    #~ #stats.print_callees()
    #~ #stats.print_callers()
    #~ print stream.getvalue()
    #~ #logging.info("Profile data:\n%s", stream.getvalue())



if __name__ == "__main__":
    test()
 