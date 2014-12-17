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

stoichiom   = re.compile(stoichiompattern,    re.IGNORECASE)
chemcomplex = re.compile(chemcomplexpattern, re.IGNORECASE)

identifier = re.compile(r"[_a-z]\w*", re.IGNORECASE)

def identifiersInExpr(_expr):
    iterator = identifier.finditer(_expr)
    return [_expr[m.span()[0]:m.span()[1]] for m in iterator]

#----------------------------------------------------------------------------
#         Utility functions
#----------------------------------------------------------------------------
def _is_sequence(arg):
    return (not hasattr(arg, "strip") and
            hasattr(arg, "__getitem__") or
            hasattr(arg, "__iter__"))

def processStoich(expr):
    match = stoichiom.match(expr)
    if not match:
        raise BadStoichError("Bad stoichiometry definition:\n"+ expr)

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
               raise BadStoichError("Bad stoichiometry definition:\n"+ expr)
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

#----------------------------------------------------------------------------
#         Model and Model component classes
#----------------------------------------------------------------------------

class ModelObject(object):
    """Base for all model components.
       
       The only common features are a name and a dictionary with metadata"""

    def __init__(self, name='?'):
        self.metadata = {}
        self.name = name
    
    def __eq__(self, other):
        if self.name != other.name:
            return False
        if len(self.metadata) != len(other.metadata):
            return False
        for k in self.metadata:
            if repr(self.metadata[k]) != repr(other.metadata[k]):
                return False
        return True

def toConstOrBounds(name, value, is_bounds=False):
    if not is_bounds:
        vv = float(value) # can raise ValueError
        return constValue(value, name=name)

    # seeking proper bounds pair
    lv = len(value) # can raise TypeError
    
    # value has len...
    #must be exactely two
    if lv != 2:
        raise TypeError(value + ' is not a pair of numbers')
    vv0 = float(value[0]) # can raise ValueError
    vv1 = float(value[1]) # can raise ValueError
    return Bounds(name, vv0, vv1)

def constValue(value = None, name = '?'):
    if isinstance(value, float) or isinstance(value, int) or isinstance(value, long):
        v = float(value)
        res = ConstValue(v)
        res.initialize(name)
    else:
        raise TypeError(value+' is not a float or int')
    return res

def _setPar(obj, name, value, is_bounds=False):
    try:
        vv = toConstOrBounds(name, value, is_bounds)
    except (TypeError, ValueError):
        if is_bounds:
            raise BadTypeComponent("Can not assign"+str(value)+"to %s.%s bounds"%(obj.name,name))
        else:
            raise BadTypeComponent("Can not assign"+str(value)+"to %s.%s"%(obj.name,name))

    c = obj.__dict__['_ownparameters']
    already_exists = name in c
    if not already_exists:
        if isinstance(vv, ConstValue):
            newvalue = vv
        else: #Bounds object
            newvalue = constValue((float(vv.min)+float(vv.max))/2.0, name=name)
            newvalue.bounds = vv
    else: #aready exists
        if isinstance(vv, ConstValue):
            newvalue = vv
            newvalue.bounds = c[name].bounds
        else: #Bounds object
            newvalue = constValue(c[name], name=name)
            newvalue.bounds = vv
    c[name] = newvalue

class _HasOwnParameters(ModelObject):
    def __init__(self, name='?', parvalues=None):
        ModelObject.__init__(self, name)
        self._ownparameters = {}
        if parvalues is None:
            parvalues = {}
        if not isinstance(parvalues, dict):
            parvalues = dict(parvalues)
        for k,v in parvalues.items():
            self._ownparameters[k] = constValue(value=v, name=k)

    def _get_parameter(self, name):
        if name in self._ownparameters:
            return self._ownparameters[name]
        else:
            raise AttributeError(name + ' is not a parameter of '+ self.name)
        
    def getp(self, name):
        o = self._get_parameter(name)
        return o
    
    def setp(self, name, value):
        _setPar(self, name, value)
            
    def set_bounds(self, name, value):
        if value is None:
            self.reset_bounds(name)
        else:
            _setPar(self, name, value, is_bounds = True)

    def get_bounds(self, name):
        o = self._get_parameter(name)
        if o.bounds is None:
            return None
        return (o.bounds.lower, o.bounds.upper)

    def reset_bounds(self, name):
        o = self._get_parameter(name)
        bb = o.bounds = None
    
    def __iter__(self):
        return iter(self._ownparameters.itervalues())
    
    def _copy_pars(self):
        ret = {}
        for k,v in self._ownparameters.items():
            ret[k] = constValue(value=v, name=k)
        return ret
  
    def __eq__(self, other):
        if not ModelObject.__eq__(self, other):
            return False
        if len(self._ownparameters) != len(other._ownparameters):
            return False
        for k in self._ownparameters:
            if not (self._ownparameters[k]) == (other._ownparameters[k]):
                return False
        return True


class _Has_Parameters_Accessor(object):
    def __init__(self, collection):
        self.__dict__['_collection'] = collection
        
    def __len__(self):
        return len(self._collection._ownparameters)

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        r = self._collection.getp(name)
        if r is not None:
            return r
        raise AttributeError(name + ' is not in %s'% collection.name)

    def __setattr__(self, name, value):
        if not name.startswith('_'):
            self._collection.setp(name, value)
        else:
            object.__setattr__(self, name, value)
        
class StateArray(_HasOwnParameters):
    def __init__(self, name, varvalues):
        _HasOwnParameters.__init__(self, name, varvalues)
    
    def reset(self):
        for k in self._ownparameters:
            self.setp(k, 0.0) 
        for k in self._ownparameters:
            self.reset_bounds(k) 
    
    def copy(self):
        new_state = StateArray(self.name, {})
        new_state._ownparameters = self._copy_pars()
        for k in self._ownparameters:
            new_state.set_bounds(k, self.get_bounds(k)) 
        return new_state
    
    def __str__(self):
        return '(%s)' % ", ".join(['%s = %s'% (x,str(float(value))) for (x,value) in  self._ownparameters.items()])

class _HasRate(_HasOwnParameters):
    
    def __init__(self, name='?', rate='0.0', parvalues=None):
        _HasOwnParameters.__init__(self, name=name, parvalues=parvalues)
        self.__rate = rate.strip()
    
    def __str__(self):
        res =  "%s:\n  rate = %s\n" % (self.name, str(self()))
        if len(self._ownparameters) > 0:
            res += "  Parameters:\n"
            for k, v in self._ownparameters.items():
                res += "    %s = %g\n"% (k, v)
        return res
    
    def __call__(self, fully_qualified = False):
        rate = self.__rate
        if fully_qualified:
            for localparname in self._ownparameters:
                fully = '%s.%s'%(self.name, localparname)
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
    
    def __init__(self, name, reagents, products, rate, 
                       parvalues=None, 
                       irreversible=False):
        
        _HasRate.__init__(self, name, rate, parvalues=parvalues)
        self._reagents = reagents
        self._products = products
        self._irreversible = irreversible
    
    def __str__(self):
        res = "%s:\n  reagents: %s\n  products: %s\n  stoichiometry: %s\n  rate    = %s\n" % (self.name, str(self._reagents), str(self._products), self.stoichiometry_string() , str(self())) 
        if len(self._ownparameters) > 0:
            res += "  Parameters:\n"
            for k, v in self._ownparameters.items():
                res += "    %s = %g\n"% (k, v)
        return res
    
    def stoichiometry_string(self):
        left = []
        for (v,c) in self._reagents:
            if c == 1:
                c = ''
            elif int(c) == c:
                c = str(int(c))
            else:
                c = str(c)
            left.append('%s %s'%(c, v))
        right = []
        for (v,c) in self._products:
            if c == 1:
                c = ''
            elif int(c) == c:
                c = str(int(c))
            else:
                c = str(c)
            right.append('%s %s'%(c, v))
        left = ' + '.join(left)
        right = ' + '.join(right)
        if self._irreversible:
            irrsign = "->"
        else:
            irrsign = "<=>"
        return ('%s %s %s' % (left, irrsign, right)).strip()
            
    def __eq__(self, other):
        if not _HasRate.__eq__(self, other):
            return False
        if (self._reagents != other._reagents) or (self._products != other._products) or (self._irreversible != other._irreversible):
            return False
        return True


class Variable_dXdt(_HasRate):
    def __init__(self, name, rate, parvalues=None):
        _HasRate.__init__(self, name, rate, parvalues = parvalues)


class Transformation(_HasRate):
    def __init__(self, name, rate, parvalues=None):
        _HasRate.__init__(self, name, rate, parvalues = parvalues)


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
        self.bounds = Bounds(get_name(self), float(pars[0]), float(pars[1]))
    
    def pprint(self):
        res = float.__str__(self)
        if self.bounds:
            res+= " ? (min = %f, max=%f)" % (self.bounds.min, self.bounds.max)
        return res
    
    def copy(self, new_name=None):
        r = constValue(self, self.name)
        if new_name is not None:
            r.name = new_name
        if self.bounds:
            r.bounds = Bounds(self.name, self.bounds.lower, self.bounds.upper)
            if new_name is not None:
                r.bounds.name = new_name
        return r
    
    def __eq__(self, other):
        if repr(self) != repr(other):
            return False
        if isinstance(other, ConstValue):
            sbounds = self.bounds is not None
            obounds = other.bounds is not None
            if sbounds != obounds:
                return False
            if self.bounds is not None:
                if (self.bounds.lower != other.bounds.lower) or (self.bounds.upper != other.bounds.upper):
                    return False
        return True
    
    def set_bounds(self, value):
        try:
            b = toConstOrBounds(self.name, value, is_bounds=True)
        except (TypeError, ValueError):
            raise BadTypeComponent("Can not assign"+str(value)+"to %s bounds"%(self.name))
        self.bounds = b
    
    def get_bounds(self):
        if self.bounds is None:
            return None
        else:
            return (self.bounds.lower, self.bounds.upper)
    
    def reset_bounds(self):
        self.bounds = None

class Bounds(ModelObject):
    def __init__(self, aname, lower = 0.0, upper = 1.0):
        ModelObject.__init__(self, name=aname)
        self.lower = lower
        self.upper = upper
    def __str__(self):
        return "(lower=%f, upper=%f)" % (self.lower, self.upper)

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
            rate = value[1]
            for numtype in (float,int,long):
                if isinstance(rate, numtype):
                    rate = massActionStr(rate, res[0])
            return Reaction(res[0], res[1], rate, {}, res[2])
    return value

class _Collection_Accessor(object):
    def __init__(self, model, collection):
        self.__dict__['model'] = model
        self.__dict__['collection'] = collection
        
    def __iter__(self):
        return iter(self.collection)
    def __len__(self):
        return len(self.collection)

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        r = self.collection.get(name)
        if r is not None:
            return r
        raise AttributeError(name + ' is not in this model')
        

class _init_Accessor(object):
    def __init__(self, model):
        self._model = model
        
    def __iter__(self):
        return self._model._init.__iter__()
    def __len__(self):
        return len(self._model._init._ownparameters)

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        return self._model._init.getp(name)
    
    def __setattr__(self, name, value):
        if not name.startswith('_'):
            self._model._init.setp(name, value)
        else:
            object.__setattr__(self, name, value)

class _Parameters_Accessor(object):
    def __init__(self, model):
        self._model = model
        self._reactions = model._Model__reactions
        self._transf = model._Model__transf
            
    def _get_iparameters(self):
        for p in self._model._ownparameters.values():
            yield p
        collections = [self._reactions, self._transf]
        for c in collections:
            for v in c:
                for iname, value in v._ownparameters.items():
                    yield value.copy(new_name = v.name + '.' + iname)

    def __iter__(self):
        return self._get_iparameters()
    
    def __len__(self):
        return len(list(self._get_iparameters()))
    
    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        o = self._reactions.get(name)
        if o:
            return _Has_Parameters_Accessor(o)
        o = self._transf.get(name)
        if o:
            return _Has_Parameters_Accessor(o)
        if name in self._model._ownparameters:
            return self._model.getp(name)
        else:
            raise AttributeError(name + ' is not a parameter of '+ self._model.name)
    
    def __setattr__(self, name, value):
        if not name.startswith('_'):
            self._model.setp(name, value)
        else:
            object.__setattr__(self, name, value)

class _With_Bounds_Accessor(_Parameters_Accessor):
    def __init__(self, model):
        _Parameters_Accessor.__init__(self, model)
        
    def _get_iparameters(self):
        for p in self._model._ownparameters.values():
            if p.bounds is not None:
                yield p
        for iname, x in self._model._init._ownparameters.items():
            if x.bounds is not None:
                yield x.copy(new_name = 'init.' + iname)
        collections = [self._reactions, self._transf]
        for c in collections:
            for v in c:
                for iname, x in v._ownparameters.items():
                    if x.bounds is not None:
                        yield x.copy(new_name = v.name + '.' + iname)

class Model(ModelObject):
    def __init__(self, title = ""):
        self.__dict__['_Model__reactions']         = QueriableList()
        self.__dict__['_Model__variables']         = []
        self.__dict__['_Model__extvariables']      = []
        self.__dict__['_ownparameters']            = {}
        self.__dict__['_Model__transf']            = QueriableList()
        self._init = StateArray('init', dict())
        ModelObject.__init__(self, name=title)
        self.__dict__['_Model__m_Parameters']      = None
        self.metadata['title'] = title
        
        self.reactions = _Collection_Accessor(self, self.__reactions)
        self.transformations = _Collection_Accessor(self, self.__transf)
        self.varnames = _Collection_Accessor(self, self.__variables)
        self.extvariables = _Collection_Accessor(self, self.__extvariables)
        self.init = _init_Accessor(self)
        self.parameters = _Parameters_Accessor(self)
        self.with_bounds = _With_Bounds_Accessor(self)
    
    def _set_in_collection(self, name, col, newobj):
        for c,elem in enumerate(col):
            if elem.name == name:
                col[c] = newobj
                return
        col.append(newobj)
        
    def set_reaction(self, name, stoichiometry, rate=0.0, pars=None):
        r,p,i = processStoich(stoichiometry)
        if isinstance(rate,float) or isinstance(rate,int) or isinstance(rate,long):
            rate = massActionStr(rate, r)
        
        newobj = Reaction(name, r, p, rate, pars, i)
        
        self._set_in_collection(name, self.__reactions, newobj) 
        self._refreshVars()
        
    def set_transformation(self, name, rate=0.0, pars=None):
        if isinstance(rate,float) or isinstance(rate,int) or isinstance(rate,long):
            rate = str(float(rate))
        
        newobj = Transformation(name, rate, pars)
        self._set_in_collection(name, self.__transf, newobj) 
        self._refreshVars()
        
    def set_variable_dXdt(self, name, rate=0.0, pars=None):
        if isinstance(rate, float) or isinstance(rate, int) or isinstance(rate,long):
            rate = str(float(rate))

        react_name = 'd_%s_dt'% name
        stoich = ' -> %s'% name
        name = react_name # hope this works...
        self.set_reaction(name, stoich, rate, pars)
    
    def setp(self, name, value):
        if '.' in name:
            alist = name.split('.')
            vn, name = alist[:2]
            # find if the model has an existing  object with that name
            # start with strict types
            o = self.__reactions.get(vn)
            if o is None:
                o = self.__transf.get(vn)
            if o is None:
                raise AttributeError('%s is not a component of this model'%vn)
        else:
            o = self
        _setPar(o, name, value)
        
    def set_bounds(self, name, value):
        if '.' in name:
            alist = name.split('.')
            vn, name = alist[:2]
            # find if the model has an existing  object with that name
            # start with strict types
            if vn == 'init':
                o = self._init
            else:
                o = self.__reactions.get(vn)
                if o is None:
                    o = self.__transf.get(vn)
                if o is None:
                    raise AttributeError('%s is not a component of this model'%vn)
        else:
            o = self
        if value is None:
            o.reset_bounds(name)
        else:
            _setPar(o, name, value, is_bounds = True)

    def reset_bounds(self, name):
        if '.' in name:
            alist = name.split('.')
            vname, name = alist[:2]
            # find if the model has an existing  object with that name
            # start with strict types
            if vname == 'init':
                o = self._init
            else:
                o = self.__reactions.get(vname)
                if o is None:
                    o = self.__transf.get(vname)
                if o is None:
                    raise AttributeError('%s is not a component of this model'%vname)
            o.reset_bounds(name)
        else:
            if name in self._ownparameters:
                self._ownparameters[name].bounds = None
            else:
                raise AttributeError(name + ' is not a parameter of '+ self.name)

    def getp(self, name):
        if '.' in name:
            alist = name.split('.')
            vname, name = alist[:2]
            # find if the model has an existing  object with that name
            # start with strict types
            o = self.__reactions.get(vname)
            if o is None:
                o = self.__transf.get(vname)
            if o is None:
                raise AttributeError('%s is not a component of this model'%vname)
            return o.getp(name)
        else:
            if name in self._ownparameters:
                return self._ownparameters[name]
            else:
                raise AttributeError(name + ' is not a parameter of '+ self.name)

    def get_bounds(self, name):
        if '.' in name:
            alist = name.split('.')
            vname, name = alist[:2]
            # find if the model has an existing  object with that name
            # start with strict types
            if vname == 'init':
                o = self._init
            else:
                o = self.__reactions.get(vname)
                if o is None:
                    o = self.__transf.get(vname)
                if o is None:
                    raise AttributeError('%s is not a component of this model'%vname)
            return o.get_bounds(name)
        else:
            if name in self._ownparameters:
                bb = self._ownparameters[name].bounds
                if bb is None:
                    return None
                return (bb.lower, bb.upper)
            else:
                raise AttributeError(name + ' is not a parameter of '+ self.name)

    def set_init(self, **varvalues):
        self._init = StateArray('init', varvalues)
    
    def reset_init(self):
        self._init.reset()
    
    def get_init(self, names=None):
        if names is None:
            return self._init._ownparameters
        if not _is_sequence(names):
            return self._init.getp(names)
        r = {}
        for n in names:
            r[n] = self._init.getp(n)
        return r

    def checkRates(self):
        for collection in (self.__reactions, self.__transf):
            for v in collection:
                resstring, value = _test_with_everything(v(),self, v)
                if resstring != "":
                    return False, '%s\nin rate of %s: %s' % (resstring, v.name, v())
        return True, 'OK'

    def __str__(self):
        return self.info()

    def info(self, no_check=False):
        if not no_check:
            check, msg = self.checkRates()
            if not check:
                raise BadRateError(msg)
        res = "%s\n"% self.metadata['title']
        res += "\nVariables: %s\n" % " ".join(self.__variables)
        if len(self.__extvariables) > 0:
            res += "External variables: %s\n" % " ".join(self.__extvariables)
        for collection in (self.__reactions, self.__transf):
            for i in collection:
                res += str(i)
        res += 'init: '+ str(self._init) + '\n'
        for p in self.parameters:
            res += p.name +' = '+ str(p) + '\n'
        for u in self.with_bounds:
            res += u.name + ' = ? (' + str(u.bounds.lower) + ', ' + str(u.bounds.upper) + ')\n'
        for k in self.metadata:
            res += "%s: %s\n"%(str(k), str(self.metadata[k]))
        return res
    
    def copy(self, new_title = None):
        m = Model(self.metadata['title'])
        for r in self.__reactions:
            m.set_reaction(r.name, r.stoichiometry_string(), r(), r._ownparameters)
        for p in self._ownparameters.values():
            m.setp(p.name, p)
        for t in self.__transf:
            m.set_transformation(t.name, t(), t._ownparameters)
        m._init = self._init.copy()
        #handle uncertainties
        for u in self.with_bounds:
            m.set_bounds(u.name, (u.bounds.lower,u.bounds.upper))
        m.metadata.update(self.metadata)
        if new_title is not None:
            m.metadata['title'] = new_title
        return m
    
    def __eq__(self, other):
        return self._is_equal_to(other, verbose=False)
    
    def _is_equal_to(self, other, verbose=False):
        if not ModelObject.__eq__(self, other):
            return False
        cnames = ('reactions', 'transf', 'init', 'pars', 'vars', 'extvars')
        collections1 = [self.__reactions, 
                        self.__transf, 
                        self._init._ownparameters, 
                        self._ownparameters, 
                        self.__variables, 
                        self.__extvariables]
        collections2 = [other.__reactions, 
                        other.__transf, 
                        other._init._ownparameters, 
                        other._ownparameters, 
                        other.__variables, 
                        other.__extvariables]
        for cname, c1,c2 in zip(cnames, collections1, collections2):
            if verbose:
                print
                print cname
            if len(c1) != len(c2):
                return False
            if isinstance(c1, dict):
                names = c1.keys()
            else:
                names = [v for v in c1]
            for ivname, vname in enumerate(names):
                if isinstance(vname, ModelObject):
                    vname = vname.name
                if hasattr(c1, 'get'):
                    r = c1.get(vname)
                    ro = c2.get(vname)
                else:
                    r = c1[ivname]
                    ro = c2[ivname]
                if not ro == r:
                    if verbose:
                        print vname, 'are not equal'
                    return False
                if verbose:
                    print vname, 'are equal'
        return True        
    
    def update(self, *p, **pdict):
        dpars = dict(*p)
        dpars.update(pdict)
        for k in dpars:
            self.setp(k,dpars[k])
        
    def set_uncertain(self, uncertainparameters):
        self.__m_Parameters = uncertainparameters

    def _genlocs4rate(self, obj = None):
        for p in self._ownparameters.items():
            yield p
        collections = [self.__reactions, self.__transf]
        for c in collections:
            for v in c:
                yield (v.name, v)
        
        if (obj is not None) and (len(obj._ownparameters) > 0):
            for p in obj._ownparameters.items():
                yield p

    def _refreshVars(self):
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
##     def __call__(self):
##         mq = ModelQuerier(self)
##         return mq

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

class QueriableList(list):
    def get(self, aname):
        for i in self:
            if i.name == aname:
                return i
        return None
        
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
    m.set_reaction('v1', "A+B <=> C", rate=3)
    m.set_reaction('v2', "    -> 4.5 A", rate=math.sqrt(4.0)/2)
    
    v3pars = (('V3',0.5),('Km', 4))
    m.set_reaction('v3', "C   ->  " , "V3 * C / (Km3 + C)", pars = v3pars)
    
    m.setp('B', 2.2)
    m.setp('V3', 0.6)
    m.setp('v3.Km', 4.4)
    
    m.set_reaction('v4', "B   ->  " , rate="2*4*step(t,at,top)")
    m.set_reaction('v5', "C ->",  "4.5*C*step(t,at,top)")
    
    m.set_transformation('t1', "A*Vt + C", dict(Vt=4))
    m.set_transformation('t2', "sqrt(2*A)")
    
    m.set_variable_dXdt('D', "-2 * D")
    
    m.setp('myconstant', 2 * m.getp('B') / 1.1) # should be 4.0
    m.setp('myconstant2', 2 * m.getp('v3.Km') / 1.1) # should be 8.0
    
    m.set_bounds('V3', (0.1, 1.0))
    m.set_bounds('v3.V3', [0.6, 0.9])
    
    m.setp('Km3', 4)
    m.set_bounds('Km3', (1,6))
    
    m.reset_bounds('v3.V3')
    m.set_bounds('v3.V3', m.get_bounds('V3'))
    
    m.set_init(A = 1.0, C = 1, D = 2)
    m.set_bounds('init.A', (0.8,3))
    
    m.setp('at', 1.0)
    m.setp('top', 2.0)
    
    m.set_transformation('input2', "4*step(t,at,top)")
    m.set_transformation('input3', "force(top, t)")
    
    m.metadata['where'] = 'in model'
    m.metadata['for what'] = 'testing'
    
    
    print '********** Testing model construction and printing **********'
    print '------- result of model construction:\n'
    print m
    print "!!!!  access to info as keys---------"
    print "m.metadata['for what'] =", m.metadata['for what']
    del(m.metadata['where'])
    print "\nafter del(m.metadata['where'])"
    print "m.metadata['where'] =", m.metadata.get('where')
    
    m2 = m.copy()
    print
    print '------- result of CLONING the model:\n'
    print m2
    
    print
    print '********** Testing equality of models *****************'
    print "m2 == m"
    print m2._is_equal_to(m, verbose=True)

    print
    m3 = m.copy(new_title = 'another model')
    print
    print '------- result of CLONING the model:\n'
    print m3
    
    print
    print '********** Testing equality of models *****************'
    print "m3 == m"
    print m3._is_equal_to(m)

    print
    print '********** Testing iteration of components *****************'
    print '!!! there are %d reactions in model'% len(m.reactions)
    print '---- iterating m.reactions'
    for v in m.reactions:
        print v.name, ':', v(), '|', v.stoichiometry_string()
    print '\n---- iterating m.reactions with fully qualified rates'
    for v in m.reactions:
        print v.name, ':', v(fully_qualified = True), '|', v.stoichiometry_string()
    print
    print '!!! there are %d transformations in model'% len(m.transformations)
    print '---- iterating m.transformations'
    for v in m.transformations:
        print v.name, ':', v()
    print '\n---- iterating m.transformations with fully qualified rates'
    for v in m.transformations:
        print v.name, ':', v(fully_qualified = True)
    print
    print '!!! there are %d variables in model'% len(m.varnames)
    print '---- iterating m.varnames:'
    for x in m.varnames:
        print x
    print
    print '!!! there are %d extvariables in model'% len(m.extvariables)
    print '---- iterating m.extvariables'
    for x in m.extvariables:
        print x

    print
    print '!!! there are %d initial values defined in model'% len(m.init)
    print '---- iterating m.init'
    for x in m.init:
        print x.name, '=', x, '\n  bounds=', x.bounds

    print
    print '!!! there are %d parameters in model'% len(m.parameters)
    print '---- iterating m.parameters'
    for p in m.parameters:
        print p.name , '=',  p, '\n  bounds=', p.bounds
    
    print
    print '!!! there are %d parameters with bounds in model'% len(m.with_bounds)
    print '---- iterating m.with_bounds'
    for p in m.with_bounds:
        print p.name, 'in (', p.bounds.lower, ',', p.bounds.upper, ')'
    
    print
    print '********** Testing component retrieval *********************'
    print 'm.getp("K3") :', m.getp('Km3')
    print 'm.getp("K3").name :', m.getp('Km3').name
    print 'm.getp("v3.Km") :', m.getp('v3.Km')
    print 'm.getp("v3.Km").name :', m.getp('v3.Km').name
    
    print 'm.get_init():',m.get_init()
    print 'm.get_init("A") :',m.get_init("A")
    print 'm.get_init("A").name :',m.get_init("A").name
    print 'm.get_init("A").bounds :',m.get_init("A").bounds
    print 'm.get_init(list("AC")) :',m.get_init(list("AC"))
    
    try:
        print 'm.get_init("B") :',m.get_init("B")
    except AttributeError:
        print 'm.init.getp("B") raised AttributeError'
    
##     print '\n---- iterating m.v3 (iterates parameters returning (name,value) tuples)'
##     for xname, x in m.v3:
##         print '\t', xname, '=', x
##     print
    
    print
    print '********** Testing faulty assignments *****************'
    print 'm.getp("myconstant") :',m.getp('myconstant')
    print '\n---- assigning "k9" to parameter myconstant'
    try:
        m.setp('myconstant', 'k9')
    except BadTypeComponent:
        print 'Failed! BadTypeComponent was caught.'
    print 'm.getp("myconstant") :',m.getp('myconstant'), '(still!)'
    print len(m.parameters), 'parameters total'
    
    print
    print "----- trying m.set_reaction('v1', 'A*B?C', 'hhhh')"
    try:
        m.set_reaction('v1', 'A*B?C', 'hhhh')
    except BadStoichError:
        print 'Failed! BadStoichError was caught.'
        
    print
    print '********** Testing update() function *****************'
    print '\niterating m.parameters'
    for p in m.parameters:
        print p.name , '=',  p

    print '\n---- after m.update([("V4",1.1),("V3",1.2),("Km3",1.3)])'
    m.update([("V4",1.1),("V3",1.2),("Km3",1.3)])
    for p in m.parameters:
        print p.name , '=',  p

    print '\n---- after m.update(V4=1.4, V3=1.5, Km3=1.6)'
    m.update(V4=1.4, V3=1.5, Km3=1.6)
    for p in m.parameters:
        print p.name , '=',  p

    print '\n---- after m.update([("V4",1.7)], V3=1.8, Km3=1.9)'
    m.update([("V4",1.7)], V3=1.8, Km3=1.9)
    for p in m.parameters:
        print p.name , '=',  p

    print '\n---- after dd={"V4":2.1, "V3":2.2, "Km3":2.3}; m.update(dd)'
    dd={"V4":2.1, "V3":2.2, "Km3":2.3}; m.update(dd)
    for p in m.parameters:
        print p.name , '=',  p
        
    print
    print '**********************************************'
    print '********** Testing accessors *****************'
    print '**********************************************'
    print
    print '--- m.reactions.v3'
    print m.reactions.v3
    print '--- m.reactions.v3()'
    print m.reactions.v3()
    print '--- m.reactions.v3.name'
    print m.reactions.v3.name
    print '--- m.reactions.v3.stoichiometry_string()'
    print m.reactions.v3.stoichiometry_string()
    print '----------------------------------------'
    print '--- m.transformations.t1'
    print m.transformations.t1
    print '--- m.transformations.t1()'
    print m.transformations.t1()
    print '----------------------------------------'
    print '--- m.parameters.V3'
    print m.parameters.V3
    print '--- m.parameters.V3 = 2.5'
    m.parameters.V3 = 2.5
    print '--- m.getp("V3")'
    print m.getp("V3")
    print '--- m.parameters.V3'
    print m.parameters.V3
    print '--- m.parameters.V3.get_bounds()'
    print m.parameters.V3.get_bounds()
    print '----------------------------------------'
    print '--- m.parameters.v3.Km'
    print m.parameters.v3.Km
    print '--- m.parameters.v3.Km = 2.5'
    m.parameters.v3.Km = 2.5
    print '--- m.getp("v3.Km")'
    print m.getp("v3.Km")
    print '--- m.parameters.v3.Km'
    print m.parameters.v3.Km
    print '--- m.parameters.v3.V3.get_bounds()'
    print m.parameters.v3.V3.get_bounds()
    print '--- m.parameters.v3.Km.get_bounds()'
    print m.parameters.v3.Km.get_bounds()
    print '--- m.parameters.v3.Km.set_bounds((0.3,0.7))'
    m.parameters.v3.Km.set_bounds((0.3,0.7))
    print '--- m.parameters.v3.Km.get_bounds()'
    print m.parameters.v3.Km.get_bounds()
    print '--- m.parameters.v3.Km.bounds'
    print m.parameters.v3.Km.bounds
    print '----------------------------------------'
    print '--- m.init.A'
    print m.init.A
    print '--- m.init.A = 0.9'
    m.init.A = 0.9
    print '--- m.init.A'
    print m.init.A
    print '--- m.init.E'
    try:
        print m.init.E
    except AttributeError:
        print 'raised AttributeError'

if __name__ == "__main__":
    test()
 