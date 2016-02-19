"""Model class and supporting functions.

This module defines the Model class, used to hold the structure and metadata
of a kinetic model.

"""
import re
import math
from kinetics import *
import dynamics
import estimation

# ----------------------------------------------------------------------------
#         Functions to check the validity of math expressions
# ----------------------------------------------------------------------------
__globs = {}
__haskinetics = {}
for k, v in globals().items():
    if hasattr(v, "is_rate"):
        __haskinetics[k] = v
__globs.update(__haskinetics)

v = vars(math)
for k in v:
    if isinstance(v[k], float) or callable(v[k]):
        __globs[k] = v[k]
#__globs.update(vars(math))


def register_kin_func(f):
    f.is_rate = True
    __globs[f.__name__] = f
    globals()[f.__name__] = f

# ----------------------------------------------------------------------------
#         Regular expressions for stoichiometry patterns
# ----------------------------------------------------------------------------
fracnumberpattern = r"[-]?\d*[.]?\d+"
realnumberpattern = fracnumberpattern + r"(e[-]?\d+)?"
fracnumber = re.compile(fracnumberpattern, re.IGNORECASE)
realnumber = re.compile(realnumberpattern, re.IGNORECASE)

stoichiompattern = r"^\s*(?P<reagents>.*)\s*(?P<irreversible>->|<=>)\s*(?P<products>.*)\s*$"
chemcomplexpattern = r"^\s*(?P<coef>("+realnumberpattern+")?)\s*(?P<variable>[_a-z]\w*)\s*$"

stoichiom = re.compile(stoichiompattern,    re.IGNORECASE)
chemcomplex = re.compile(chemcomplexpattern, re.IGNORECASE)

# ----------------------------------------------------------------------------
#         Utility functions
# ----------------------------------------------------------------------------


def _is_sequence(arg):
    return (not hasattr(arg, "strip") and
            hasattr(arg, "__getitem__") or
            hasattr(arg, "__iter__"))


def _is_number(a):
    return (isinstance(a, float) or
            isinstance(a, int) or
            isinstance(a, long))


def processStoich(expr):
    """Split a stoichiometry string into reagents, products and irreversible flag.

    This function accepts a string that conforms to a pattern like

    2 A + B -> 3.5 C

    and splits into reagents, products and a boolean flag for irreversibility.

    Parameters
    ----------
    expr : str
        A stoichiomety pattern.

    Returns
    -------
    tuple as (reagents, products, irreversible)
        `reagents` and `products` are lists of
        (`name`: str, `coefficient`:float)
        describing the 'complexes' of the stoichiometry.

        `irreversible` (bool) True if '->' is the separator, False if '<=>'
        is the separator.

    Raises
    ------
    BadStoichError
        If `expr` is not a properly formatted stoichiometry string.

    """

    match = stoichiom.match(expr)
    if not match:
        raise BadStoichError("Bad stoichiometry definition:\n" + expr)

    # process irreversible
    irrsign = match.group('irreversible')
    irreversible = irrsign == "->"
    reagents = []
    products = []

    # process stoichiometry
    fields = [(reagents, 'reagents'), (products, 'products')]
    for target, f in fields:
        complexesstring = match.group(f).strip()
        if len(complexesstring) == 0:  # empty complexes allowed
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
                if coef == 0.0:
                    continue  # a coef equal to zero means ignore
                target.append((var, coef))
            else:
                raise BadStoichError("Bad stoichiometry definition:\n" + expr)
    return reagents, products, irreversible


def _massActionStr(k=1.0, reagents=[]):
    res = str(float(k))
    factors = []
    for var, coef in reagents:
        if coef == 1.0:
            factor = '%s' % var
        else:
            factor = '%s**%f' % (var, coef)
        factors.append(factor)
    factors = '*'.join(factors)
    if factors != '':
        res = res + '*' + factors
    return res

# ----------------------------------------------------------------------------
#         Model and Model component classes
# ----------------------------------------------------------------------------


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
        vv = float(value)  # can raise ValueError
        return constValue(value, name=name)

    # seeking proper bounds pair
    lv = len(value)  # can raise TypeError

    # value has len...
    # must be exactely two
    if lv != 2:
        raise TypeError(value + ' is not a pair of numbers')
    vv0 = float(value[0])  # can raise ValueError
    vv1 = float(value[1])  # can raise ValueError
    return Bounds(name, vv0, vv1)


def constValue(value=None, name='?'):
    if _is_number(value):
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
            raise BadTypeComponent("Can not assign"+str(value)+"to %s.%s bounds" % (obj.name, name))
        else:
            raise BadTypeComponent("Can not assign"+str(value)+"to %s.%s" % (obj.name, name))

    c = obj.__dict__['_ownparameters']
    already_exists = name in c
    if not already_exists:
        if isinstance(vv, ConstValue):
            newvalue = vv
        else:  # Bounds object
            newvalue = constValue((float(vv.lower)+float(vv.upper))/2.0, name=name)
            newvalue.bounds = vv
    else:  # aready exists
        if isinstance(vv, ConstValue):
            newvalue = vv
            newvalue.bounds = c[name].bounds
        else:  # Bounds object
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
        for k, v in parvalues.items():
            self._ownparameters[k] = constValue(value=v, name=k)

    def _get_parameter(self, name):
        if name in self._ownparameters:
            return self._ownparameters[name]
        else:
            raise AttributeError(name + ' is not a parameter of ' + self.name)

    def getp(self, name):
        o = self._get_parameter(name)
        return o

    def setp(self, name, value):
        _setPar(self, name, value)

    def set_bounds(self, name, value):
        if value is None:
            self.reset_bounds(name)
        else:
            _setPar(self, name, value, is_bounds=True)

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
        
    @property
    def parameters(self):
        return [(p.name, p) for p in list(self._ownparameters.values())]

    def _copy_pars(self):
        ret = {}
        for k, v in self._ownparameters.items():
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
    def __init__(self, haspar_obj):
        self.__dict__['_haspar_obj'] = haspar_obj

    def __len__(self):
        return len(self._haspar_obj._ownparameters)

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        r = self._haspar_obj.getp(name)
        if r is not None:
            return r
        raise AttributeError(name + ' is not in %s' % self._haspar_obj.name)

    def __setattr__(self, name, value):
        if not name.startswith('_'):
            self._haspar_obj.setp(name, value)
        else:
            object.__setattr__(self, name, value)

    def __contains__(self, name):
        try:
            r = self._haspar_obj.getp(name)
        except AttributeError:
            return False
        return True


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
        return '(%s)' % ", ".join(['%s = %s' % (x, str(float(value))) for (x, value) in self._ownparameters.items()])


class _HasRate(_HasOwnParameters):

    def __init__(self, name='?', rate='0.0', parvalues=None):
        _HasOwnParameters.__init__(self, name=name, parvalues=parvalues)
        self.__rate = rate.strip()

    def __str__(self):
        res = "%s:\n  rate = %s\n" % (self.name, str(self()))
        if len(self._ownparameters) > 0:
            res += "  Parameters:\n"
            for k, v in self._ownparameters.items():
                res += "    %s = %g\n" % (k, v)
        return res

    def __call__(self, fully_qualified=False):
        rate = self.__rate
        if fully_qualified:
            for localparname in self._ownparameters:
                fully = '%s.%s' % (self.name, localparname)
                rate = re.sub(r"(?<!\.)\b%s\b(?![.\[])" % localparname, fully, rate)
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
        res = ['%s:' % self.name,
               '  reagents: %s' % str(self._reagents),
               '  products: %s' % str(self._products),
               '  stoichiometry: %s' % self.stoichiometry_string,
               '  rate = %s' % str(self())]
        res = '\n'.join(res) + '\n'
        
        if len(self._ownparameters) > 0:
            resp = ["  Parameters:"]
            for k, v in self._ownparameters.items():
                resp.append("    %s = %g" % (k, v))
            res = res + '\n'.join(resp) + '\n'
        return res

    @property
    def reagents(self):
        """The reagents of the reaction."""
        return self._reagents
    
    @property
    def products(self):
        """The products of the reaction."""
        return self._products
    
    @property
    def stoichiometry(self):
        """The stoichiometry of the reaction.
           
           This is just a list of (coefficient, name) pairs with
           reagents with negative coefficients"""
        res = [(v, -c) for (v, c) in self._reagents]
        res.extend([(v, c) for (v, c) in self._products])
        return res
    
    def _stoichiometry_string(self):
        """Generate a canonical string representation of stoichiometry"""
        left = []
        for (v, c) in self._reagents:
            if c == 1:
                c = ''
            elif int(c) == c:
                c = str(int(c))
            else:
                c = str(c)
            left.append('%s %s' % (c, v))
        right = []
        for (v, c) in self._products:
            if c == 1:
                c = ''
            elif int(c) == c:
                c = str(int(c))
            else:
                c = str(c)
            right.append('%s %s' % (c, v))
        left = ' + '.join(left)
        right = ' + '.join(right)
        if self._irreversible:
            irrsign = "->"
        else:
            irrsign = "<=>"
        return ('%s %s %s' % (left, irrsign, right)).strip()
    
    stoichiometry_string = property(_stoichiometry_string)

    def __eq__(self, other):
        if not _HasRate.__eq__(self, other):
            return False
        if (self._reagents != other._reagents) or (self._products != other._products) or (self._irreversible != other._irreversible):
            return False
        return True


class Transformation(_HasRate):
    def __init__(self, name, rate, parvalues=None):
        _HasRate.__init__(self, name, rate, parvalues=parvalues)


class ConstValue(float, ModelObject):

    def __new__(cls, value):
        return float.__new__(cls, value)

    def initialize(self, aname='?'):
        ModelObject.__init__(self, aname)
        self.bounds = None

    def pprint(self):
        res = float.__str__(self)
        if self.bounds:
            res += " ? (min = %f, max=%f)" % (self.bounds.min, self.bounds.max)
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
        if value is None:
            self.reset_bounds()
            return
        try:
            b = toConstOrBounds(self.name, value, is_bounds=True)
        except (TypeError, ValueError):
            msg = "Can not assign %s to %s.bounds" % (str(value), self.name)
            raise BadTypeComponent(msg)
        self.bounds = b

    def get_bounds(self):
        if self.bounds is None:
            return None
        else:
            return (self.bounds.lower, self.bounds.upper)

    def reset_bounds(self):
        self.bounds = None


class Bounds(ModelObject):
    def __init__(self, aname, lower=0.0, upper=1.0):
        ModelObject.__init__(self, name=aname)
        self.lower = lower
        self.upper = upper

    def __str__(self):
        return "(lower=%f, upper=%f)" % (self.lower, self.upper)


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

    def __contains__(self, item):
        r = self.collection.get(item)
        if r is not None:
            return True
        else:
            return None


class _Simple_Collection_Accessor(object):
    def __init__(self, model, collection):
        self.__dict__['model'] = model
        self.__dict__['collection'] = collection

    def __iter__(self):
        return iter(self.collection)

    def __len__(self):
        return len(self.collection)

    def __contains__(self, item):
        if item in self.collection:
            return True
        else:
            return None


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

    def __contains__(self, item):
        if item in self._model._init._ownparameters:
            return True
        else:
            return None


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
                    yield value.copy(new_name=v.name + '.' + iname)

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
            report = (name, self._model.name)
            raise AttributeError('%s is not a parameter of %s' % report)

    def __contains__(self, item):
        try:
            o = self._model.getp(item)
        except AttributeError:
            return False
        return True

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
                yield x.copy(new_name='init.' + iname)
        collections = [self._reactions, self._transf]
        for c in collections:
            for v in c:
                for iname, x in v._ownparameters.items():
                    if x.bounds is not None:
                        yield x.copy(new_name=v.name + '.' + iname)


class Model(ModelObject):
    """The class that holds the description of a kinetic model.

    This class holds several members describing the data associated with
    the description of a kinetic model.

    A model is comprised of:

    - reactions
    - parameters
    - initial values
    - transformations
    - external variables

    Attributes
    ----------
    reactions : Accessor
        The processes ("reactions") in the model.
    transformations : Accessor to collection
        The transformations in the model.
    parameters : Accessor
        The parameters in the model.
    varnames : list[str]
        The names of the variables defined in the model. This list should be
        treated as read-only and is internally refreshed everytime a reaction
        is added or changed.
    extvariables : list[str]
        The names of the external variables defined in the model. This list
        should be treated as read-only and is internally refreshed everytime a
        reaction or parameter is added or changed.
    init : Accessor
        The initial state of the model.
    with_bounds : Accessor
        Iterates through the parameters in the model for which Bounds were
        assigned.
    """

    def __init__(self, title=""):
        """Construct an empty model object.


        Parameters
        ----------
        title : str
            The title of the model.
        """
        
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
        self.varnames = self.__variables
        self.extvariables = self.__extvariables
        self.init = _init_Accessor(self)
        self.parameters = _Parameters_Accessor(self)
        self.with_bounds = _With_Bounds_Accessor(self)

    def _set_in_collection(self, name, col, newobj):
        for c, elem in enumerate(col):
            if elem.name == name:
                col[c] = newobj
                return
        col.append(newobj)

    def set_reaction(self, name, stoichiometry, rate=0.0, pars=None):
        """Insert or modify a reaction in the model.


        Parameters
        ----------
        name : str
            The name of the reaction.
        stoichiometry : str
            The stoichiometry of the reaction.
        rate : str or int or float.
            The kinetic function of the reaction. If it is a number,
            a mass-action rate will be assumed.
        pars : dict of iterable of (name, value) pairs
            The 'local' parameters of the reaction.
        """
        reagents, products, irrv = processStoich(stoichiometry)
        if _is_number(rate):
            rate = _massActionStr(rate, reagents)

        newobj = Reaction(name, reagents, products, rate, pars, irrv)

        self._set_in_collection(name, self.__reactions, newobj)
        self._refreshVars()

    def set_transformation(self, name, rate=0.0, pars=None):
        """Insert or modify a transformation in the model.


        Parameters
        ----------
        name : str
            The name of the transformation.
        rate : str or int or float.
            The rate of the transformation. If it is a number,
            a constant rate will be assumed.
        pars : dict of iterable of (name, value) pairs
            The 'local' parameters of the transformation.

        """

        if _is_number(rate):
            rate = str(float(rate))

        newobj = Transformation(name, rate, pars)
        self._set_in_collection(name, self.__transf, newobj)
        self._refreshVars()

    def set_variable_dXdt(self, name, rate=0.0, pars=None):
        """Insert or modify a dx/dt equation in the model.


        Parameters
        ----------
        name : str
            The name of the transformation.
        rate : str or int or float.
            The rhs of the equation. If it is a number,
            a constant rate will be assumed.
        pars : dict of iterable of (name, value) pairs
            The 'local' parameters of the equation.

        """
        if _is_number(rate):
            rate = str(float(rate))

        react_name = 'd_%s_dt' % name
        stoich = ' -> %s' % name
        name = react_name  # hope this works...
        self.set_reaction(name, stoich, rate, pars)

    def setp(self, name, value):
        """Insert or modify a parameter of the model.


        Parameters
        ----------
        name : str
            The name of the parameter. "Dot" access to parameters of reactions
            or transformations, for example ``model.setp('v1.k', 2)`` is allowed.
        value : number or str that can be transformed to a float.
            The value of the parameter.
        """
        if '.' in name:
            alist = name.split('.')
            vn, name = alist[:2]
            # find if the model has an existing  object with that name
            # start with strict types
            o = self._get_reaction_or_transf(vn)
        else:
            o = self
        _setPar(o, name, value)
        self._refreshVars()


    def getp(self, name):
        """Retrieve a parameter of the model.


        Parameters
        ----------
        name : str
            The name of the parameter. "Dot" access to parameters of reactions
            or transformations, for example ``model.setp('v1.k', 2)`` is allowed.
        Returns
        -------
        float
            The value of the parameter

        """
        if '.' in name:
            alist = name.split('.')
            vname, name = alist[:2]
            o = self._get_reaction_or_transf(vname)
            return o.getp(name)
        else:
            if name in self._ownparameters:
                return self._ownparameters[name]
            else:
                raise AttributeError(name + ' is not a parameter of ' + self.name)


    def _get_reaction_or_transf(self, name):
        o = self.__reactions.get(name)
        if o is None:
            o = self.__transf.get(name)
        if o is None:
            raise AttributeError('%s is not a component of this model' % name)
        return o

    def _findComponent(self, name):
        if name in self._ownparameters:
            return -1, 'parameters'
        o = self.__reactions.get(name)
        if o is not None:
            return o, 'reactions'
        try:
            c = self.__variables.index(name)
            return c, 'variables'
        except:
            pass
        o = self.__transf.get(name)
        if o is not None:
            return o, 'transf'
        raise AttributeError('%s is not a component in this model' % name)

    def set_bounds(self, name, value):
        if '.' in name:
            alist = name.split('.')
            vn, name = alist[:2]
            # find if the model has an existing  object with that name
            # start with strict types
            if vn == 'init':
                o = self._init
            else:
                o = self._get_reaction_or_transf(vn)
        else:
            o = self
        if value is None:
            o.reset_bounds(name)
        else:
            _setPar(o, name, value, is_bounds=True)

    def reset_bounds(self, name):
        if '.' in name:
            alist = name.split('.')
            vname, name = alist[:2]
            # find if the model has an existing  object with that name
            # start with strict types
            if vname == 'init':
                o = self._init
            else:
                o = self._get_reaction_or_transf(vname)
            o.reset_bounds(name)
        else:
            if name in self._ownparameters:
                self._ownparameters[name].bounds = None
            else:
                raise AttributeError(name + ' is not a parameter of ' + self.name)

    def get_bounds(self, name):
        if '.' in name:
            alist = name.split('.')
            vname, name = alist[:2]
            # find if the model has an existing  object with that name
            # start with strict types
            if vname == 'init':
                o = self._init
            else:
                o = self._get_reaction_or_transf(vname)
            return o.get_bounds(name)
        else:
            if name in self._ownparameters:
                bb = self._ownparameters[name].bounds
                if bb is None:
                    return None
                return (bb.lower, bb.upper)
            else:
                raise AttributeError(name + ' is not a parameter of ' + self.name)

    def set_init(self, *p, **pdict):
        dpars = dict(*p, **pdict)
        for k in dpars:
            self._init.setp(k, dpars[k])
        self._refreshVars()

    def reset_init(self):
        self._init.reset()

    def get_init(self, names=None, default=0.0):
        if names is None:
            return self._init._ownparameters
        if not _is_sequence(names):
            try:
                p = self._init.getp(names)
            except AttributeError:
                p = default
            return p
        r = {}
        for n in names:
            try:
                p = self._init.getp(n)
            except AttributeError:
                p = default
            r[n] = p
        return r

    def checkRates(self):
        self._refreshVars()
        for collection in (self.__reactions, self.__transf):
            for v in collection:
                resstring, value = _test_with_everything(v(), self, v)
                if resstring != "":
                    return False, '%s\nin rate of %s: %s' % (resstring, v.name, v())
        return True, 'OK'

    def __str__(self):
        return self.info()

    def info(self, no_check=False):
        """Generate a string with a description of the model.
        
        Used when a string describing a model is needed,
        for example in `print(model)`.


        Parameters
        ----------
        no_check : boolean
            Whether a check of the validity of reaction rates is performed.
        Returns
        -------
        str
            A string with a description of the model.

        """
        self._refreshVars()
        if not no_check:
            check, msg = self.checkRates()
            if not check:
                raise BadRateError(msg)
        
        res = [self.metadata['title']]
        res.append("\nVariables: %s\n" % " ".join(self.__variables))
        if len(self.__extvariables) > 0:
            res.append("External variables: %s\n" % " ".join(self.__extvariables))
        for collection in (self.__reactions, self.__transf):
            for i in collection:
                res.append(str(i))
        res.append('init: %s\n' % str(self._init))

        for p in self._ownparameters.values():
            res.append('%s = %g' % (p.name, p))

        if len(self.with_bounds) > 0:
            res.append("\nWith bounds:")

            for u in self.with_bounds:
                res.append('%s = ?(%g, %g)' % (u.name,
                                               u.bounds.lower,
                                               u.bounds.upper))
        for k in self.metadata:
            o = self.metadata[k]
            # skip title and empty container metadata
            if k == 'title' or (hasattr(o, '__len__') and len(o)==0):
                continue
            res.append("\n%s: %s" % (str(k), str(o)))
        return '\n'.join(res)

    def copy(self, new_title=None):
        """Retrieves a deep copy of a model.
        
        Parameters
        ----------
        new_title : str
            A new title can be provided.
        Returns
        -------
        Model
            A deep copy of a model.

        """
        m = Model(self.metadata['title'])
        for r in self.__reactions:
            m.set_reaction(r.name, r.stoichiometry_string, r(), r._ownparameters)
        for p in self._ownparameters.values():
            m.setp(p.name, p)
        for t in self.__transf:
            m.set_transformation(t.name, t(), t._ownparameters)
        m._init = self._init.copy()
        # handle uncertainties
        for u in self.with_bounds:
            m.set_bounds(u.name, (u.bounds.lower, u.bounds.upper))
        m.metadata.update(self.metadata)
        if new_title is not None:
            m.metadata['title'] = new_title
        self._refreshVars()
        return m

    def __eq__(self, other):
        return self._is_equal_to(other, verbose=False)

    def _is_equal_to(self, other, verbose=False):
        if not ModelObject.__eq__(self, other):
            return False
        self._refreshVars()
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
        for cname, c1, c2 in zip(cnames, collections1, collections2):
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
        dpars = dict(*p, **pdict)
        for k in dpars:
            self.setp(k, dpars[k])

    def solve(self, **kwargs):
        return dynamics.solve(self, **kwargs)

    def scan(self, plan, **kwargs):
        return dynamics.scan(self, plan, **kwargs)

    def estimate(self, timecourses=None, **kwargs):
        return estimation.s_timate(self, timecourses=timecourses, **kwargs)

    def set_uncertain(self, uncertainparameters):
        self.__m_Parameters = uncertainparameters

    def _genlocs4rate(self, obj=None):
        for p in self._ownparameters.items():
            yield p
        collections = [self.__reactions, self.__transf]
        for c in collections:
            for v in c:
                yield (v.name, _Has_Parameters_Accessor(v))

        if (obj is not None) and (len(obj._ownparameters) > 0):
            for p in obj._ownparameters.items():
                yield p

    def _refreshVars(self):
        # can't use self.__variables=[] Triggers __setattr__
        del(self.__variables[:])
        del(self.__extvariables[:])
        for v in self.__reactions:
            for rp in (v._reagents, v._products):
                for (vname, coef) in rp:
                    if vname in self.__variables:
                        continue
                    else:
                        if vname in self._ownparameters:
                            if vname not in self.__extvariables:
                                self.__extvariables.append(vname)
                        else:
                            self.__variables.append(vname)


def _test_with_everything(valueexpr, model, obj):
    locs = {}
    for (name, value) in model._genlocs4rate(obj):
        locs[name] = value
##     print "\nvalueexpr:", valueexpr
##     print "---locs"
##     for k in locs:
##         if isinstance(locs[k], _Has_Parameters_Accessor):
##             print k, 'has parameters and rate'
##         else:
##             print k, '=', locs[k]
##     print "---__globs"
##     for k in __globs:
##         hasattr(v, "is_rate")
##         if hasattr(__globs[k], "is_rate"):
##             print k, ' is rate'
##         else:
##             print k, '-->', type(__globs[k])

    # part 1: nonpermissive, except for NameError
    try:
        value = float(eval(valueexpr, __globs, locs))
    except NameError:
        pass
    except Exception, e:
##        print 'returned on first pass'
        return ("%s : %s" % (str(e.__class__.__name__), str(e)), 0.0)
    # part 2: permissive, with dummy values (1.0) for vars
    vardict = {}
    for i in model._Model__variables:
        vardict[i] = 1.0
    vardict['t'] = 1.0
    locs.update(vardict)
    try:
        value = float(eval(valueexpr, __globs, locs))
    except (ArithmeticError, ValueError):
        pass  # might fail but we don't know the values of vars
    except Exception, e:
##        print 'returned on second pass'
        return ("%s : %s" % (str(e.__class__.__name__), str(e)), 0.0)
    return "", value


class QueriableList(list):
    def get(self, aname):
        for o in self:
            if o.name == aname:
                return o
        return None

    def get_i(self, aname):
        for i, o in enumerate(self):
            if o.name == aname:
                return i
        return None


class BadStoichError(Exception):
    """Flags a wrong stoichiometry expression"""


class BadRateError(Exception):
    """Flags a wrong rate expression"""


class BadTypeComponent(Exception):
    """Flags an assignement of a model component to a wrong type object"""
