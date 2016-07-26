"""Model class and supporting functions.

This module defines the Model class, used to hold the structure and metadata
of a kinetic model.

"""
from __future__ import print_function, absolute_import
import re
import math
from collections import OrderedDict
from stimator.utils import _args_2_dict, _is_sequence, _is_number
import stimator.kinetics as kinetics
import stimator.dynamics as dynamics
import stimator.estimation as estimation

# ----------------------------------------------------------------------------
#         Functions to check the validity of math expressions
# ----------------------------------------------------------------------------

def get_allowed_f():
    fdict = {}

    # from module kinetics
    v = vars(kinetics)
    haskinetics = {}
    for k in v:
        if hasattr(v[k], "is_rate"):
            haskinetics[k] = v[k]
    fdict.update(haskinetics)

    # from module math
    v = vars(math)
    math_f = {}
    for k in v:
        if isinstance(v[k], float) or callable(v[k]):
            math_f[k] = v[k]
    fdict.update(math_f)
    return fdict


# ----------------------------------------------------------------------------
#         Regular expressions for stoichiometry patterns
# ----------------------------------------------------------------------------
FRAC_NUMBER_P = r"[-]?\d*[.]?\d+"
REAL_NUMBER_P = FRAC_NUMBER_P + r"(e[-]?\d+)?"
FRAC_NUMBER = re.compile(FRAC_NUMBER_P, re.IGNORECASE)
REAL_NUMBER = re.compile(REAL_NUMBER_P, re.IGNORECASE)

STOICHIOM_P = r"^\s*(?P<reagents>.*)\s*(?P<irreversible>->|<=>)\s*(?P<products>.*)\s*$"
COMPLEX_P = r"^\s*(?P<coef>("+REAL_NUMBER_P+r")?)\s*(?P<variable>[_a-z]\w*)\s*$"

STOICHIOM = re.compile(STOICHIOM_P, re.IGNORECASE)
CHEMCOMPLEX = re.compile(COMPLEX_P, re.IGNORECASE)

# ----------------------------------------------------------------------------
#         Utility functions
# ----------------------------------------------------------------------------

def process_stoich(expr):
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

    match = STOICHIOM.match(expr)
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
            m = CHEMCOMPLEX.match(c)
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


def _mass_action_str(k, reagents):
    res = str(float(k))
    factors = []
    for var, coef in reagents:
        if coef == 1.0:
            factor = '%s' % var
        else:
            factor = '%s**%f' % (var, coef)
        factors.append(factor)
    strfactors = '*'.join(factors)
    if strfactors != '':
        res = res + '*' + strfactors
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


def to_const_or_bounds(name, value, is_bounds=False):
    if value is None:
        # just return None to caller
        return None
    if not is_bounds:
        vv = float(value)  # can raise ValueError
        return create_const_value(vv, name=name)

    # seeking proper bounds pair
    lv = len(value)  # can raise TypeError

    # value has len...
    # must be exactely two
    if lv != 2:
        raise TypeError('{} is not a pair of numbers'.format(value))
    vv0 = float(value[0])  # can raise ValueError
    vv1 = float(value[1])  # can raise ValueError
    return Bounds(name, vv0, vv1)


def create_const_value(value=None, name='?', bounds=None):
    if _is_number(value):
        v = float(value)
        res = ConstValue(v, name, bounds)
    else:
        raise TypeError('{} is not a number'.format(value))
    return res


def _set_par(obj, name, value, is_bounds=False):
    try:
        vv = to_const_or_bounds(name, value, is_bounds)
    except (TypeError, ValueError):
        ms = "Can not assign {} to {}.{}"
        if is_bounds:
            raise BadTypeComponent(ms.format(value, obj.name, name) + 'bounds')
        else:
            raise BadTypeComponent(ms.format(value, obj.name, name))

    c = obj.__dict__['_ownparameters']
    already_exists = name in c

    if not already_exists:
        if vv is None:
            raise BadTypeComponent("Can not set parameter %s to None" % name)
        if isinstance(vv, ConstValue):
            newvalue = vv
        else:  # Bounds object
            nvalue = (float(vv.lower)+float(vv.upper))/2.0
            newvalue = create_const_value(nvalue, name=name, bounds=vv)
    else:  # aready exists
        if vv is None:
            if is_bounds:
                c[name].set_bounds(vv)
            else:
                del c[name]
            return
        if isinstance(vv, ConstValue):
            newvalue = vv
            newvalue.set_bounds(c[name].bounds)
        else:  # Bounds object
            newvalue = create_const_value(c[name], name=name, bounds=vv)
    c[name] = newvalue


class ConstValue(float, ModelObject):

    def __new__(cls, value, aname='?', bounds=None):
        return float.__new__(cls, value)

    def __init__(self, value, aname='?', bounds=None):
        ModelObject.__init__(self, aname)
        self.bounds = bounds

    def copy(self, new_name=None):
        name = self.name
        if new_name is not None:
            name = new_name
        r = create_const_value(self, name)
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
                if self.bounds.lower != other.bounds.lower:
                    return False
                if self.bounds.upper != other.bounds.upper:
                    return False
        return True

    def set_bounds(self, value):
        if value is None:
            self.reset_bounds()
            return
        if isinstance(value, Bounds):
            b = value
        else:
            try:
                b = to_const_or_bounds(self.name, value, is_bounds=True)
            except (TypeError, ValueError):
                msg = "Can not use %s in %s.bounds" % (str(value), self.name)
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


class _HasOwnParameters(ModelObject):
    def __init__(self, name='?', parvalues=None):
        ModelObject.__init__(self, name)
        self._ownparameters = OrderedDict()
        if parvalues is None:
            parvalues = {}
        if not isinstance(parvalues, dict):
            parvalues = dict(parvalues)
        for k, v in parvalues.items():
            self._ownparameters[k] = create_const_value(value=v, name=k)

    def _get_parameter(self, name):
        if name in self._ownparameters:
            return self._ownparameters[name]
        else:
            raise AttributeError(name + ' is not a parameter of ' + self.name)

    def getp(self, name):
        o = self._get_parameter(name)
        return o

    def setp(self, name, value):
        _set_par(self, name, value)

    def set_bounds(self, name, value):
        _set_par(self, name, value, is_bounds=True)

    def get_bounds(self, name):
        o = self._get_parameter(name)
        if o.bounds is None:
            return None
        return (o.bounds.lower, o.bounds.upper)

    def reset_bounds(self, name):
        o = self._get_parameter(name)
        o.bounds = None

    def __iter__(self):
        return iter(self._ownparameters.values())

    @property
    def parameters(self):
        return {p.name: p for p in self._ownparameters.values()}

    def _copy_pars(self):
        ret = {}
        for k, v in self._ownparameters.items():
            ret[k] = create_const_value(value=v, name=k, bounds=v.bounds)
        return ret

    def __eq__(self, other):
        if not ModelObject.__eq__(self, other):
            return False
        these_pars = self.parameters
        other_pars = other.parameters

        if len(these_pars) != len(other_pars):
            return False
        for k in these_pars:
            if not (these_pars[k]) == (other_pars.get(k)):
                return False
        return True


class _Has_Parameters_Accessor(object):
    def __init__(self, haspar_obj):
        self._haspar_obj = haspar_obj

    def __len__(self):
        return len(self._haspar_obj._ownparameters)

    def __getattr__(self, name):
        return self._haspar_obj.getp(name)

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
        tlist = [(k, str(float(v))) for (k, v) in self._ownparameters.items()]
        elist = ['{} = {}'.format(k, v) for (k, v) in tlist]
        return '({})'.format(", ".join(elist))


class _HasRate(_HasOwnParameters):

    def __init__(self, name='?', rate='0.0', parvalues=None):
        _HasOwnParameters.__init__(self, name=name, parvalues=parvalues)
        self.__rate = rate.strip()
        self._value = None

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
        rel = ['%s:' % self.name,
               '  reagents: %s' % str(self._reagents),
               '  products: %s' % str(self._products),
               '  stoichiometry: %s' % self.stoichiometry_string,
               '  rate = %s' % str(self())]
        res = '\n'.join(rel) + '\n'

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
        if self._reagents != other._reagents:
            return False
        if self._products != other._products:
            return False
        if self._irreversible != other._irreversible:
            return False
        return True


class Transformation(_HasRate):
    def __init__(self, name, rate, parvalues=None):
        _HasRate.__init__(self, name, rate, parvalues=parvalues)


class Input_Variable(_HasRate):
    def __init__(self, name, rate, parvalues=None):
        _HasRate.__init__(self, name, rate, parvalues=parvalues)


class _Collection_Accessor(object):
    def __init__(self, model, collection):
        self.__dict__['model'] = model
        self.__dict__['collection'] = collection

    def __iter__(self):
        return iter(self.collection)

    def __len__(self):
        return len(self.collection)

    def __getattr__(self, name):
##         if name in self.__dict__:
##             return self.__dict__[name]
        r = self.collection.get(name)
        if r is not None:
            return r
        raise AttributeError(name + ' is not in this model')

    def __contains__(self, item):
        r = self.collection.get(item)
        return r is not None


def _set_in_collection(name, col, newobj):
    for c, elem in enumerate(col):
        if elem.name == name:
            col[c] = newobj
            return
    col.append(newobj)


class _init_Accessor(object):
    def __init__(self, model):
        self._model = model

    def __iter__(self):
        return self._model._init.__iter__()

    def __len__(self):
        return len(self._model._init._ownparameters)

    def __getattr__(self, name):
##         if name in self.__dict__:
##             return self.__dict__[name]
        return self._model._init.getp(name)

    def __setattr__(self, name, value):
        if not name.startswith('_'):
            self._model._init.setp(name, value)
        else:
            object.__setattr__(self, name, value)

    def __contains__(self, item):
        return item in self._model._init._ownparameters


class _Parameters_Accessor(object):
    def __init__(self, model):
        self._model = model
        self._reactions = model._Model__reactions
        self._transf = model._Model__transf
        self._invars = model._Model__invars

    def _get_iparameters(self):
        for p in self._model._ownparameters.values():
            yield p
        collections = [self._reactions, self._transf, self._invars]
        for c in collections:
            for v in c:
                for iname, value in v._ownparameters.items():
                    yield value.copy(new_name=v.name + '.' + iname)

    def __iter__(self):
        return self._get_iparameters()

    def __len__(self):
        return len(list(self._get_iparameters()))

    def __getattr__(self, name):
##         if name in self.__dict__:
##             return self.__dict__[name]
        o = self._reactions.get(name)
        if o:
            return _Has_Parameters_Accessor(o)
        o = self._transf.get(name)
        if o:
            return _Has_Parameters_Accessor(o)
        o = self._invars.get(name)
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
        collections = [self._reactions, self._transf, self._invars]
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

        self.__reactions = QueriableList()
        self.__variables = []
        self.__extvariables = []
        self._ownparameters = OrderedDict()
        self.__transf = QueriableList()
        self.__invars = QueriableList()
        self._init = StateArray('init', dict())
        ModelObject.__init__(self, name=title)
        self.__m_Parameters = None
        self.metadata['title'] = title

        self.reactions = _Collection_Accessor(self, self.__reactions)
        self.transformations = _Collection_Accessor(self, self.__transf)
        self.input_variables = _Collection_Accessor(self, self.__invars)
        self.varnames = self.__variables
        self.extvariables = self.__extvariables
        self.init = _init_Accessor(self)
        self.parameters = _Parameters_Accessor(self)
        self.with_bounds = _With_Bounds_Accessor(self)
        self._usable_functions = get_allowed_f()

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
        reagents, products, irrv = process_stoich(stoichiometry)
        if _is_number(rate):
            rate = _mass_action_str(rate, reagents)

        newobj = Reaction(name, reagents, products, rate, pars, irrv)

        _set_in_collection(name, self.__reactions, newobj)
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
        _set_in_collection(name, self.__transf, newobj)
        self._refreshVars()

    def set_input_var(self, name, rate=0.0, pars=None):
        """Insert or modify an input variable in the model.


        Parameters
        ----------
        name : str
            The name of the input variable.
        rate : str or int or float.
            The rate of the input variable. If it is a number,
            a constant rate will be assumed.
        pars : dict of iterable of (name, value) pairs
            The 'local' parameters of the input variable.

        """

        if _is_number(rate):
            rate = str(float(rate))

        newobj = Input_Variable(name, rate, pars)
        _set_in_collection(name, self.__invars, newobj)
        if name in self.parameters:
            self.setp(name, None)
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

    def setp(self, *p, **pdict):
        """Insert or modify a parameter of the model.


        Parameters
        ----------
        name : str
            The name of the parameter. "Dot" access to parameters of reactions
            or transformations, for example ``model.setp('v1.k', 2)`` is allowed.
        value : number or str that can be transformed to a float.
            The value of the parameter.
        """
        dpars = _args_2_dict(*p, **pdict)

        for name, value in dpars.items():
            if '.' in name:
                alist = name.split('.')
                vn, name = alist[:2]
                # find if the model has an existing  object with that name
                # start with strict types
                o = self._get_obj_withpars(vn)
            else:
                o = self
                if value is not None and name in self.input_variables:
                    # delete name in collection self.input_variables
                    self.__invars.delete(name)
            _set_par(o, name, value)
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
            o = self._get_obj_withpars(vname)
            return o.getp(name)
        else:
            if name in self._ownparameters:
                return self._ownparameters[name]
            else:
                raise AttributeError(name + ' is not a parameter of ' + self.name)


    def _get_obj_withpars(self, name):
        o = self.__reactions.get(name)
        if o is None:
            o = self.__transf.get(name)
        if o is None:
            o = self.__invars.get(name)
        if o is None:
            raise AttributeError('%s is not a component of this model' % name)
        return o

    def set_bounds(self, name, value):
        if '.' in name:
            alist = name.split('.')
            vn, name = alist[:2]
            # find if the model has an existing  object with that name
            # start with strict types
            if vn == 'init':
                o = self._init
            else:
                o = self._get_obj_withpars(vn)
        else:
            o = self
        _set_par(o, name, value, is_bounds=True)

    def reset_bounds(self, name):
        if '.' in name:
            alist = name.split('.')
            vname, name = alist[:2]
            # find if the model has an existing  object with that name
            # start with strict types
            if vname == 'init':
                o = self._init
            else:
                o = self._get_obj_withpars(vname)
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
                o = self._get_obj_withpars(vname)
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
        dpars = _args_2_dict(*p, **pdict)
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
        for collection in (self.__reactions, self.__transf, self.__invars):
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
            if k == 'title' or (hasattr(o, '__len__') and len(o) == 0):
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
        for i in self.__invars:
            m.set_input_var(i.name, i(), i._ownparameters)
        m._init = self._init.copy()
        # handle uncertainties
        for u in self.with_bounds:
            m.set_bounds(u.name, (u.bounds.lower, u.bounds.upper))
        m.metadata.update(self.metadata)
        m._usable_functions.update(self._usable_functions)

        if new_title is not None:
            m.metadata['title'] = new_title

        self._refreshVars()
        return m

    def __eq__(self, other):
        return self._is_equal_to(other, verbose=False)

    def __ne__(self, other):
        return not self.__eq__(other)

    def _is_equal_to(self, other, verbose=False):
        if not ModelObject.__eq__(self, other):
            if verbose:
                print ('ModelObjects are not the same')
            return False
        self._refreshVars()
        cnames = ('reactions',
                  'transf', 'invars', 'init',
                  'pars', 'vars', 'extvars')
        collections1 = [self.__reactions,
                        self.__transf,
                        self.__invars,
                        self._init._ownparameters,
                        self._ownparameters,
                        self.__variables,
                        self.__extvariables]
        collections2 = [other.__reactions,
                        other.__transf,
                        other.__invars,
                        other._init._ownparameters,
                        other._ownparameters,
                        other.__variables,
                        other.__extvariables]
        for cname, c1, c2 in zip(cnames, collections1, collections2):
            if verbose:
                print ('\n', cname)
            if len(c1) != len(c2):
                if verbose:
                    print (cname, 'lenghts are not equal')
                return False
            if isinstance(c1, dict):
                names = c1.keys()
            else:
                names = list(c1)
                # names = [v for v in c1]
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
                        print (vname, 'are not equal')
                    return False
                if verbose:
                    print (vname, 'are equal')
        return True

    def solve(self, **kwargs):
        return dynamics.solve(self, **kwargs)

    def scan(self, plan, **kwargs):
        return dynamics.scan(self, plan, **kwargs)

    def estimate(self, timecourses=None, **kwargs):
        return estimation.s_timate(self, timecourses=timecourses, **kwargs)

    def set_uncertain(self, uncertainparameters):
        self.__m_Parameters = uncertainparameters

    def register_kin_func(self, f):
        f.is_rate = True
        self._usable_functions[f.__name__] = f
        #globals()[f.__name__] = f

    def _refreshVars(self):
        # can't use self.__variables=[] Triggers __setattr__
        del self.__variables[:]
        del self.__extvariables[:]
        for v in self.__reactions:
            for rp in (v._reagents, v._products):
                for (vname, _) in rp:
                    if vname in self.__variables:
                        continue
                    else:
                        if vname in self.input_variables or vname in self._ownparameters:
                            if vname not in self.__extvariables:
                                self.__extvariables.append(vname)
                        else:
                            self.__variables.append(vname)

    def checkRates(self):
        self._refreshVars()
        # Reset input variables
        for v in self.__invars:
            v._value = None
        for collection in (self.__invars, self.__reactions, self.__transf):
            for v in collection:
                msg, value = self._test_with_everything(v(), v)
                if msg != "":
                    return False, '%s\nin rate of %s: %s' % (msg, v.name, v())
                else:
                    v._value = value
        return True, 'OK'

    def _genlocs4rate(self, obj=None):
        # global model parameters
        for p in self._ownparameters.items():
            yield p

        # values of input variables
        for v in self.__invars:
            yield (v.name, v._value)

        # parameters own by reactions or transformations
        collections = [self.__reactions, self.__transf]
        for c in collections:
            for v in c:
                yield (v.name, _Has_Parameters_Accessor(v))

        # own parameters of obj
        # this may overide (correctely) other parameters with the same name
        if (obj is not None) and (len(obj._ownparameters) > 0):
            for p in obj._ownparameters.items():
                yield p


    def _test_with_everything(self, expr, obj):
        locs = dict(self._genlocs4rate(obj))

##         print '\nChecking {}, expr = {}'.format(obj.name, expr)
##         print "---locs"
##         for k in locs:
##             if k in self.input_variables:
##                 pf = '{} is a {}, value = {}'
##                 print (pf.format(k, 'Input var', locs[k]))
##             elif isinstance(locs[k], _Has_Parameters_Accessor):
##                 pf = '{} is a {}'
##                 if k in self.reactions:
##                     ttt = 'Reaction'
##                 elif k in self.transformations:
##                     ttt = 'Transformation'
##                 else:
##                     ttt = 'Something with parameters'
##                 print (pf.format(k, ttt))
##             else:
##                 print k, '=', locs[k]

##         print '\nfirst pass...'

        # part 1: nonpermissive, except for NameError
        try:
            value = float(eval(expr, self._usable_functions, locs))
        except NameError:
            pass
        except TypeError as e:
            return ("Invalid use of a rate in expression", 0.0)
        except Exception as e:
##             print 'failed on first pass'
            return ("%s : %s" % (str(e.__class__.__name__), str(e)), 0.0)
##         print 'second pass...'
        # part 2: permissive, with dummy values (1.0) for vars
        vardict = {}
        for i in self.__variables:
            vardict[i] = 1.0
        vardict['t'] = 1.0
        locs.update(vardict)
        try:
            value = float(eval(expr, self._usable_functions, locs))
        except (ArithmeticError, ValueError):
            pass  # might fail but we don't know the values of vars
        except Exception as e:
##             print 'failed on second pass...'
            return ("%s : %s" % (str(e.__class__.__name__), str(e)), 0.0)
##         print 'VALUE = ', value
        return "", value


class QueriableList(list):
    def get(self, aname):
        for o in self:
            if o.name == aname:
                return o
        return None

    def iget(self, aname):
        for i, o in enumerate(self):
            if o.name == aname:
                return i
        return None

    def delete(self, aname):
        for i, o in enumerate(self):
            if o.name == aname:
                del self[i]
                break


class BadStoichError(Exception):
    """Flags a wrong stoichiometry expression"""


class BadRateError(Exception):
    """Flags a wrong rate expression"""


class BadTypeComponent(Exception):
    """Flags an assignment of a model component to a wrong type object"""
