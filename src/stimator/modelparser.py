"""This module contains code to parse a model definition text,

The class ('StimatorParser') parses text representing a valid model.
The result is a Model object.
The parsing loop relies on regular expressions."""

import re
from io import StringIO

import sympy
from sympy.parsing.sympy_parser import T, parse_expr

import stimator.model as model
from stimator.model import (dotmap_insert, dotmap_contains)

from stimator.kinetics import step

def allowed_sympy_funcs():
    funcs = {"sin": sympy.sin,
             "cos": sympy.cos,
             "tan": sympy.tan,
             "cot": sympy.cot,
             "log": sympy.log,
             "exp": sympy.exp,
             "sign": sympy.sign,
             "abs": sympy.Abs,
             "root": sympy.root,
             "sqrt": sympy.sqrt,
             "pi":sympy.pi,
             "step": step}
    return funcs


globals2useinparse = allowed_sympy_funcs()

def parse_const(model, parsed_name, valueexpr):
    """Uses builtin eval function to check for the validity
    of a math expression.

        Constants previously defined can be used"""

    print(f"\n---- parsing {parsed_name}\nwith expression {valueexpr}")
    print("constants namespace ======")
    model._all_constants.pprint()
    print("==========================")
    try:
        # value = float(eval(valueexpr,
        #                    model._usable_functions,
        #                    dict(model._all_constants)))
        value = parse_expr(
            f"float({valueexpr})",
            dict(model._all_constants),
            global_dict=globals2useinparse,
            transformations=T[4],
        )
        print("Resulting value:")
        print(value)
        print("--------------------------")
    except Exception as e:
        excpt_type = str(e.__class__.__name__)
        excpt_msg = str(e)
        if excpt_type == "SyntaxError":
            excpt_msg = "Bad expression"
        return ("%s : %s" % (excpt_type, excpt_msg), 0.0)
    return ("", value)


# ----------------------------------------------------------------------------
#         Regular expressions for grammar elements and dispatchers
# ----------------------------------------------------------------------------

identifierpattern = r"[_a-z]\w*"
reppattern = r"([_a-z]\w*|>{1,2}|~|->|\.{1,3})"
multdotidspattern = r"[_a-z]\w*(\.[_a-z]\w*)*"
fracnumberpattern = r"[-]?\d*[.]?\d+"
realnumberpattern = fracnumberpattern + r"(e[-]?\d+)?"

emptylinepattern = r"^\s*(?:#.*)?$"
constdefpattern = (
    r"^\s*(?P<name>"
    + identifierpattern
    + r")\s*=\s*(?P<value>[^#]*)(?:\s*#.*)?$"
)
varlistpattern = (
    r"^\s*variables\s*(?::\s*)?(?P<names>("
    + identifierpattern
    + r"\s*)+)(?:#.*)?$"
)
finddefpattern = (
    r"^\s*(?:find)\s+(?P<name>"
    + multdotidspattern
    + r")\s*in\s*(\[|\()\s*(?P<lower>.*)\s*,\s*(?P<upper>.*)\s*(\]|\))\s*(?:#.*)?$"
)
ratedefpattern = (
    r"^\s*(?:reaction\s+)?(?P<name>"
    + identifierpattern
    + r")\s*(:|=)\s*(?P<stoich>.*\s*(->|<=>)\s*[^,]*)\s*,(?:\s*rate\s*=)?\s*(?P<rate>[^#]+)(?:#.*)?$"
)
tcdefpattern = r"^\s*timecourse\s+?(?P<filename>[^#]+)(?:#.*)?$"
atdefpattern = (
    r"^\s*@\s*(?P<timevalue>[^#]*)\s+(?P<name>"
    + identifierpattern
    + r")\s*=\s*(?P<value>[^#]*)(?:\s*#.*)?$"
)
titlepattern = r"^\s*title\s*(?::\s*)?(?P<title>[^#]+)(?:#.*)?$"
tfpattern = r"^\s*tf\s*(?::\s*)?(?P<tf>[^#]+)(?:#.*)?$"
replistpattern = (
    r"^\s*!!\s*(?::\s*)?(?P<names>(" + reppattern + r"\s*)+)(?:#.*)?$"
)
statepattern = (
    r"^\s*(?P<name>"
    + identifierpattern
    + r")\s*=\s*(?P<value>state[^#]*)(?:\s*#.*)?$"
)
initpattern = r"^\s*(?P<name>init)\s*:\s*(?P<value>[^#]*)(?:\s*#.*)?$"
dxdtpattern = (
    r"^\s*(?P<name>"
    + identifierpattern
    + r")\s*'\s*=\s*(?P<value>[^#]*)(?:\s*#.*)?$"
)
transfpattern = (
    r"^\s*(transf|~)\s*(?P<name>"
    + identifierpattern
    + r")\s*=\s*(?P<value>[^#]*)(?:\s*#.*)?$"
)
inputvarpattern = (
    r"^\s*(input|in|->)\s*(?P<name>"
    + identifierpattern
    + r")\s*=\s*(?P<value>[^#]*)(?:\s*#.*)?$"
)

stoichpattern = r"^\s*(?P<coef>\d*)\s*(?P<variable>[_a-z]\w*)\s*$"

nameErrorpattern = r"NameError : name '(?P<name>\S+)' is not defined"
inRateErrorpattern = r".*in rate of (?P<name>\w+):"
syntaxErrorpattern = r"SyntaxError.*(?P<inrate>in rate of.*)"

identifier = re.compile(identifierpattern, re.IGNORECASE)
fracnumber = re.compile(fracnumberpattern, re.IGNORECASE)
realnumber = re.compile(realnumberpattern, re.IGNORECASE)

emptyline = re.compile(emptylinepattern)
constdef = re.compile(constdefpattern, re.IGNORECASE)
varlist = re.compile(varlistpattern, re.IGNORECASE)
finddef = re.compile(finddefpattern, re.IGNORECASE)
ratedef = re.compile(ratedefpattern, re.IGNORECASE)
statedef = re.compile(statepattern, re.IGNORECASE)
initdef = re.compile(initpattern, re.IGNORECASE)
tcdef = re.compile(tcdefpattern)
atdef = re.compile(atdefpattern)
titledef = re.compile(titlepattern)
tfdef = re.compile(tfpattern)
replistdef = re.compile(replistpattern, re.IGNORECASE)
dxdtdef = re.compile(dxdtpattern, re.IGNORECASE)
transfdef = re.compile(transfpattern, re.IGNORECASE)
invardef = re.compile(inputvarpattern, re.IGNORECASE)

stoichmatch = re.compile(stoichpattern, re.IGNORECASE)

nameErrormatch = re.compile(nameErrorpattern)
inRateErrormatch = re.compile(inRateErrorpattern, re.DOTALL)
syntaxErrormatch = re.compile(syntaxErrorpattern, re.DOTALL)

dispatchers = [
    (emptyline, "emptyLineParse"),
    (ratedef, "rateDefParse"),
    (varlist, "varListParse"),
    (finddef, "findDefParse"),
    (tcdef, "tcDefParse"),
    (atdef, "atDefParse"),
    (statedef, "stateDefParse"),
    (initdef, "initDefParse"),
    (dxdtdef, "dxdtDefParse"),
    (transfdef, "transfDefParse"),
    (invardef, "invarDefParse"),
    (constdef, "constDefParse"),
    (titledef, "titleDefParse"),
    (tfdef, "tfDefParse"),
    (replistdef, "repListDefParse"),
]

hascontpattern = r"^.*\\$"
hascontinuation = re.compile(hascontpattern)


def logicalLines(textlines):
    """generator that parses logical lines of input text.

    The backslash is a continuation character
    """
    linenum = -1
    currloc = 0
    llen = 0

    continuing = False
    for line in textlines:
        lthisline = len(line)
        llen += lthisline
        line = line.rstrip()
        tocontinue = False
        if hascontinuation.match(line):
            line = line[:-1]
            tocontinue = True
            line = line.ljust(lthisline)
        if not continuing:
            start = currloc
            linenum += 1
            logline = line
        else:
            logline += line
        if tocontinue:
            continuing = True
        else:
            end = currloc + llen
            llen = 0
            currloc = end
            yield (logline, linenum, start, end)
            continuing = False
    return


class _Physical_Line(object):
    def __init__(
        self,
        start,
        startline,
        nstartline,
        startlinepos,
        end,
        endline,
        nendline,
        endlinepos,
    ):
        self.start = start  # start pos, relative to whole text
        self.nstartline = nstartline  # start line number
        self.startline = startline  # start line
        self.startlinepos = startlinepos  # start pos, relative to start line
        self.end = end  # end relative to whole text
        self.nendline = nendline  # end line number
        self.endline = endline  # end line
        self.endlinepos = endlinepos  # end pos, relative to end line


class _Logical_Line(object):
    def __init__(self, nline, start, end, linestart, lineend):
        self.nline = nline  # logical line
        self.start = start  # start relative to logical line
        self.end = end  # end relative to logical line
        self.linestart = linestart  # start of logical line
        self.lineend = lineend  # end of logical line


def _line_from_logical_line(textlines, logpos):
    """Computes a physical line from a logical line."""

    textlines = _get_text_as_file(textlines)

    physstart = logpos.start + logpos.linestart
    physend = logpos.end + logpos.linestart

    tot = 0
    start_found = False
    for iline, line in enumerate(textlines):
        line_start_pos = tot
        tot += len(line)
        if tot > physstart and not start_found:
            nstartline = iline
            startline = line
            startlinepos = physstart - line_start_pos
            start_found = True
        if tot > physend:
            nendline = iline
            endline = line
            endlinepos = physend - line_start_pos
            return _Physical_Line(
                physstart,
                startline,
                nstartline,
                startlinepos,
                physend,
                endline,
                nendline,
                endlinepos,
            )
    return None


class StimatorParserError(Exception):
    def __init__(self, value, physloc, logloc):
        self.value = value
        self.physloc = physloc
        self.logloc = logloc

    def __str__(self):
        return str(self.value)


def _get_text_as_file(text):
    # try to open in case it is a pathname
    try:
        f = open(text)
        text = f.read()
        f.close()
    except (IOError, OSError):
        pass

    textlines = StringIO(str(text))
    return textlines


def read_model(text):
    parser = StimatorParser()
    parser.parse(text)
    if parser.error is None:
        if len(parser.tc["filenames"]) > 0:
            parser.model.metadata["timecourses"] = parser.tc["filenames"]
        if "defaultnames" in parser.tc:
            parser.model.metadata["defaultnames"] = parser.tc["defaultnames"]
        parser.model.metadata["optSettings"] = parser.optSettings
        return parser.model
    logloc = parser.errorloc
    ppos = _line_from_logical_line(text, logloc)
    raise StimatorParserError(parser.error, ppos, logloc)


# ----------------------------------------------------------------------------
#         The core StimatorParser class
# ----------------------------------------------------------------------------


class StimatorParser(object):
    def __init__(self):
        self.reset()

    def reset(self):
        self.textlines = None
        self.problemname = ""  # the name of the problem

        self.error = None  # different of None if an error occurs
        self.errorloc = None

        self.model = model.Model()
        # self.model._all_constants = dict()
        self.model._clear_constants()
        self.tc = {"filenames": []}
        # optimizer configuration
        self.optSettings = {}

        # location of timecourse def lines for error reporting
        self.tclines = []
        self.vname = []
        # location of rate def for error reporting, a list of _Logical_Line's
        self.rateloc = []

    def parse(self, text):
        "Parses a model definition text line by line"

        self.reset()

        self.textlines = _get_text_as_file(text)

        # parse the lines of text using matches
        # and dispatch to *Parse functions
        for line, nline, start, end in logicalLines(self.textlines):
            # package _Logical_Line
            loc = _Logical_Line(nline, 0, len(line), start, end)

            matchfound = False
            for d in dispatchers:
                matchresult = d[0].match(line)
                if matchresult:
                    output_function = getattr(StimatorParser, d[1])
                    output_function(self, line, loc, matchresult)
                    if self.error:
                        return  # quit on first error. Needs revision!
                    matchfound = True
                    break  # do not try any more patterns
            if not matchfound:
                self.setError("Invalid syntax:\n%s" % line, loc)
                return

        # check the validity of rate laws
        check, msg = self.model.checkRates()
        if not check:
            m = syntaxErrormatch.match(msg)
            if m:
                inrate = m.group("inrate")
                msg = "Syntax Error: bad math expression \n%s" % inrate
            # get name of transformation with offending rate
            m = inRateErrormatch.match(msg)
            if m:
                vn = m.group("name")
                indx = self.vname.index(vn)
                if vn.startswith("d_") and vn.endswith("_dt"):
                    msg = msg.replace("in rate of", "in the definition of")

                tt = self.model._get_obj_withpars(vn)
                if isinstance(tt, model.Transformation):
                    msg = msg.replace("in rate of", "in the definition of")

                self.setError(msg, self.rateloc[indx])
                rateexpr = tt()
                self.setIfNameError(msg, rateexpr, self.errorloc)
                return

    def setError(self, text, errorloc):
        self.error = text
        self.errorloc = errorloc

    def setIfNameError(self, text, exprtext, loc):
        m = nameErrormatch.match(text)
        if m:
            undefname = m.group("name")
            pos = self.errorloc.start + exprtext.find(undefname)
            loc.start = pos
            loc.end = pos + len(undefname)
            self.setError(text, loc)

    def _process_consts_in_rate(self, rate, loc):
        pardict = {}
        decls = rate.split(",")
        n_localpars = 0
        for dindex in range(len(decls) - 1, 0, -1):
            d = decls[dindex].strip()
            match = constdef.match(d)
            if match:
                name = match.group("name")
                valueexpr = match.group("value").rstrip()

                resstring, value = parse_const(self.model, name, valueexpr)
                if resstring != "":
                    loc.start = loc.start + rate.index(valueexpr)
                    loc.end = loc.start + len(valueexpr)
                    self.setError(resstring, loc)
                    self.setIfNameError(resstring, valueexpr, loc)
                    return (None, None)

                pardict[name] = value
                n_localpars += 1
            else:
                break
        rate = ",".join(decls[: len(decls) - n_localpars])
        return rate, pardict

    def rateDefParse(self, line, loc, match):
        # process name
        name = match.group("name")
        if name in self.model.reactions:  # repeated declaration
            self.setError("Repeated declaration", loc)
            return
        # process rate
        rate = match.group("rate").strip()
        stoich = match.group("stoich").strip()
        rate_loc = _Logical_Line(
            loc.nline,
            match.start("rate"),
            match.start("rate") + len(rate),
            loc.linestart,
            loc.lineend,
        )

        rate, pardict = self._process_consts_in_rate(rate, rate_loc)
        if rate is None:
            return

        if rate.endswith(".."):
            rate = rate.rstrip("..")
            resstring, value = parse_const(self.model, name, rate)
            if resstring != "":
                loc.start = match.start("rate")
                loc.end = match.start("rate") + len(rate)
                self.setError(resstring, loc)
                self.setIfNameError(resstring, rate, loc)
                return
            else:
                # it will be a float and mass action kinetics will be assumed
                rate = float(value)

        try:
            for parname, parvalue in pardict.items():
                dotmap_insert(
                    self.model._all_constants, f"{name}.{parname}", parvalue
                )
            pardict = {n: float(v) for (n, v) in pardict.items()}
            self.model.set_reaction(name, stoich, rate, pars=pardict)
        except model.BadStoichError:
            loc.start = match.start("stoich")
            loc.end = match.end("stoich")
            self.setError(
                f"'{stoich}' is an invalid stoichiometry expression", loc
            )
            return
        loc.start = match.start("rate")
        loc.end = match.end("rate")
        self.rateloc.append(loc)
        self.vname.append(name)

    def dxdtDefParse(self, line, loc, match):
        name = match.group("name")
        dxdtname = "d_%s_dt" % name
        if dxdtname in self.model.reactions:  # repeated declaration
            self.setError("Repeated declaration", loc)
            return
        expr = match.group("value").strip()
        rate_loc = _Logical_Line(
            loc.nline,
            match.start("value"),
            match.start("value") + len(expr),
            loc.linestart,
            loc.lineend,
        )
        expr, pardict = self._process_consts_in_rate(expr, rate_loc)
        if expr is None:
            return
        # setattr(self.model, name, model.variable(expr, pars=pardict))
        self.model.set_variable_dXdt(name, expr, pars=pardict)
        loc.start = match.start("value")
        loc.end = match.end("value")
        self.rateloc.append(loc)
        self.vname.append(dxdtname)

    def transfDefParse(self, line, loc, match):
        name = match.group("name")
        if name in self.model.transformations:  # repeated declaration
            self.setError("Repeated declaration", loc)
            return
        expr = match.group("value").strip()
        rate_loc = _Logical_Line(
            loc.nline,
            match.start("value"),
            match.start("value") + len(expr),
            loc.linestart,
            loc.lineend,
        )
        expr, pardict = self._process_consts_in_rate(expr, rate_loc)
        if expr is None:
            return
        for parname, parvalue in pardict.items():
            dotmap_insert(
                self.model._all_constants, f"{name}.{parname}", parvalue
            )
        pardict = {n: float(v) for (n, v) in pardict.items()}
        self.model.set_transformation(name, expr, pars=pardict)
        loc.start = match.start("value")
        loc.end = match.end("value")
        self.rateloc.append(loc)
        self.vname.append(name)

    def invarDefParse(self, line, loc, match):
        name = match.group("name")
        if name in self.model.input_variables:  # repeated declaration
            self.setError("Repeated declaration", loc)
            return
        expr = match.group("value").strip()
        rate_loc = _Logical_Line(
            loc.nline,
            match.start("value"),
            match.start("value") + len(expr),
            loc.linestart,
            loc.lineend,
        )
        expr, pardict = self._process_consts_in_rate(expr, rate_loc)
        if expr is None:
            return
        self.model.set_input_var(name, expr, pars=pardict)
        loc.start = match.start("value")
        loc.end = match.end("value")
        self.rateloc.append(loc)
        self.vname.append(name)

    def emptyLineParse(self, line, loc, match):
        pass

    def tcDefParse(self, line, loc, match):
        filename = match.group("filename").strip()
        self.tclines.append(loc.nline)
        self.tc["filenames"].append(filename)

    def stateDefParse(self, line, loc, match):
        name = match.group("name")
        state = match.group("value")
        state = state.replace("state", "self.model.set_init")

        try:
            _ = eval(state)
        except Exception:
            self.setError("Bad '%s' state definition" % name, loc)
            return

    def initDefParse(self, line, loc, match):
        name = match.group("name")
        state = match.group("value")
        if state[0] == "(" and state[-1] == ")":
            state = state[1:-1]
        state = "self.model.set_init(%s)" % state

        try:
            _ = eval(state)
        except Exception:
            self.setError("Bad '%s' state definition" % name, loc)
            return

    def constDefParse(self, line, loc, match):
        name = match.group("name")
        valueexpr = match.group("value").rstrip()

        if name in self.model.parameters:  # repeated declaration
            self.setError("Repeated declaration", loc)
            return

        resstring, value = parse_const(self.model, name, valueexpr)
        if resstring != "":
            loc.start = match.start("value")
            loc.end = match.start("value") + len(valueexpr)
            self.setError(resstring, loc)
            self.setIfNameError(resstring, valueexpr, loc)
            return

        if name in ("generations", "maxgenerations"):
            self.optSettings["generations"] = int(value)
            self.optSettings["max_generations"] = int(value)
        elif name in ("genomesize", "popsize"):
            self.optSettings["genomesize"] = int(value)
            self.optSettings["pop_size"] = int(value)
        else:
            self.model.setp(name, float(value))
            dotmap_insert(self.model._all_constants, name, value)

    def atDefParse(self, line, nline, match):
        pass  # for now

    def varListParse(self, line, loc, match):
        if "defaultnames" in self.tc:  # repeated declaration
            self.setError("Repeated declaration", loc)
            return

        names = match.group("names")
        names = names.strip()
        self.tc["defaultnames"] = names.split()

    def findDefParse(self, line, loc, match):
        name = match.group("name")

        lulist = ["lower", "upper"]
        flulist = []
        for k in lulist:
            valueexpr = match.group(k)
            resstring, v = parse_const(self.model, name, valueexpr)
            if resstring != "":
                loc.start = match.start(k)
                loc.end = match.end(k)
                self.setError(resstring, loc)
                self.setIfNameError(resstring, valueexpr, loc)
                return
            flulist.append(v)
        self.model.set_bounds(name, (flulist[0], flulist[1]))
        if name not in self.model._all_constants:
            half = sum((flulist[0], flulist[1])) / 2.0
            if not dotmap_contains(self.model._all_constants, name):
                dotmap_insert(self.model._all_constants, name, half)

    def titleDefParse(self, line, loc, match):
        title = match.group("title")
        self.model.metadata["title"] = title

    def tfDefParse(self, line, loc, match):
        title = match.group("tf")
        self.model.metadata["tf"] = title

    def repListDefParse(self, line, loc, match):
        title = match.group("names")
        self.model.metadata["!!"] = title


# ----------------------------------------------------------------------------
#         TESTING CODE
# ----------------------------------------------------------------------------

model_text = """
#This is an example of a valid model:
title: A model to test parsing.
variables: X1 X2 X3

r1 : X2  + X3 -> X1, rate = Vmax1*X2*X3 / ((k1_1+X3)*(KmX2+X2)), \\
    k1_1 = sqrt(1e-2)
leak : X3 -> 4.2 X3out, 10 ..
reaction React2 : X1 ->  2  OutVar,  \\
    step(t, 2.0, Vmax2*X1 / (Km2 + X1)) #reaction 2
kout_global = 3.14 * sign(5)
export: OutVar ->, kout * OutVar, kout = sqrt(4.0)/2.0 * kout_global,\\
    k9 = 2 * r1.k1_1

in i1 = 20 - X2
-> i2 = i1 * 15
input i3 = i1 + i2

~ totX = X2 + X1
~ OutVarmult = mult * OutVar,  mult = (kout_global/export.kout) * 2
pi   = 3.1416 * step(6, 1)
pi2  = 2*pi
pypi = pi**2  #this is pi square

another_const = r1.k1_1 ** 3

Vmax = 0.0001
#Vmax = pi*1e100**10000
#Vmax1 = 0.0001**10000
#Vmax = 100**10000
#Vmax = 100**1000**1000
find Vmax1 in [1e-9, 1e-3]
find   KmX3  in [1e-5, 1]
find KmX2 in [1e-5, pi/pi]

find Km2   in [1e-5, 1]
find Vmax2 in (1e-9, 1e-3)

find export.kout in (3,4)

@ 3.4 pi = 2*pi
x' = X3/2
#init  = state(X2 = 0.1, X3 = 0.63655, X1 = 0.0, x = 0)
init: X2 = 0.1, X3 = 0.63655, X1 = 0.0, x = 0

genomesize = 50 #should be enough
generations = 400
popsize = 20

timecourse my file.txt  # this is a timecourse filename
timecourse anotherfile.txt
#timecourse stillanotherfile.txt
tf: 10
!! X1 > X2 -> ~ ..

"""


def try2read_model(text):
    try:
        m = read_model(text)
        print("try2read_model")
        print("-------------- MODEL TEXT")
        print(text)
        print(f"\n-- Model '{m.metadata["title"]}' successfuly read:\n")
        print(m)
        print("-----------------------------------\n")
        return
    except StimatorParserError as expt:
        print("\n*****************************************")

        if expt.physloc.nstartline == expt.physloc.nendline:
            locmsg = (
                f"Error in line {expt.physloc.nendline} of model definition"
            )
        else:
            locmsg = "Error in lines %d-%d of model definition" % (
                expt.physloc.nstartline,
                expt.physloc.nendline,
            )
        print(locmsg)

        ppos = expt.physloc
        if ppos.nstartline != ppos.nendline:
            caretline = [" "] * (len(ppos.startline) + 1)
            caretline[ppos.startlinepos] = "^"
            caretline = "".join(caretline)
            value = "%s\n%s\n" % (ppos.startline.rstrip(), caretline)
            caretline = [" "] * (len(ppos.endline) + 1)
            caretline[ppos.endlinepos] = "^"
            caretline = "".join(caretline)
            value = "%s\n%s\n%s" % (value, ppos.endline.rstrip(), caretline)
        else:
            caretline = [" "] * (len(ppos.startline) + 1)
            caretline[ppos.startlinepos] = "^"
            caretline[ppos.endlinepos] = "^"
            caretline = "".join(caretline)
            value = "%s\n%s" % (ppos.startline.rstrip(), caretline)
        print(value)

        print(expt)


if __name__ == "__main__":
    try2read_model(model_text)
