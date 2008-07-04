#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator Parser class
# Antonio Ferreira may 2006
#----------------------------------------------------------------------------
"""This module contains code to parse a model definition text,

The most important class ('StimatorParser') holds data to represent a valid model.
The parsing loop relies on regular expressions."""

import re
import math
from numpy import *

#----------------------------------------------------------------------------
#         Functions to check the validity of math expressions
#----------------------------------------------------------------------------
def test_with_consts(valueexpr, consts={}):
    """Uses builtin eval function to check for the validity of a math expression.

       Constants previously defined can be used"""

    locs = consts.copy()
    try :
       value = float(eval(valueexpr, vars(math), locs))
    except Exception, e:
       return ("%s : %s"%(str(e.__class__.__name__),str(e)), 0.0)
    return ("", value)


def test_with_everything(valueexpr, consts={}, varlist=[]):
    locs = consts.copy()
    #part 1: nonpermissive, except for NameError
    try :
       value = float(eval(valueexpr, vars(math), locs))
    except NameError:
       pass
    except Exception, e:
       return ("%s : %s"%(str(e.__class__.__name__),str(e)), 0.0)
    #part 2: permissive, with dummy values (1.0) for vars and parameters
    vardict = {}
    for i in varlist:
        vardict[i]=1.0
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
identifierpattern = r"[_a-z]\w*"
fracnumberpattern = r"[-]?\d*[.]?\d+"
realnumberpattern = fracnumberpattern + r"(e[-]?\d+)?"

emptylinepattern  = r"^\s*(?:#.*)?$"
constdefpattern   = r"^\s*(?P<name>"+identifierpattern+r")\s*=\s*(?P<value>[^#]*)(?:\s*#.*)?$"
varlistpattern    = r"^\s*variables\s*(?::\s*)?(?P<names>("+identifierpattern+r"\s*)+)(?:#.*)?$"
finddefpattern    = r"^\s*(?:find)\s+(?P<name>"+identifierpattern+r")\s*in\s*\[\s*(?P<lower>.*)\s*,\s*(?P<upper>.*)\s*\]\s*(?:#.*)?$"
ratedefpattern    = r"^\s*(?:reaction\s+)?(?P<name>"+identifierpattern+r")\s*:\s*(?P<reagents>.*)\s*(?P<irreversible>->|<=>)\s*(?P<products>.*)\s*,\s*rate\s*=\s*(?P<rate>[^#]+)(?:#.*)?$"
tcdefpattern      = r"^\s*timecourse\s+?(?P<filename>[^#]+)(?:#.*)?$"
atdefpattern      = r"^\s*@\s*(?P<timevalue>[^#]*)\s+(?P<name>"+identifierpattern+r")\s*=\s*(?P<value>[^#]*)(?:\s*#.*)?$"

stoichpattern = r"^\s*(?P<coef>\d*)\s*(?P<variable>[_a-z]\w*)\s*$"

nameErrorpattern = r"NameError : name '(?P<name>\S+)' is not defined"

identifier = re.compile(identifierpattern, re.IGNORECASE)
fracnumber = re.compile(fracnumberpattern, re.IGNORECASE)
realnumber = re.compile(realnumberpattern, re.IGNORECASE)

emptyline = re.compile(emptylinepattern)
constdef  = re.compile(constdefpattern,    re.IGNORECASE)
varlist   = re.compile(varlistpattern,     re.IGNORECASE)
finddef   = re.compile(finddefpattern,     re.IGNORECASE)
ratedef   = re.compile(ratedefpattern,     re.IGNORECASE)
tcdef     = re.compile(tcdefpattern)
atdef     = re.compile(atdefpattern)

stoichmatch = re.compile(stoichpattern, re.IGNORECASE)

nameErrormatch = re.compile(nameErrorpattern)

dispatchers = [(emptyline, "emptyLineParse"),
               (constdef,  "constDefParse"),
               (varlist,   "varListParse"),
               (finddef,   "findDefParse"),
               (ratedef,   "rateDefParse"),
               (tcdef,     "tcDefParse"),
               (atdef,     "atDefParse")]

#----------------------------------------------------------------------------
#         The core StimatorParser class
#----------------------------------------------------------------------------
class StimatorParser:
    def __init__(self):
        self.reset()
    
    def reset(self):
        self.problemname = ""   # the name of the problem
        self.error = None      # different of None if an error occurs
        self.errorlinetext = "" # if error, the text of the offending line
        self.errorline   = -1   # line number of currently parsed line. If error, the line number of the offending line,
        self.errorstart  = -1   # if error, the starting position of the offending expression
        self.errorend    = -1   # if error, one past the ending position of the offending expression

        self.generations = 200  # default differential evolution num of generations
        self.genomesize  = 10   # default differential evolution population size

        self.rates       = []  # a list of {'name', 'reagents', 'products', 'rate', 'irreversible' and rate locations}
        self.variables   = []  # a list of names of variables (order matters)
        self.constants   = {}  # a 'name':value dictionary
        self.parameters  = []  # a list of (name,min,max)
        self.timecourses = []  # a list of filenames
        self.atdefs      = []  # a list of (time,name,newvalue)

        self.stoichmatrixrows = []  #sparse, using reactionname:coef dictionaries

    
    
    def parse (self,modeltext):
          "Parses a model definition text line by line"
          
          self.reset()

          #parse the lines of text using matches and dispatch to *Parse functions
          self.errorline = 0
          for line in modeltext:
              matchfound = False
              for d in dispatchers:
                  matchresult = d[0].match(line)
                  if matchresult:
                     #print line
                     #print "Dispatcher:", d[1]
                     output_function = getattr(StimatorParser, d[1])
                     output_function(self, line, matchresult)
                     if self.error :
                         self.errorlinetext = line
                         return #quit on first error. Needs revision!
                     matchfound = True
                     break #do not try any more patterns
              if not matchfound:
                  self.errorlinetext = line
                  self.setError("Invalid syntax", 0, len(line))
                  return
              self.errorline = self.errorline + 1

          # build stoichiometry matrix row-wise
          nvars = len(self.variables)
          for i in range(nvars):
              self.stoichmatrixrows.append({})
          for v in self.rates:
              fields = [('reagents',-1.0),('products',1.0)]
              for f, rpf in fields:
                  for c in v[f]:
                      coef = c[1]*rpf
                      var  = c[0]
                      if self.constants.has_key(var):
                            continue # there are no rows for constants in stoich. matrix
                      if not(var in self.variables):
                            self.setError("variable %s has not been declared anywhere"%var, -1, len(var))
                            return
                      ivariable = self.variables.index(var) # error handling here
                      self.stoichmatrixrows[ivariable][v['name']] = coef

          # check the validity of rate laws
          varlist=self.variables[:]
          for i in self.parameters:
              varlist.append(i[0])
          for v in self.rates:
              resstring, value = test_with_everything(v['rate'], self.constants, varlist)
              if resstring != "":
                  self.errorline = v['rateline']
                  self.errorlinetext = modeltext[self.errorline]
                  self.setError(resstring, v['ratestart'], v['rateend'])
                  self.setIfNameError(resstring, v['rate'])
                  return

    def setError(self, text, start, end):
        self.error = text
        self.errorstart = start
        self.errorend = end

    def setIfNameError(self, text, exprtext):
        m = nameErrormatch.match(text)
        if m:
            undefname = m.group('name')
            pos = self.errorstart + exprtext.find(undefname)
            self.setError(text, pos, pos+len(undefname))


    def rateDefParse(self, line, match):
        entry={}
        #process name
        name = match.group('name')
        found = False
        for k in self.rates:
             if k['name'] == name:
                found = True
                break
        if found :#repeated declaration
             self.setError("Repeated declaration", 0, len(line))
             return

        entry['name'] = name

        #process rate
        entry['rate']     = match.group('rate').strip()
        entry['rateline'] = self.errorline
        entry['ratestart']= match.start('rate')
        entry['rateend']  = match.end('rate')

        #process irreversible
        irrsign = match.group('irreversible')
        entry['irreversible']= irrsign == "->"

        #process stoichiometry
        fields = ['reagents','products']
        for f in fields:
            entry[f]=[]
            complexesstring = match.group(f).strip()
            if len(complexesstring)==0:  #empty complexes allowed
                continue
            complexcomps = complexesstring.split("+")
            for c in complexcomps:
                m = stoichmatch.match(c)
                if m:
                   coef = m.group('coef')
                   var = m.group('variable')
                   if coef == "":
                      coef = 1.0
                   else:
                      coef = float(coef)
                   if coef == 0.0: continue # a coef equal to zero means ignore

                   entry[f].append((var,coef))
                else:
                   self.setError("'%s' is an invalid stoichiometry expression"%complexesstring, match.start(f), match.end(f) )
                   return


        #append to list of rates
        self.rates.append(entry)

    def emptyLineParse(self, line, match):
        pass #do nothing

    def tcDefParse(self, line, match):
        filename = match.group('filename').strip()
        self.timecourses.append(filename)

    def constDefParse(self, line, match):
        name      = match.group('name')
        valueexpr = match.group('value').rstrip()

        if self.constants.has_key(name) :#repeated declaration
             self.setError("Repeated declaration", 0, len(line))
             return

        resstring, value = test_with_consts(valueexpr, self.constants)
        if resstring != "":
           self.setError(resstring, match.start('value'), match.start('value')+len(valueexpr) )
           self.setIfNameError(resstring, valueexpr)
           return

        if name == "generations":
            self.generations = int(value)
        elif name == "genomesize":
            self.genomesize = int(value)
        else:
            self.constants[name]=value

    def atDefParse(self, line, match):
        name      = match.group('name')
        valueexpr = match.group('value').rstrip()
        timeexpr = match.group('timevalue').rstrip()

        if not self.constants.has_key(name) :#constant has not been defined
             self.setError("Wrong @: constant %s has not been defined", match.start('name'), len(name))
             return

        resstring, value = test_with_consts(valueexpr, self.constants)
        if resstring != "":
           self.setError(resstring, match.start('value'), match.start('value')+ len(valueexpr) )
           self.setIfNameError(resstring, valueexpr)
           return

        resstring, timevalue = test_with_consts(timeexpr, self.constants)
        if resstring != "":
           self.setError(resstring, match.start('timevalue'), match.start('timevalue')+ len(timeexpr) )
           self.setIfNameError(resstring, timeexpr)
           return

        self.atdefs.append((timevalue,name,value))

    def varListParse(self, line, match):
        if len(self.variables) > 0  :#repeated declaration
             self.setError("Repeated declaration", 0, len(line))
             return

        names = match.group('names')
        names = names.strip()
        self.variables = names.split()

    def findDefParse(self, line, match):
        name = match.group('name')
        found = False
        for k in self.parameters:
             if k[0] == name:
                 found = True
                 break
        if found :#repeated declaration
             self.setError("Repeated declaration", 0, len(line))
             return

        lulist = ['lower', 'upper']
        flulist = []
        for k in lulist:
            valueexpr = match.group(k)
            resstring, v = test_with_consts(valueexpr, self.constants)
            if resstring != "":
               self.setError(resstring, match.start(k), match.end(k))
               self.setIfNameError(resstring, valueexpr)
               return
            flulist.append(v)
        entry = tuple([name,flulist[0],flulist[1]])
        self.parameters.append(entry)

    def rateCalcString(self, rateString):
        if self.error:
           return ""
        nvars = len(self.variables)
        # replace varnames
        for i in range(nvars):
          rateString = re.sub(r"\b"+ self.variables[i]+r"\b", "variables[%d]"%i, rateString)
        # replace parameters
        for i in range(len(self.parameters)):
          rateString = re.sub(r"\b"+ self.parameters[i][0]+r"\b", "m_Parameters[%d]"%i, rateString)
        # replace constants
        for const in self.constants.keys():
          rateString = re.sub(r"\b"+ const + r"\b", "%e"% self.constants[const], rateString)
        return rateString
        

    def ODEcalcString(self, scale=1.0):
        if self.error:
           return ""
        nvars = len(self.variables)
        result = "def calcDerivs(variables,t):\n\tglobal m_Parameters\n"
        result +="\tderivatives = empty(%d)\n" % nvars

        #write @ definitions    #TODO!!!
        #~ for k in parser.atdefs:
              #~ vline = "\tif(solution_time*scale >= %g) %s = %g;\n" % (k[0], k[1], k[2])
              #~ result = result + vline

        #write rates
        for k in self.rates:
              vline = k['rate']
              # replace varnames
              for i in range(nvars):
                  vline = re.sub(r"\b"+ self.variables[i]+r"\b", "variables[%d]"%i, vline)
              # replace parameters
              for i in range(len(self.parameters)):
                  vline = re.sub(r"\b"+ self.parameters[i][0]+r"\b", "m_Parameters[%d]"%i, vline)
              # replace constants
              for const in self.constants.keys():
                  vline = re.sub(r"\b"+ const + r"\b", "%e"% self.constants[const], vline)
              vline = "\tv_%s = %s\n" %(k['name'],vline)
              result = result + vline

        #write differential equations
        for k in range(nvars):
              row = self.stoichmatrixrows[k]
              eqline = "\tderivatives[%d] = " % k
              for r in row.keys():
                  if row[r]>0:
                      eqline = eqline + "+"
                  eqline = eqline + ("%g * v_%s " %(float(row[r]), r))
              result = result + eqline + "\n"
              

        return result + '\treturn derivatives*%f' % scale

#----------------------------------------------------------------------------
#         TIME COURSE READING FUNCTION
#----------------------------------------------------------------------------

def readTimeCourseFromFile(filename):
    """Reads a time course from file.
    
    Returns a header with variable names (possibly empty) and a 2D numpy array with data.
    """
    
    header = []
    nvars = 0
    rows = []
    headerFound = False
    t0found = False
    f = open (filename)
    for line in f:
        line = line.strip()
        if len(line) == 0:continue          #empty lines are skipped
        if line.startswith('#'): continue   #comment lines are skipped
        
        items = line.split()
        
        if identifier.match(items[0]):
            if not headerFound and not t0found:
                header = filter (identifier.match, items)
                headerFound = True
            else:
                continue
        elif not realnumber.match(items[0]):
            continue
        else:
            if not t0found:
                nvars = len(items)
                t0found = True
            temprow = [nan]*nvars
            for i in enumerate(items):
                if realnumber.match(i[1]):
                    temprow[i[0]] = float(i[1])
            rows.append(temprow)
    f.close()
    
    return header, array(rows)


#----------------------------------------------------------------------------
#         TESTING CODE
#----------------------------------------------------------------------------

def printParserResults(parser):
    if parser.error:
         print
         print "*****************************************"
         print "ERROR in line %d:" % (parser.errorline)
         print parser.errorlinetext
         caretline = [" "]*(len(parser.errorlinetext)+1)
         caretline[parser.errorstart] = "^"
         caretline[parser.errorend] = "^"
         caretline = "".join(caretline)
         print caretline
         print parser.error
         return

    print
    print "the variables are" , parser.variables
    print
    print "the constants are"
    for k in parser.constants.keys():
           print "%s = %g" % (k, parser.constants[k])
    print
    print "the parameters to find are"
    for k in parser.parameters:
          print k[0],"from", k[1], "to", k[2]
    print
    print "the timecourses to load are", parser.timecourses
    print
    print "the reactions are"
    for k in parser.rates:
          irrstring = ""
          if k['irreversible']: irrstring = "(irreversible)"
          print k['name'], irrstring, ":"
          print " reagents:", k['reagents']
          print " products:", k['products']
          print " rate =", k['rate']
    print
    print "the @ definitions are"
    for k in parser.atdefs:
          print "@", k[0], k[1], "=", k[2]
    print
    print "the rows of the stoichiometry matrix are"
    for k in range(len(parser.variables)):
          name = parser.variables[k]
          row = parser.stoichmatrixrows[k]
          print "for", name, ":"
          for r in row.keys():
              print "%g %s" % (row[r], r)
    print

if __name__ == "__main__":
    modelText = """
#This is an example of a valid model:

variables: SDLTSH TSH2 MG

Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG)*(KmTSH2+TSH2))

reaction Glx2 : SDLTSH ->  , rate = Vmax2*SDLTSH / (Km2 + SDLTSH) #reaction 2

pi   = 3.1416
pi2  = 2*pi
pipi = pi**2  #this is pi square

find Vmax1 in [1e-9, 1e-3]
find   KmMG  in [1e-5, 1]
find KmTSH2 in [1e-5, pi/pi]

find Km2   in [1e-5, 1]
find Vmax2 in [1e-9, 1e-3]

@ 3.4 pi = 2*pi

genomesize = 50 #should be enough
generations = 400

timecourse my file.txt  # this is a timecourse filename
timecourse anotherfile.txt
#timecourse stillanotherfile.txt

"""

    parser = StimatorParser()
    textlines = modelText.split("\n")

    parser.parse(textlines)
    printParserResults(parser)
    sss = parser.ODEcalcString()
    print 'ODEcalcString:'
    print '------------------------------------------------'
    print sss
    print '------------------------------------------------'
    cc = compile(sss, 'bof.log','exec')
    print 'Result from compile:', cc
    exec cc
    print 'Result from exec:', calcDerivs
    
    m_Parameters = (1,1,1,1,1)
    Xpoint = (0,0,0)
    print 'calcDerivs(',Xpoint, ',0) with parameters =', m_Parameters, 'equals',
    print calcDerivs((0,0,0),0)
    
    print '\n======================================================'
    print 'The following modifications produce errors............'

    textlines.insert(12,'pipipi = pois  #this is an error')

    parser.parse(textlines)
    printParserResults(parser)

    del(textlines[12])
    textlines.insert(6,'find pois in [1e-5, 2 + kkk]  #this is an error')

    parser.parse(textlines)
    printParserResults(parser)

    del(textlines[6])
    textlines.insert(6,'pipipi = pi*1e100**10000  #this is an overflow')

    parser.parse(textlines)
    printParserResults(parser)

    del(textlines[6])
    textlines.insert(12,'Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG)*(KmTSH2+TSH2))')

    parser.parse(textlines)
    printParserResults(parser)

    del(textlines[12])
    del(textlines[5])
    textlines.insert(5,'Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG2)*(KmTSH2+TSH2))')

    parser.parse(textlines)
    printParserResults(parser)

    del(textlines[5])
    textlines.insert(5,'bolas !! not good')

    parser.parse(textlines)
    printParserResults(parser)


    print '\n======================================================'
    print 'Testing reading a time course...........................'

    demodata = """
#this is demo data with a header
t x y z
0       1 0         0
0.1                  0.1

  0.2 skip 0.2 skip this
nothing really usefull here
- 0.3 0.3 this line should be skipped
#0.4 0.4
0.5  - 0.5 - -
0.6 0.6 0.8 0.9

"""

    w = open('bof.txt', 'w')
    w.write(demodata)
    w.close()
    
    
    h, d =  readTimeCourseFromFile('bof.txt')   
    print 'source:\n------------------------------------'
    print demodata
    print '------------------------------------'
    print '\nheader:'
    print h
    print '\ndata'
    print d
    
    h, d =  readTimeCourseFromFile('examples\\TSH2b.txt')   
    print '\n\n================================================'
    print '\n\nData from TSH2b.txt'
    print 'header:'
    print h
    print '\ndata'
    print d
    print 'dimensions are %d by %d'% d.shape

    print
    raw_input("Press ENTER to finish...")
