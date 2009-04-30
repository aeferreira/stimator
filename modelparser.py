#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator Parser class
# Copyright António Ferreira 2006-2009
#----------------------------------------------------------------------------
"""This module contains code to parse a model definition text,

The most important class ('StimatorParser') holds data to represent a valid model.
The parsing loop relies on regular expressions."""

import os
import os.path
import re
import math
import model
import timecourse
import utils

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
        self.errorLoc = {
                "linetext" : "", # if error, the text of the offending line
                "line"     : -1, # if error, the line number of the offending line,
                "start"    : -1, # if error, the starting position of the offending expression
                "end"      : -1} # if error, one past the ending position of the offending expression

        # default Differential Evolution num of generations and population size
        
        self.model = model.Model()
        self.tc          = timecourse.TimeCourseCollection()
        self.optSettings = {'generations':200, 'genomesize' :10}

        self.currentline = -1
        self.tclines     = []  #location of timecourse def lines for error reporting
        self.rateloc     = []  #location of rate def for error reporting, a list of {'rateline', 'ratestart', 'rateend'}
    
    
    def parse (self,modeltext):
        "Parses a model definition text line by line"

        self.reset()

        #parse the lines of text using matches and dispatch to *Parse functions
        self.currentline = 0
        for line in modeltext:
            matchfound = False
            for d in dispatchers:
                matchresult = d[0].match(line)
                if matchresult:
                    output_function = getattr(StimatorParser, d[1])
                    output_function(self, line, self.currentline, matchresult)
                    if self.error :
                        return #quit on first error. Needs revision!
                    matchfound = True
                    break #do not try any more patterns
            if not matchfound:
                self.setError("Invalid syntax", 0, len(line), self.currentline, line)
                return
            self.currentline += 1

        # build list of variables
        for v in self.model.rates:
          for rORp in ('reagents','products'):
              for c in v[rORp]:
                  var = c[0]
                  if self.model.constants.has_key(var):
                        continue # this is a constant ("external variable")
                  if not var in self.model.variables:
                      self.model.variables.append(var)

        # build list of ints with the order of variables (time is at pos 0)
        if not self.tc.variablesorder:
            self.tc.intvarsorder = range(len(self.model.variables)+1)
        else:
            self.tc.intvarsorder = [self.model.variables.index(name)+1 for name in self.tc.variablesorder]
            self.tc.intvarsorder = [0] + self.tc.intvarsorder

        # check the validity of rate laws
        varlist = self.model.variables[:]
        for i in self.model.parameters:
            varlist.append(i[0])
        for i, v in enumerate(self.model.rates):
            resstring, value = test_with_everything(v['rate'], self.model.constants, varlist)
            if resstring != "":
                errorloc = self.rateloc[i]
                self.setError(resstring, errorloc['ratestart'], errorloc['rateend'], errorloc['rateline'], modeltext[errorloc['rateline']])
                self.setIfNameError(resstring, v['rate'])
                return

    def loadTimeCourses (self,modeltext,filedir):

        if len(self.tc.filenames) == 0 :
           self.setError("No time courses to load!\nPlease indicate some time courses with 'timecourse <filename>'", -1, -1, -1, "")
           return None
        
        # check and load timecourses
        self.tc.basedir = filedir
        os.chdir(self.tc.basedir)
        pathlist = [os.path.abspath(k) for k in self.tc.filenames]

        print "-------------------------------------------------------"
        self.tc.data = []
        for (i,filename) in zip(self.tclines, pathlist):
            if not os.path.exists(filename) or not os.path.isfile(filename):
                self.setError("Time course file \n%s\ndoes not exist"% filename, 0, len(modeltext[i]), i, modeltext[i])
                return None
            h,d = timecourse.readTimeCourseFromFile(filename, atindexes=self.tc.intvarsorder)
            if d.shape == (0,0):
                self.setError("File\n%s\ndoes not contain valid time-course data"% filename, 0, len(modeltext[i]), i, modeltext[i])
                return None
            else:
                print "%d time points for %d variables read from file %s" % (d.shape[0], d.shape[1]-1, filename)
                self.tc.headers.append(h)
                self.tc.data.append(d)
        for i,d in enumerate(self.tc.data):
            if d.shape[1] != len(self.model.variables)+1:
                self.setError("There are %i initial values in time course %s but model has %i variables"%(d.shape[1]-1,
                               self.tc.filenames[i],len(self.model.variables)),-1, -1, -1, "")
                return None
        self.tc.shapes     = [i.shape for i in self.tc.data]
        self.tc.shortnames = [os.path.split(filename)[1] for filename in pathlist]
        return self.tc

    def setError(self, text, start, end, nline=None, line=None):
        self.error = text
        self.errorLoc['start'] = start
        self.errorLoc['end'] = end
        if nline != None: self.errorLoc['line'] = nline
        if line != None:  self.errorLoc['linetext'] = line

    def setIfNameError(self, text, exprtext):
        m = nameErrormatch.match(text)
        if m:
            undefname = m.group('name')
            pos = self.errorLoc['start'] + exprtext.find(undefname)
            self.setError(text, pos, pos+len(undefname))

    def rateDefParse(self, line, nline, match):
        entry={}
        entryloc = {}
        #process name
        name = match.group('name')
        found = False
        for k in self.model.rates:
            if k['name'] == name:
                found = True
                break
        if found :#repeated declaration
            self.setError("Repeated declaration", 0, len(line), nline, line)
            return

        entry['name'] = name

        #process rate
        entry['rate']     = match.group('rate').strip()
        entryloc['rateline'] = self.currentline
        entryloc['ratestart']= match.start('rate')
        entryloc['rateend']  = match.end('rate')

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
                   self.setError("'%s' is an invalid stoichiometry expression"%complexesstring, 
                                 match.start(f), match.end(f), nline, line)
                   return


        #append to list of rates
        self.model.rates.append(entry)
        self.rateloc.append(entryloc)

    def emptyLineParse(self, line, nline, match):
        pass #do nothing

    def tcDefParse(self, line, nline, match):
        filename = match.group('filename').strip()
        self.tclines.append(nline)
        self.tc.filenames.append(filename)

    def constDefParse(self, line, nline, match):
        name      = match.group('name')
        valueexpr = match.group('value').rstrip()

        if self.model.constants.has_key(name) :#repeated declaration
            self.setError("Repeated declaration", 0, len(line), nline, line)
            return

        resstring, value = test_with_consts(valueexpr, self.model.constants)
        if resstring != "":
            self.setError(resstring, match.start('value'), match.start('value')+len(valueexpr), nline, line)
            self.setIfNameError(resstring, valueexpr)
            return

        if name == "generations":
            self.optSettings['generations'] = int(value)
        elif name == "genomesize":
            self.optSettings['genomesize'] = int(value)
        else:
            self.model.constants[name]=value

    def atDefParse(self, line, nline, match):
        name      = match.group('name')
        valueexpr = match.group('value').rstrip()
        timeexpr = match.group('timevalue').rstrip()

        if not self.model.constants.has_key(name) :#constant has not been defined
            self.setError("Wrong @: constant %s has not been defined", match.start('name'), len(name), nline, line)
            return

        resstring, value = test_with_consts(valueexpr, self.model.constants)
        if resstring != "":
            self.setError(resstring, match.start('value'), match.start('value')+ len(valueexpr), nline, line)
            self.setIfNameError(resstring, valueexpr)
            return

        resstring, timevalue = test_with_consts(timeexpr, self.model.constants)
        if resstring != "":
            self.setError(resstring, match.start('timevalue'), match.start('timevalue')+ len(timeexpr), nline, line)
            self.setIfNameError(resstring, timeexpr)
            return

        self.model.atdefs.append((timevalue,name,value))

    def varListParse(self, line, nline, match):
        if self.tc.variablesorder: #repeated declaration
            self.setError("Repeated declaration", 0, len(line), nline, line)
            return

        names = match.group('names')
        names = names.strip()
        self.tc.variablesorder = names.split()

    def findDefParse(self, line, nline, match):
        name = match.group('name')
        found = False
        for k in self.model.parameters:
            if k[0] == name:
                found = True
                break
        if found :#repeated declaration
            self.setError("Repeated declaration", 0, len(line), nline, line)
            return

        lulist = ['lower', 'upper']
        flulist = []
        for k in lulist:
            valueexpr = match.group(k)
            resstring, v = test_with_consts(valueexpr, self.model.constants)
            if resstring != "":
                self.setError(resstring, match.start(k), match.end(k), nline, line)
                self.setIfNameError(resstring, valueexpr)
                return
            flulist.append(v)
        entry = tuple([name,flulist[0],flulist[1]])
        self.model.parameters.append(entry)

#----------------------------------------------------------------------------
#         TESTING CODE
#----------------------------------------------------------------------------

def printParserResults(parser):
    
    if parser.error:
        #candy syntax: upgrade dict to dot lookup style dict
        errorLoc = utils.DictDotLookup(parser.errorLoc)
        print
        print "*****************************************"
        if errorLoc.line != -1:
            print "ERROR in line %d:" % (errorLoc.line)
            print errorLoc.linetext
            caretline = [" "]*(len(errorLoc.linetext)+1)
            caretline[errorLoc.start] = "^"
            caretline[errorLoc.end] = "^"
            caretline = "".join(caretline)
            print caretline
        print parser.error
        return
    parser.model.pprint()
    print "the timecourses to load are", parser.tc.filenames
    print
    print "the order of variables in timecourses is", parser.tc.variablesorder
    print
    

def real_main():
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
    print '------------- test model -----------------------'
    print modelText
    print '------------------------------------------------'
    parser = StimatorParser()
    textlines = modelText.split("\n")

    parser.parse(textlines)
    printParserResults(parser)
    
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
    textlines.insert(5,'Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG)*(KmTSH2+TSH2))')
    textlines.insert(6,'bolas !! not good')

    parser.parse(textlines)
    printParserResults(parser)
    del(textlines[6])

    del(textlines[25])  # delete timecourse declarations
    del(textlines[25])

    parser.parse(textlines)
    printParserResults(parser)

def profile_main():
 # This is the main function for profiling 
 # We've renamed our original main() above to real_main()
 import cProfile, pstats, StringIO
 prof = cProfile.Profile()
 prof = prof.runctx("real_main()", globals(), locals())
 stream = StringIO.StringIO()
 stats = pstats.Stats(prof, stream=stream)
 stats.sort_stats("time")  # Or cumulative
 stats.print_stats(80)  # 80 = how many to print
 # The rest is optional.
 # stats.print_callees()
 # stats.print_callers()
 print stream.getvalue()
 #logging.info("Profile data:\n%s", stream.getvalue())


if __name__ == "__main__":
    real_main()
 
 
 