#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator Parser class
# Copyright António Ferreira 2006-2009
#----------------------------------------------------------------------------
"""This module contains code to parse a model definition text,

The class ('StimatorParser') parses text representing a valid model.
The result is a Model object. 
The parsing loop relies on regular expressions."""

import os
import os.path
import StringIO
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
ratedefpattern    = r"^\s*(?:reaction\s+)?(?P<name>"+identifierpattern+r")\s*(:|=)\s*(?P<stoich>.*\s*(->|<=>)\s*.*)\s*,(?:\s*rate\s*=)?\s*(?P<rate>[^#]+)(?:#.*)?$"
tcdefpattern      = r"^\s*timecourse\s+?(?P<filename>[^#]+)(?:#.*)?$"
atdefpattern      = r"^\s*@\s*(?P<timevalue>[^#]*)\s+(?P<name>"+identifierpattern+r")\s*=\s*(?P<value>[^#]*)(?:\s*#.*)?$"
titlepattern      = r"^\s*title\s*(?::\s*)?(?P<title>[^#]+)(?:#.*)?$"
statepattern      = r"^\s*(?P<name>"+identifierpattern+r")\s*=\s*(?P<value>state[^#]*)(?:\s*#.*)?$"

stoichpattern = r"^\s*(?P<coef>\d*)\s*(?P<variable>[_a-z]\w*)\s*$"

nameErrorpattern = r"NameError : name '(?P<name>\S+)' is not defined"
inRateErrorpattern = r".*in rate of (?P<name>\S+)\s*"

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
titledef  = re.compile(titlepattern)

stoichmatch = re.compile(stoichpattern, re.IGNORECASE)

nameErrormatch   = re.compile(nameErrorpattern)
inRateErrormatch = re.compile(inRateErrorpattern, re.DOTALL)

dispatchers = [(emptyline, "emptyLineParse"),
               (ratedef,   "rateDefParse"),
               (constdef,  "constDefParse"),
               (varlist,   "varListParse"),
               (finddef,   "findDefParse"),
               (tcdef,     "tcDefParse"),
               (atdef,     "atDefParse"),
               (titledef,  "titleDefParse")]

hascontpattern  = r"^.*\\$"
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

class PhysicalLoc(object):
    def __init__(self, start, startline, nstartline, startlinepos, end, endline, nendline, endlinepos):
        self.start         = start        # start pos, relative to whole text
        self.nstartline    = nstartline   # start line number
        self.startline     = startline    # start line
        self.startlinepos  = startlinepos # start pos, relative to start line
        self.end           = end          # end relative to whole text
        self.nendline      = nendline     # end line number
        self.endline       = endline      # end line
        self.endlinepos    = endlinepos   # end pos, relative to end line

class LogicalLoc(object):
    def __init__(self, nline, start, end, linestart, lineend):
        self.nline     = nline # logical line
        self.start     = start # start relative to logical line
        self.end       = end   # end relative to logical line
        self.linestart = linestart # start of logical line
        self.lineend   = lineend   # end of logical line

def getLinesFromText(text):
    if isinstance(text,list):
        return text
    textlines = StringIO.StringIO(text)
    return textlines

def getPhysicalLineData(textlines, logpos):
    textlines = getLinesFromText(textlines)
    
    physstart = logpos.start + logpos.linestart
    physend   = logpos.end   + logpos.linestart
    
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
            return PhysicalLoc(physstart, startline, nstartline, startlinepos, physend, endline, nendline, endlinepos)
    return None

class StimatorParserError(Exception):
    def __init__(self, value, physloc, logloc):
        self.value = value
        self.physloc = physloc
        self.logloc = logloc
    def __str__(self):
        return str(self.value)
    
def read_model(text, otherdata = False):
    parser = StimatorParser()
    parser.parse(text)
    if parser.error is None:
        if otherdata:
            return (parser.model, parser.tc, parser.optSettings)
        else:
            return parser.model
    logloc = parser.errorloc
    ppos = getPhysicalLineData(text, logloc)
    raise StimatorParserError(parser.error, ppos, logloc)

def try2read_model(text):
    try:
        m, tc, os = read_model(text, True)
        print m
        print "the timecourses to load are", tc.filenames
        print
        print "the order of variables in timecourses is", tc.variablesorder
        print
        return
    except StimatorParserError, expt:
        print
        print "*****************************************"
        
        if expt.physloc.nstartline == expt.physloc.nendline:
            locmsg = "Error in line %d of model definition" % (expt.physloc.nendline)
        else:
            locmsg = "Error in lines %d-%d of model definition" % (expt.physloc.nstartline,expt.physloc.nendline)
        print locmsg
        
        ppos = expt.physloc
        if ppos.nstartline != ppos.nendline:
            caretline = [" "]*(len(ppos.startline)+1)
            caretline[ppos.startlinepos] = "^"
            caretline = ''.join(caretline)
            value = "%s\n%s\n" % (ppos.startline.rstrip(), caretline)
            caretline = [" "]*(len(ppos.endline)+1)
            caretline[ppos.endlinepos] = "^"
            caretline = ''.join(caretline)
            value = "%s\n%s\n%s" % (value, ppos.endline.rstrip(), caretline)
        else:
            caretline = [" "]*(len(ppos.startline)+1)
            caretline[ppos.startlinepos] = "^"
            caretline[ppos.endlinepos] = "^"
            caretline = ''.join(caretline)
            value = "%s\n%s" % (ppos.startline.rstrip(), caretline)
        print value

        print expt

#----------------------------------------------------------------------------
#         The core StimatorParser class
#----------------------------------------------------------------------------
class StimatorParser:
    def __init__(self):
        self.reset()
    
    def reset(self):
        self.textlines = None
        self.problemname = ""   # the name of the problem

        self.error = None      # different of None if an error occurs
        #self.errorLogLineloc = None
        self.errorloc = None
        
        # default Differential Evolution num of generations and population size
        self.model       = model.Model()
        self.tc          = timecourse.TimeCourseCollection()
        self.optSettings = {'generations':200, 'genomesize' :10}

        self.tclines     = []  #location of timecourse def lines for error reporting
        self.vname       = []
        self.rateloc     = []  #location of rate def for error reporting, a list of LogicalLoc's
    
    
    def parse (self,text):
        "Parses a model definition text line by line"

        self.reset()
        
        self.textlines = getLinesFromText(text)

        #parse the lines of text using matches and dispatch to *Parse functions
        for (line, nline,start,end) in logicalLines(self.textlines):
            #package LogicalLoc
            loc = LogicalLoc(nline, 0, len(line), start, end)
                     
            matchfound = False
            for d in dispatchers:
                matchresult = d[0].match(line)
                if matchresult:
                    output_function = getattr(StimatorParser, d[1])
                    output_function(self, line, loc, matchresult)
                    if self.error :
                        return #quit on first error. Needs revision!
                    matchfound = True
                    break #do not try any more patterns
            if not matchfound:
                self.setError("Invalid syntax", loc)
                return

        # build list of ints with the order of variables (time is at pos 0)
        if not self.tc.variablesorder:
            self.tc.intvarsorder = range(len(self.model.variables)+1)
        else:
            varnames = [x.name for x in self.model.variables]
            self.tc.intvarsorder = [varnames.index(name)+1 for name in self.tc.variablesorder]
            self.tc.intvarsorder = [0] + self.tc.intvarsorder

        # check the validity of rate laws
        check, msg = self.model.checkRates()
        if not check:
            #get name of transformation with offending rate
            m = inRateErrormatch.match(msg)
            if m:
                vn = m.group('name')
                indx = self.vname.index(vn)
                self.setError(msg, self.rateloc[indx])
                expr = getattr(self.model, vn).rate
                #~ print 'expr =',expr
                self.setIfNameError(msg, expr, self.errorloc)
                return

    def setError(self, text, errorloc):
        self.error = text
        self.errorloc = errorloc
            

    def setIfNameError(self, text, exprtext,loc):
        m = nameErrormatch.match(text)
        if m:
            undefname = m.group('name')
            pos = self.errorloc.start + exprtext.find(undefname)
            loc.start = pos
            loc.end = pos+len(undefname)
            self.setError(text, loc)

    def rateDefParse(self, line, loc, match):
        entry={}
        entryloc = {}
        #process name
        name = match.group('name')
        if model.findWithName(name, self.model.reactions): #repeated declaration
            self.setError("Repeated declaration", loc)
            return
        #process rate
        rate = match.group('rate').strip()
        stoich = match.group('stoich').strip()
        try:
            setattr(self.model, name, model.react(stoich, rate))
        except model.BadStoichError:
            loc.start = match.start(f)
            loc.end   = match.end(f)
            self.setError("'%s' is an invalid stoichiometry expression"% stoich, loc)
            return
        loc.start = match.start('rate')
        loc.end   = match.end('rate')
        self.rateloc.append(loc)
        self.vname.append(name)

    def emptyLineParse(self, line, loc, match):
        pass #do nothing

    def tcDefParse(self, line, loc, match):
        filename = match.group('filename').strip()
        self.tclines.append(loc.nline)
        self.tc.filenames.append(filename)

    def constDefParse(self, line, loc, match):
        name      = match.group('name')
        valueexpr = match.group('value').rstrip()

        if model.findWithName(name, self.model.parameters): #repeated declaration
            self.setError("Repeated declaration", loc)
            return
        
        localsdict = dict([(p.name, p) for p in self.model.parameters])

        resstring, value = test_with_consts(valueexpr, localsdict)
        if resstring != "":
            loc.start = match.start('value')
            loc.end   = match.start('value')+len(valueexpr)
            self.setError(resstring, loc)
            self.setIfNameError(resstring, valueexpr, loc)
            return
        
        if name == "generations":
            self.optSettings['generations'] = int(value)
        elif name == "genomesize":
            self.optSettings['genomesize'] = int(value)
        else:
            setattr (self.model, name, value)
        
    def atDefParse(self, line, nline, match):
        pass # for now

    def varListParse(self, line, loc, match):
        if self.tc.variablesorder: #repeated declaration
            self.setError("Repeated declaration", loc)
            return

        names = match.group('names')
        names = names.strip()
        self.tc.variablesorder = names.split()

    def findDefParse(self, line, loc, match):
        name = match.group('name')
        found = False

        localsdict = dict([(p.name, p) for p in self.model.parameters])

        lulist = ['lower', 'upper']
        flulist = []
        for k in lulist:
            valueexpr = match.group(k)
            resstring, v = test_with_consts(valueexpr, localsdict)
            if resstring != "":
                loc.start = match.start(k)
                loc.end   = match.end(k)
                self.setError(resstring, loc)
                self.setIfNameError(resstring, valueexpr, loc)
                return
            flulist.append(v)
        setattr(self.model, name, (flulist[0],flulist[1]))

    def titleDefParse(self, line, loc, match):
        title = match.group('title')
        setattr(self.model, 'title', title)

#----------------------------------------------------------------------------
#         TESTING CODE
#----------------------------------------------------------------------------

def test():
    
    modelText = """
#This is an example of a valid model:
title: Glyoxalase system in L. infantum
variables: SDLTSH TSH2 MG

Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG)*(KmTSH2+TSH2))

reaction Glx2 : SDLTSH ->  ,  \\
    Vmax2*SDLTSH / (Km2 + SDLTSH) #reaction 2

pi   = 3.1416
pi2  = 2*pi
pipi = pi**2  #this is pi square

Vmax1 = 0.0001
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
    #~ f = StringIO.StringIO(modelText)

    #~ for l,n, st, nd in logicalLines(f):
        #~ print '%d (%d,%d):'%(n,st,nd)
        #~ print "%s|"%l
        #~ print modelText[st:nd]
        #~ print len(l), len(modelText[st:nd])
    #~ print '#####################################################################'

    print '------------- test model -----------------------'
    print modelText
    print '------------------------------------------------'
    
    try2read_model(modelText)
    
    textlines = modelText.split('\n')
    print '\n======================================================'
    print 'Testing error handling...'

    textlines.insert(12,'pipipi = pois  #this is an error')
    modelText = '\n'.join(textlines)
    try2read_model(modelText)

    del(textlines[12])
    textlines.insert(6,'find pois in [1e-5, 2 + kkk]  #this is an error')
    modelText = '\n'.join(textlines)
    try2read_model(modelText)

    del(textlines[6])
    textlines.insert(6,'pipipi = pi*1e100**10000  #this is an overflow')
    modelText = '\n'.join(textlines)
    try2read_model(modelText)

    del(textlines[6])
    #test repeated declaration
    textlines.insert(12,'Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG)*(KmTSH2+TSH2))')
    modelText = '\n'.join(textlines)
    try2read_model(modelText)

    del(textlines[12])
    del(textlines[5])
    #text bad rate
    textlines.insert(5,'Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG2)*(KmTSH2+TSH2))')
    modelText = '\n'.join(textlines)
    try2read_model(modelText)

    del(textlines[5])
    #text bad rate
    textlines.insert(5,'Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG2))*(KmTSH2+TSH2))')
    modelText = '\n'.join(textlines)
    try2read_model(modelText)

    del(textlines[5])
    textlines.insert(5,'Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG)*(KmTSH2+TSH2))')
    del(textlines[8])
    #text bad rate
    textlines.insert(8,'    Vmax2*SDLTSH / (Km2 + SDLTSH)) #reaction 2')
    modelText = '\n'.join(textlines)
    try2read_model(modelText)

    del(textlines[8])
    textlines.insert(8,'    Vmax2*SDLTSH / (Km2 + SDLTSH) #reaction 2')
    textlines.insert(6,'bolas !! not good')
    modelText = '\n'.join(textlines)
    try2read_model(modelText)

    del(textlines[6])

    #~ del(textlines[27])  # delete timecourse declarations
    #~ del(textlines[27])

    #~ modelText = '\n'.join(textlines)

def profile_main():
 # This is the main function for profiling 
 # We've renamed our original main() above to real_main()
 import cProfile, pstats, StringIO
 prof = cProfile.Profile()
 prof = prof.runctx("test()", globals(), locals())
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
    test()
 
 
 