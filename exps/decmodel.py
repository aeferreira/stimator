""" second experiment with function introspection and function decoration """
import sys
import os.path

#append parent directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir))

import inspect
import tokenize
import StringIO
from stimator import *

def model(func):
    m = Model()
    ll = inspect.getsourcelines(func)[0]
    llall = ''.join(ll)
    # tokenize the string
    tt = [ (tokenize.tok_name[toknum], tokval,a,b) for toknum, tokval, a, b, _  in tokenize.generate_tokens(StringIO.StringIO(llall).readline)]
            
    insideIndent = False
    currtks = []
    for t in tt:
        if t[0] == 'INDENT':
            insideIndent = True
            continue
        if t[0] == 'DEDENT':
            insideIndent = False
            continue
        if not insideIndent:
            continue
        if t[0] in ('COMMENT','NL'):
            continue
        if t[0] != 'NEWLINE':
            currtks.append(t)
            continue
        #~ # print the logical line
        #~ startl, startc = currtks[0][2]
        #~ endl, endc = currtks[-1][3]
        #~ if startl == endl:
            #~ text = ll[startl-1][startc:endc]
        #~ else:
            #~ text = ll[startl-1][startc:]
            #~ i = startl+1
            #~ while i < endl:
                #~ text += ll[i-1]
                #~ i +=1
            #~ text += ll[i-1][:endc]
        #~ print text
        
        #~ print currtks
        
        #process tokens
        if len(currtks) < 2:
            print "ILLEGAL LINE"
            currtks = []
            continue
        tk1,op2 = currtks[0][0], currtks[1][1]
        if tk1 != 'NAME' or op2 != '=':
            print "ILLEGAL LINE"
            currtks = []
            continue
        
        name = currtks[0][1]
        res = [(tk[0], tk[1]) for tk in currtks[2:]]
        if len(res) == 1:
            rhs = tokenize.untokenize(res)
        else:
            r1,r2 = res[0], res[1]
            if r1[0] == 'STRING' and r2[1] == ',':
                rate = tokenize.untokenize(res[2:])
                if len(res) > 3:
                    rate = '"%s"' % rate
                rhs = '%s, %s' % (r1[1],rate)
            else:
                rhs = tokenize.untokenize(res)
        #~ print name
        #~ print rhs
        #~ print "---------------------"
        setattr(m, name, eval(rhs))
            
        currtks = []    
    return m
