# -*- coding: utf-8 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator kinetics functions
# Copyright António Ferreira 2006-2015
#----------------------------------------------------------------------------

def step (t, at, top=1.0):
    if t < at:
        return 0.0
    else:
        return top

#mark step as a rate law
step.is_rate = True

def sqrpulse (t, aton, atoff, top=1.0):
    if t < aton:
        return 0.0
    elif t >= aton and t <= atoff:
        return top
    else:
        return 0.0

sqrpulse.is_rate = True