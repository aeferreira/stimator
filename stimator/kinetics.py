#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# S-timator kinetics functions
# Copyright António Ferreira 2006-2011
#----------------------------------------------------------------------------

def step (t, at, top=1.0):
    if t < at:
        return 0.0
    else:
        return top

#mark step as a rate law
step.is_rate = True