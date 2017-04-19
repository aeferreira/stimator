import os
import sys

_THIS_DIR, _ = os.path.split(os.path.abspath(__file__))
#print(_THIS_DIR)
_UPPER, _ = os.path.split(_THIS_DIR)
#print(_UPPER)
sys.path.insert(0, _UPPER)

import tests

def run_tests():
    tests.run_tests()
