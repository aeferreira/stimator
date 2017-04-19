import os
import sys
import pytest

_curr_dir = os.getcwd()
_THIS_DIR, _ = os.path.split(os.path.abspath(__file__))

def run_tests():
    os.chdir(_THIS_DIR)
    pytest.main()
    os.chdir(_curr_dir)
