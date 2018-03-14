import os
import pytest

_THIS_DIR, _ = os.path.split(os.path.abspath(__file__))

def run_tests():
    pytest.main([_THIS_DIR])
