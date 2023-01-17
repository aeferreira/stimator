import math
from itertools import cycle
from string import Template
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from cycler import cycler, Cycler

from packaging.version import Version

from stimator.utils import _is_string, _is_sequence

def use(*args):
    """Validates args and calls pyplot.style.use()."""

    new_args = []
    for arg in args:
        if _is_string(arg) and arg.starts_with('seaborn'):
            if not arg in plt.style.available:
                if 'v0_8' in arg:
                    arg = arg.replace('-v0_8', '')
                else:
                    arg = arg.replace('seaborn', 'seaborn-v0_8')
        new_args.append(arg)


