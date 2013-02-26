"""S-timator package"""

from model import Model, react, transf, state, variable, register_kin_func, get_name, set_name
from analysis import solve, plot, ModelSolver
from timecourse import readTCs, Solutions, TimeCourses
from modelparser import read_model
from version_info import __version__
import examples
