"""S-timator package"""
from pathlib import Path

from . import style
from .dynamics import scan, solve
from .estimation import s_timate
from .model import Model
from .modelparser import read_model as read_model
from .timecourse import (
    Solution,
    Solutions,
    TimeCourses,
    read_tc,
    readTCs,
)

Model.estimate = s_timate  # type: ignore
Model.solve = solve  # type: ignore
Model.scan = scan  # type: ignore

__version__ = '1.0.1'


def get_examples_path():
    import stimator.examples.timecourses as expl_tcs
    return Path(expl_tcs.__path__[0])


class VersionObj:
    def __init__(self):
        self.version = __version__
        self.fullversion = self.version
        self.date = "Dec 2026"

    def __str__(self):
        return self.version


__full_version__ = VersionObj()

__all__ = [
    "Model",
    "Solution",
    "Solutions",
    "TimeCourses",
    "__version__",
    "readTCs",
    "read_model",
    "read_tc",
    "style"
]

if __name__ == '__main__':
    print(__version__)
