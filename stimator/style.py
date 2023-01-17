from pathlib import Path
from matplotlib import pyplot as plt

from stimator.utils import _is_string

basic = {'xaxis.labellocation': 'right',
         'legend.frameon': True,
         'legend.facecolor': 'white'}


def _process_style_args(args):

    if _is_string(args) or isinstance(args, dict) or isinstance(args, Path):
        args = [args]

    # add basic mods to style if starts with st-

    new_args = []
    for arg in args:
        if _is_string(arg) and arg.startswith('st-'):
            arg = arg.replace('st-', '')
            new_args.extend([arg, basic])
        else:
            new_args.append(arg)

    # fix change of names of seaaborn styles in v3.6
    final_args = []
    for arg in new_args:
        if _is_string(arg) and arg.startswith('seaborn'):
            if arg not in plt.style.available:
                if 'v0_8' in arg:
                    arg = arg.replace('-v0_8', '')
                else:
                    arg = arg.replace('seaborn', 'seaborn-v0_8')
        final_args.append(arg)

    if len(final_args) == 1:
        return final_args[0]
    else:
        return final_args


def use(args):
    """Validates args and calls pyplot.style.use()."""
    final_args = _process_style_args(args)
    plt.style.use(final_args)


def context(args):
    final_args = _process_style_args(args)
    return plt.style.context(final_args)


if __name__ == '__main__':
    print(_process_style_args('seaborn'))
    print(_process_style_args('seaborn-whitegrid'))
    print(_process_style_args('seaborn-v0_8-whitegrid'))
    print(_process_style_args('bmh'))
    print(_process_style_args(['seaborn-whitegrid', basic]))
    print(_process_style_args({'figure.figsize': (10, 8)}))
    print(_process_style_args(['seaborn-darkgrid',
                               'st-seaborn-whitegrid',
                               {'figure.figsize': (10, 8)}]))
    use(['seaborn-darkgrid', 'st-seaborn-whitegrid',
                             {'figure.figsize': (10, 8)}])
