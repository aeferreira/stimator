from pathlib import Path

from matplotlib import pyplot as plt

from stimator.utils import _is_string


def _get_stimator_settings():
    return {'xaxis.labellocation': 'right',
            'legend.frameon': True,
            'legend.facecolor': 'white'}

def _process_style_args(args):

    basic = _get_stimator_settings()

    if _is_string(args) or isinstance(args, (dict, Path)):
        args = [args]

    # add basic mods to style if starts with st-

    new_args = []
    for arg in args:
        if _is_string(arg) and arg.startswith('st-'):
            arg = arg.replace('st-', '')
            new_args.extend([arg, basic])
        else:
            new_args.append(arg)

    # fix change of names of seaborn styles in v3.6
    final_args = []
    for arg in new_args:
        is_seaborn = _is_string(arg) and arg.startswith('seaborn')
        if is_seaborn and (arg not in plt.style.available):
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
    basic = _get_stimator_settings()
    print(_process_style_args('seaborn'))
    print(_process_style_args('seaborn-whitegrid'))
    print(_process_style_args(['st-seaborn-whitegrid']))
    original = {s:v for s,v in plt.rcParams.items()}
    for s, v in plt.rcParams.items():
        print(f"{s} : {v}")
    use('st-seaborn-whitegrid')
    print("----------- after use st-seaborn-whitegrid")
    for s, v in plt.rcParams.items():
        if v != original[s]:
            print(f"{s} : {v}")
    use('seaborn-talk')
    print("----------- after use seaborn-talk")
    for s, v in plt.rcParams.items():
        if v != original[s]:
            print(f"{s} : {v}")
    print(_process_style_args('seaborn-v0_8-whitegrid'))
    print(_process_style_args('bmh'))
    print(_process_style_args({'figure.figsize': (10, 8)}))
    print(_process_style_args(['seaborn-darkgrid',
                               'st-seaborn-whitegrid',
                               {'figure.figsize': (10, 8)}]))
    use(['seaborn-darkgrid', 'st-seaborn-whitegrid',
                             {'figure.figsize': (3, 2)}])
    print("----------- after use of a mix")
    for s, v in plt.rcParams.items():
        if v != original[s]:
            print(f"{s} : {v}")
