import sympy

def step(t, at, top=1.0):
    if t < at:
        return 0.0
    else:
        return top


step.is_rate = True


def sqrpulse(t, aton, atoff, top=1.0):
    if t < aton:
        return 0.0
    elif t >= aton and t <= atoff:
        return top
    else:
        return 0.0


sqrpulse.is_rate = True


def stairway(t, times, values):
    if len(times) == 0:
        return 0.0
    if t < times[0]:
        return 0.0
    value = 0.0
    for i, time in enumerate(times):
        if t > time:
            value = values[i]
    return value


stairway.is_rate = True

def allowed_sympy_funcs():
    funcs = {"sin": sympy.sin,
             "cos": sympy.cos,
             "log": sympy.log,
             "exp": sympy.exp}
    return funcs

def main():
    allowed = allowed_sympy_funcs()
    for n, v in allowed.items():
        print(f'{n:>10} ----> {v} , type: {type(v)}')

if __name__ == '__main__':
    main()