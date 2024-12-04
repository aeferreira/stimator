import sympy
from sympy.abc import x

def step(t, at, top=1.0):
    if t < at:
        return 0.0
    else:
        return top


def newstep(t, at, top=1.0):
    piece_step = sympy.Piecewise((0.0, x < at), (top, x >= at))
    return piece_step.subs(x, t)

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
             "tan": sympy.tan,
             "cot": sympy.cot,
             "log": sympy.log,
             "exp": sympy.exp,
             "sign": sympy.sign,
             "abs": sympy.Abs,
             "root": sympy.root,
             "sqrt": sympy.sqrt,
             "step": step}
    return funcs

def main():
    allowed = allowed_sympy_funcs()
    for n, v in allowed.items():
        print(f'{n:>10} ----> {v} , type: {type(v)}')

if __name__ == '__main__':
    main()