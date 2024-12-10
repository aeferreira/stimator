import math
import sympy
from sympy.codegen.cfunctions import log10 as symlog10
# from sympy.abc import t, a, x
# from sympy import lambdify, Piecewise

def cotangent(x):
    return 1.0 / math.tan(x)


def msign(x):
    return math.copysign(1.0, x)


def step(t, at, top=1.0):
    return top if t >= at else 0.0


step.is_rate = True


def sympystep(t, at, top=1.0):
    return sympy.Piecewise((0.0, t < at), (top, t >= at))


def sqrpulse(t, aton, atoff, top=1.0):
    if t < aton:
        return 0.0
    elif t >= aton and t <= atoff:
        return top
    else:
        return 0.0


sqrpulse.is_rate = True


def sympysqrpulse(t, aton, atoff, top=1.0):
    return sympy.Piecewise((0.0, t < aton),
                           (top, t <= atoff),
                           (0.0, t > atoff))


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


def sympystairway(t, times, values):
    args = [(0.0, t < times[0])]
    for i, time in enumerate(times[1:]):
        args.append((values[i+1], t < time))
    return sympy.Piecewise(*args)


allowed_sympy_funcs = { "sin": sympy.sin,
                        "cos": sympy.cos,
                        "tan": sympy.tan,
                        "cot": sympy.cot,
                        "log": sympy.log,
                        "log10": symlog10,
                        "exp": sympy.exp,
                        "sign": sympy.sign,
                        "abs": sympy.Abs,
                        "sqrt": sympy.sqrt,
                        "pi": sympy.pi,
                        "e": sympy.E,
                        "step": sympystep,
                        "stairway": sympystairway,
                        "sqrpulse": sympysqrpulse,
                        }


allowed_math_funcs = { "sin": math.sin,
                        "cos": math.cos,
                        "tan": math.tan,
                        "cot": cotangent,
                        "log": math.log,
                        "log10": symlog10,
                        "exp": math.exp,
                        "sign": msign,
                        "abs": abs,
                        "sqrt": math.sqrt,
                        "pi": math.pi,
                        "e": math.e,
                        "step": step,
                        "stairway": stairway,
                        "sqrpulse": sqrpulse,
                        }


def main():
    allowed = allowed_sympy_funcs
    for n, v in allowed.items():
        print(f'{n:>10} ----> {v} , type: {type(v)}')

if __name__ == '__main__':
    main()