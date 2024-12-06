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


# def newstep(t, at, top=1.0):
#     piece_step = sympy.Piecewise((0.0, x < at), (top, x >= at))
#     return piece_step.subs(x, t)

# step = lambdify([t, a], Piecewise((0, t < a), (1, t >= a)))
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
                        "pi": math.pi,
                        "e": math.e,
                        "step": step,
                        "stairway": stairway,
                        "sqrpulse": sqrpulse,
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