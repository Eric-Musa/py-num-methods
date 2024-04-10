import numpy as np
import matplotlib.pyplot as plt
from root_finding import bisection1d


def forward_euler(y0, F, t0=0, h=1e-3, eps=1e-3, verbose_step=100):
    ys = [y0]
    ts = [t0]
    k = 0
    while True:
        k += 1
        ts.append(ts[-1] + h)
        ys.append(ys[-1] + h * F(ts[-1], ys[-1]))
        if (eps < 1. and np.linalg.norm(ys[-1] - ys[-2]) < eps) or \
            (eps > 1. and k > eps):
            break
        if verbose_step and k % verbose_step == 0: print(f'Step {k}: t={ts[-1]}, Y={ys[-1]}')
    return ys, ts


def backward_euler(y0, F, t0=0, h=1e-3, eps=1e-3, verbose_step=100, root=bisection1d):
    ys = [y0]
    ys = [t0]
    k = 0
    while True:
        k += 1
        ts.append(ts[-1] + h)
        root_func = lambda x : ys[-1] - x + h*F(ts[-1]+h, x)
        y_new, _ = root(root_func, )
        ys.append(y_new)
        if (eps < 1. and np.linalg.norm(ys[-1] - ys[-2]) < eps) or \
            (eps > 1. and k > eps):
            break
        if verbose_step and k % verbose_step == 0: print(f'Step {k}: t={ts[-1]}, Y={ys[-1]}')
    return ys, ts




if __name__ == '__main__':
    k = 2
    T0 = 100
    Tenv = 25


    def newtons_cooling(t, T):
        return -k*(T - Tenv)

    
    def newtons_cooling_exact(t):
        return (T0 - Tenv) * np.exp(-k*t) + Tenv

    
    ys_fe, ts = forward_euler(T0, newtons_cooling)
