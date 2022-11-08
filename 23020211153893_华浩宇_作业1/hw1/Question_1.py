import numpy as np


def f(x):
    return np.array((2 * x - 1) ** 3 + 4 * (4 - 1024 * x) ** 3)


def secant(f, x0, e=0.00001):
    assert (len(x0) == 2)
    x1, x0 = x0[1], x0[0]
    x = []
    iter = 1
    while True:
        x2 = x1 - (x1 - x0) * f(x1) / (f(x1) - f(x0))
        x.append(x2)
        iter = iter + 1
        tol = abs(x2 - x1)
        if tol < e * abs(x1):
            return x2, x
        else:
            x0 = x1
            x1 = x2


if __name__ == '__main__':
    x0 = [0, 1]
    x, _ = secant(f, x0)
    print("方程的根为{}".format(x))
