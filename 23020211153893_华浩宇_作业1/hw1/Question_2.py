import random
import numpy as np


def linesearch_secant(f, df, d, x, alpham, rho, t):
    '''
    线性搜索子函数，函数f，导数df，当前迭代点x和当前搜索方向d
    '''
    flag = 0

    a = 0
    b = alpham
    fk = f(x)
    gk = df(x)

    phi0 = fk
    dphi0 = np.dot(gk, d)
    alpha = b * random.uniform(0, 1)

    while (flag == 0):
        newfk = f(x + alpha * d)
        phi = newfk
        if (phi - phi0) <= (rho * alpha * dphi0):
            if (phi - phi0) >= ((1 - rho) * alpha * dphi0):
                flag = 1
            else:
                a = alpha
                b = b
                if (b < alpham):
                    alpha = (a + b) / 2
                else:
                    alpha = t * alpha
        else:
            a = a
            b = alpha
            alpha = (a + b) / 2
    return alpha
