import numpy as np
from fractions import Fraction as f


def getinput():
    global m, n
    m = int(input('输入约束的个数：'))
    n = int(input('输入变量的个数：'))
    a = []
    for i in range(m):
        a.append(
            [eval(input('a({}{})='.format(i + 1, j + 1))) for j in range(n)]
        )
    e = [list(i) for i in np.diag([1] * m)]
    b = [eval(input('b({})='.format(i + 1))) for i in range(m)]
    r_1 = []
    for i in range(n):
        r_1.append(sum([a[j][i] for j in range(m)]))
    r_1 += [0] * m
    r_1.append(sum(b))
    r = [eval(input('r({})='.format(i + 1))) for i in range(n)] + [0] * (m + 1)
    vect = [n + i + 1 for i in range(m)]
    a = [i + j + [k] for i, j, k in zip(a, e, b)] + [r_1] + [r]
    return a, vect


def judge(matrix):
    if max(matrix[-2][:-1]) <= 0:
        flag = False
    else:
        flag = True
    return flag


def trans():
    global a, matrix, vect
    maxi = max(matrix[-2][:-1])
    index = a[-2][:-1].index(maxi)
    d = {}
    for i in a[:-2]:
        if i[index] > 0:
            d[i[-1] / i[index]] = a.index(i)
    pivot = d[min(d)]
    matrix[pivot] = matrix[pivot] / matrix[pivot][index]
    for i in range(m + 2):
        if i != pivot:
            matrix[i] = matrix[i] - matrix[i][index] * matrix[pivot]
    vect[pivot] = index + 1
    a = [list(i) for i in matrix]


def pr(matrix, vect):
    print('XB', end='\t\t')
    for i in range(n + m):
        print('X_{}'.format(i + 1), end='\t\t')
    print('b')
    for i in range(m + 2):
        if i <= m - 1:
            if vect[i] > m:
                print('(x_{})'.format(vect[i]), end='\t\t')
            else:
                print('x_{}'.format(vect[i]), end='\t\t')
        elif i == m:
            print('r_1', end='\t\t')
        elif i == m + 1:
            print('r', end='\t\t')
        for j in matrix[i]:
            print(f(str(j)).limit_denominator(), end='\t\t')  # 输出分数形式
        print(end='\n')


def print_solution(matrix):
    for i in range(n + m):
        print('x{}*={}'.format(i + 1, matrix[-2][i]), end='，')
    print('\nz*={}'.format(-matrix[-2][-1]))


def main():
    global a, matrix, vect
    a, vect = getinput()
    matrix = np.array(a, dtype=np.float64)
    pr(matrix, vect)
    while judge(matrix):
        trans()
        pr(matrix, vect)
    print_solution(matrix)


if __name__ == '__main__':
    main()
