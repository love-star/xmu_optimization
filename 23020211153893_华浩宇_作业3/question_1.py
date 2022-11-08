import numpy as np


def simplex(c, A, b):
    if c.shape[0] != A.shape[1]:
        print("A和C形状不匹配")
        return 0
    if b.shape[0] != A.shape[0]:
        print("A和b形状不匹配")
        return 0
    num = A.shape[1] - A.shape[0]
    N_indexs = np.arange(0, num)
    B_indexs = np.arange(num, A.shape[1])
    N = A[:, N_indexs]
    B = A[:, B_indexs]
    c_N = c[N_indexs, :]
    c_B = c[B_indexs, :]
    x_B = np.matmul(np.linalg.inv(B), b)
    j = 0
    while True:
        R = (c_N.T - np.matmul(np.matmul(c_B.T, np.linalg.inv(B)), N)).flatten()
        if all(R >= 0):
            print("最优点：")
            print("非基变量的索引：", N_indexs)
            print("基变量的索引：", B_indexs)
            print("最优点：", x_B.flatten())
            print("最优值：", -(z.flatten()[0]))
            return 0
        else:
            j = j + 1
            print("步骤", j)
            N_in_index = np.argmin(R)
            N_in = N[:, N_in_index]
            negd = np.matmul(np.linalg.inv(B), N_in)
            with np.errstate(divide='ignore'):
                y = x_B.flatten() / negd
            for i in range(len(y)):
                if negd[i] > 0:
                    index = i
                    break
            for i in range(len(y)):
                if negd[i] <= 0:
                    continue
                elif y[i] < y[index]:
                    index = i
            B_out_index = index
            B_out = B[:, B_out_index]
            in_index = N_indexs[N_in_index]
            out_index = B_indexs[B_out_index]
            temp = N_indexs[N_in_index]
            N_indexs[N_in_index] = B_indexs[B_out_index]
            B_indexs[B_out_index] = temp
            N[:, N_in_index] = B_out
            B[:, B_out_index] = N_in
            c_N = c[N_indexs, :]
            c_B = c[B_indexs, :]
            N = A[:, N_indexs]
            B = A[:, B_indexs]
            x_B = np.matmul(np.linalg.inv(B), b)
            z = np.matmul(c_B.T, x_B)
            print("出基变量索引：", out_index)
            print("入基变量索引：", in_index)
            print("值：", z.flatten()[0])


if __name__ == '__main__':
    c = np.array([-2, -5, 0, 0, 0]).reshape(-1, 1)
    A = np.array([[1, 0, 0, 1, 0], [0, 1, 0, 0, 1], [1, 1, 1, 0, 0]])
    b = np.array([4, 6, 8]).reshape(-1, 1)
    simplex(c, A, b)
