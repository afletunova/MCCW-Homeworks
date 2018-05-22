import numpy as np
import sympy as sp
from sympy.matrices import *
import math

eps = (np.float32)(10e-10)

def vector_norm(v):
    n = v.shape[0]

    w = np.zeros(n, dtype=np.float32)
    for i in range(n):
        w[i] = abs(v[i])

    return w.max()


def scalar_vector_multi(alpha, v):
    n = v.shape[0]

    w = np.zeros(n)
    for i in range(n):
        w[i] = alpha * v[i]
    return w


def power_method(A):
    n = A.shape[0]
    y_0 = np.ones((n, 1), dtype=np.float32)
    y_1 = y_0.copy()
    lam = (np.float32)(0.0)
    count = 0
    k = 0
    while count < n:
        k += 1
        y_0 = y_1.copy() / vector_norm(y_1)
        y_1 = A * y_0
        lam = (np.float32)(y_1[0] / y_0[0])
        count = 0
        for i in range(n):
            count += 1
            abs_val = abs((np.float32)(y_1[i] / y_0[i]) - lam)
            if abs_val > eps:
                break

    return k, lam, y_1


def scalar_method(A):
    n = A.shape[0]
    y_0 = Matrix(np.ones((n, 1), dtype=np.float32))
    y_1 = y_0.copy()
    lam_0 = (np.float32)(0.0)
    lam_1 = eps + 1
    k = 0

    while abs(lam_1 - lam_0) > eps:
        lam_0 = lam_1
        k += 1
        y_0 = (y_1.copy() / vector_norm(y_1))
        y_1 = A * y_0
        dot_0 = y_1.dot(y_0)
        dot_1 = y_0.dot(y_0)
        lam_1 = (np.float32)(dot_0 / dot_1)

    return k, lam_1, y_1


def opp_spectr_bound(A, lam):
    n = A.shape[0]
    B = A - lam * np.eye(n)

    k, lam_B, x_1 = power_method(B)

    lam_A = (np.float32)(lam_B + lam)

    return lam_A

def Jakobi_method(A, epsi):
    n = A.shape[0]
    A_new = A.copy()
    k = 0
    U = np.eye(n)

    while True:
        k += 1
        i_max = -1
        j_max = -1
        a_max = -1
        for i in range(0, n):
            for j in range(i + 1, n):
                if a_max < abs(A_new[i, j]):
                    a_max = abs(A_new[i, j])
                    i_max = i
                    j_max = j

        print("Step: {}".format(k))
        print("A: ")
        sp.pprint(A_new)

        if a_max < epsi:
            break

        T = np.eye(n)
        d = math.sqrt((A_new[i_max, i_max] - A_new[j_max, j_max]) ** 2 + 4.0 * A_new[i_max, j_max] ** 2)
        c = math.sqrt(0.5 * (1.0 + abs(A_new[i_max, i_max] - A_new[j_max, j_max]) / d))
        s = np.sign(A_new[i_max, j_max] * (A_new[i_max, i_max] - A_new[j_max, j_max])) * math.sqrt(0.5 * (1.0 - (abs(A_new[i_max, i_max] - A_new[j_max, j_max]) / d)))

        T[i_max, i_max] = c
        T[j_max, j_max] = c
        T[i_max, j_max] = -s
        T[j_max, i_max] = s

        #for i in range(n):
        #    if i not in [i_max, j_max]:
        #        A_new[i, i_max] = A_new[i_max, i] = c * A_new[i, i_max] + s * A_new[i, j_max]
        #        A_new[i, j_max] = A_new[j_max, i] = (-1) * s * A_new[i, i_max] + c * A_new[i, j_max]
        #        A_new[i_max, i_max] = c ** 2 * A_new[i_max, i_max] + 2 * c * s * A_new[i_max, j_max] + s ** 2 * A_new[j_max, j_max]
        #        A_new[j_max, j_max] = s ** 2 * A_new[i_max, i_max] - 2 * c * s * A_new[i_max, j_max] + c ** 2 * A_new[j_max, j_max]
        #       A_new[i_max, j_max] = A_new[j_max, i_max] = 0

        A_new = A_new * T
        A_new = T.T * A_new
        U = U * T

    return A_new, U

A = Matrix([[-1.536984, -0.1999070, 0.958551],
            [-0.199070, 1.177416, 0.069925],
            [0.958551, 0.069925, -1.551506]], dtype = np.float32)

n = A.shape[0]
print('Power method:')
k, lam, x_1 = power_method(A)

sp.pprint(A)

print('N: {}, lambda_1 = {}'.format(k + 1, lam))
print('X_1 =')
sp.pprint(x_1)
print('||A * X_1 - lambda_1 * X_1|| =', vector_norm(A * x_1 - lam * x_1))
p_lam = lam

input('If you want to continue, press Enter:')

print('Scalar multiply method: ')
k, lam, x_1 = scalar_method(A)

print('N: {}, lambda_1 = {}'.format(k + 1, lam))
print('X_1 =')
sp.pprint(x_1)
print('||A * X_1 - lambda_1 * X_1|| =', vector_norm(A * x_1 - lam * x_1))

lam_n = opp_spectr_bound(A, p_lam)

print('Opposite spectrum boundary:', lam_n)

if n == 3:
    summ = 0
    for i in range(n):
        summ += A[i, i]
    lam_3 = summ - p_lam - lam_n
    print('Another eigenvalue:', lam_3)

epsi = (np.float32)(10e-05)
Lam, U = Jakobi_method(A, epsi)

print('lambda_1 = {}, lambda_2 = {}, lambda_3 = {}'.format(Lam[0, 0], Lam[1, 1], Lam[2, 2]))
print('x_1:')
sp.pprint(Matrix(U[:, 0]))
print('x_2:')
sp.pprint(Matrix(U[:, 1]))
print('x_3:')
sp.pprint(Matrix(U[:, 2]))