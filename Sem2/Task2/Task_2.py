import numpy as np
import sympy as sp

from math import log, sqrt, pi, cos, ceil
from sympy.matrices import *

epsilon = 10e-07
tmp = (np.float32)(0.0)

def vector_norm(v):
    n = v.shape[0]

    w = np.zeros(n, dtype=np.float32)
    for i in range(n):
        w[i] = abs(v[i])

    return w.max()


def matrix_norm(A):
    (n,m) = A.shape

    if n != m:
        return -1
    else:
        summs = np.zeros(n, dtype=np.float32)
        for i in range(n):
            for j in range(n):
                summs[i] += abs(A[i, j])
        return summs.max()


def find_main_elem(A, k, det=1):
    row = A.shape[0]
    irow = k
    max = A[k, k]
    for i in range(k, row):
        if max < A[i, k]:
            max = A[i, k]
            irow = i

    A.row_swap(k, irow)
    if k != irow:
        det *= (-1)

    return det


def straight_step(A, k, det=1):
    row = A.shape[0]
    col = A.shape[1]
    tmp = A[k, k]
    det *= tmp

    if abs(tmp) < epsilon:
        print("Warning Element on {0} step is lesser than epsilon!".format(k + 1))
    for j in range(k, col):
        A[k, j] /= tmp

    for i in range(k + 1, row):
        tmp = A[i, k]
        for j in range(k, col):
            A[i, j] -= A[k, j] * tmp

    return det


def gauss_main_elem(A, det_A=1):
    row = A.shape[0]
    for k in range(0, row):
        det_A = find_main_elem(A, k, det_A)
        det_A = straight_step(A, k, det_A)

    return det_A


def reverse_step(A):
    row = A.shape[0]
    col = A.shape[1]
    x = sp.zeros(row, 1)

    for i in range(row - 1, -1, -1):
        summ = 0
        for j in range(i + 1, row):
            summ += A[i, j] * x[j]
        x[i] = A[i, col - 1] - summ

    return x


def check_x(A, x, b):
    v = A * x - b
    n = vector_norm(v)

    print("||Ax - b|| = ", n)
    return n < epsilon


def borders_Gersh(A):
    m = (np.float32)(max(A))
    M = (np.float32)(0.0)
    row = A.shape[0]

    for i in range(0, row):
        summ = (np.float32)(0.0)
        for j in range(0, row):
            if i!= j:
                summ += abs(A[i, j])
        if m > (A[i, i] - summ):
            m = (A[i, i] - summ)
        if M < (A[i, i] - summ):
            M = (A[i, i] + summ)

    return m, M

def odd_numbers(l):
    Q = []
    for i in range(l):
        Q.append([])

    Q[0] = [1]

    for i in range(1, l):
        Q[i] = np.zeros(2 ** i)

        for j in range(2 ** (i - 1)):
            Q[i][2 * j] = Q[i - 1][j]
            Q[i][2 * j + 1] = 2 ** (i + 1) - Q[i][2 * j]

    return Q[l // 4]


A = Matrix([[2.20219, 0.33266, 0.16768, 2.17523],
            [0.33266, 3.17137, 0.54055, 6.58335],
            [0.16768, 0.54055, 4.92343, 6.36904]], dtype = np.float32)


row = A.shape[0]
col = A.shape[1]


print("Matrix A|b:")
sp.pprint(A)

A_clear = A[:, 0 : col - 1]
b = A[:, col - 1]


#1 Gauss method

print('Single-fission scheme:')

det_A=1
A_1 = A.copy()
det_A = gauss_main_elem(A_1, det_A)
print("A in a triangular form: ")
sp.pprint(A_1)

x = reverse_step(A_1)

if(check_x(A_clear, x, b)):
    print("Solution: ")
    sp.pprint(x)

#2 Findin m and M and alpha

m, M = borders_Gersh(A)
alpha = (np.float32)(2 / (m + M))

print('alpha = ', alpha)

#3 System Conversion

E = Matrix(np.eye(row))
B = E - (alpha * A_clear)
c = alpha * b

print("B: ")
sp.pprint(B)

print("c: ")
sp.pprint(c)

B_norm = (np.float32)(matrix_norm(B))
print("||B|| = ", B_norm)

if B_norm >= 1:
    print("The method of simple iterations can not converge!")

#4 The a priori estimate of k

print("Enter epsilon: ")
epsilon = (np.float32)(input())

x_0 = Matrix(np.zeros(row, dtype = np.float32))

print("Initial approximation x_0: ")
sp.pprint(x_0)

x_1 = Matrix(np.array(B * x_0 + c))

k = (int)(log(epsilon * (1 - B_norm) / vector_norm(x_1 - x_0), B_norm))

print("A priori estimate of k: ", k)

print("||x* - x^N|| = ", vector_norm(x - x_1))

#5 Method of simple iteration with accuracy

x_0 = Matrix(np.zeros(row, dtype = np.float32))

x_1 = Matrix(np.array(B * x_0 + c))

first_apr = vector_norm(x_1 - x_0)

print("Step 1:")
apriori = (np.float32)((np.float32)(B_norm / (1 - B_norm)) * first_apr)
print("A priori: ", apriori)
aposteriori = (np.float32)((np.float32)(B_norm) / (np.float32)(1 - B_norm))
print("A posteriori: ", aposteriori  * first_apr)
print("||x* - x ^", 1, "|| = ", vector_norm(x - x_1), sep='')

k_iter = 1

while vector_norm(x_1 - x_0) > epsilon:
    k_iter += 1
    print("Step {0}: ".format(k_iter))
    x_0 = x_1.copy()
    x_1 = Matrix(np.array(B * x_0 + c))
    apriori = (np.float32)((np.float32)(B_norm ** k_iter/ (1 - B_norm)) * first_apr)
    print("A priori: ", apriori)
    aposteriori = (np.float32)((np.float32)(B_norm) / (np.float32)(1 - B_norm) * vector_norm(x_1 - x_0))
    print("A posteriori: ", aposteriori)
    print("||x* - x ^", k_iter, "|| = ", vector_norm(x - x_1), sep='')

print("The iterative method required {0} iterations".format(k_iter))

#6 Method of simple iterations with Chebyshev parameters

print("Enter epsilon: ")
epsilon = (np.float32)(input())

k = (int)(1 / 2 * sqrt(pi / m) * log(2 / epsilon))

p = ceil(log(k, 2))
k = 2 ** p

tetas = odd_numbers(k)

t = cos(pi * tetas[0] / (2 * k))
tao = (np.float32)(2 / (M + m - (M - m) * t))

x_0 = Matrix(np.zeros(row))
x_1 = x_0 + tao * (b - A_clear * x_0)

first_apr = vector_norm(x_1 - x_0)
k_iter = 1

print("Step 1:")
apriori = (np.float32)((np.float32)(B_norm / (1 - B_norm)) * first_apr)
print("A priori: ", apriori)
aposteriori = (np.float32)((np.float32)(B_norm) / (np.float32)(1 - B_norm))
print("A posteriori: ", aposteriori  * first_apr)
print("||x* - x ^", 1, "|| = ", vector_norm(x - x_1), sep='')


while vector_norm(x_1 - x_0) > epsilon:
    x_0 = x_1.copy()
    k_iter += 1
    t = cos(pi * tetas[k_iter - 1] / (2 * k))
    tao = (np.float32)(2 / (M + m - (M - m) * t))
    x_1 = x_0 + tao * (b - A_clear * x_0)

    print("Step {0}: ".format(k_iter))
    apriori = (np.float32)((np.float32)(B_norm ** k_iter / (1 - B_norm)) * first_apr)
    print("A priori: ", apriori)
    aposteriori = (np.float32)((np.float32)(B_norm) / (np.float32)(1 - B_norm) * vector_norm(x_1 - x_0))
    print("A posteriori: ", aposteriori)
    print("||x* - x^", k_iter, "|| = ", vector_norm(x - x_1), sep='')

print("The iterative method with Chebyshev parameters required {0} iterations".format(k_iter))

