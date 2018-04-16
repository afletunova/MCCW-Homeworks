import numpy as np
import sympy as sp
from sympy.matrices import *

epsilon = 10e-07

def vector_norm(v):
    n = v.shape[0]

    w = np.zeros(n)
    for i in range(n):
        w[i] = abs(v[i])

    return w.max()


def matrix_norm(A):
    (n,m) = A.shape

    if n != m:
        return -1
    else:
        summs = np.zeros(n)
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


def find_not_null_elem(A, k, det=1):
    row = A.shape[0]
    irow = k
    for i in range(k, row):
        if A[i, k] != 0:
            irow = i
            break

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


def gauss(A, det_A=1):
    row = A.shape[0]
    for k in range(0, row):
        det_A = find_not_null_elem(A, k, det_A)
        det_A = straight_step(A, k, det_A)

    return det_A


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


A = Matrix([[6.8704e-06, -7.3592e-03, 4.34112, 3.90080],
            [6.2704e-03, -0.73592, 1.57112, 1.70505],
            [0.90704, 0.86408, 1.02112, 2.04310]])


row = A.shape[0]
col = A.shape[1]


print("Matrix A|b:")
sp.pprint(A)

A_clear = A[:, 0 : col - 1]
b = A[:, col - 1]

#1 Gauss

print('Single-fission scheme:')

det_A=1
A_1 = A.copy()
det_A = gauss(A_1, det_A)

print("A in a triangular form: ")
sp.pprint(A_1)

x = reverse_step(A_1)

if(check_x(A_clear, x, b)):
    print("Solution: ")
    sp.pprint(x)

#print("det(A):", det_A)

#2 Gauss

print('Gauss method with the selection of the main element:')

det_A=1
A_1 = A.copy()
det_A = gauss_main_elem(A_1, det_A)
print("A in a triangular form: ")
sp.pprint(A_1)

x = reverse_step(A_1)

if(check_x(A_clear, x, b)):
    print("Solution: ")
    sp.pprint(x)


#3 A determinator

print("det(A):", det_A)


#4 The inverse matrix by the Gauss method with the choice

A_inv = np.zeros((row, row))
for i in range (row):
    A_clear = A[:, 0: col - 1]
    e = np.zeros(row)
    e[i] = 1
    A_clear = Matrix(np.column_stack((A_clear, e)))
    print(i)
    sp.pprint(A_clear)
    gauss_main_elem(A_clear)
    z = reverse_step(A_clear)
    for j in range(row):
        A_inv[j][i] = z[j]

print('A inverse:')
sp.pprint(A_inv)

#5 Searching for a solution using an inverse matrix

x = A_inv * b
A_clear = A[:, 0: col - 1]

if(check_x(A_clear, x, b)):
    print("Solution: ")
    sp.pprint(x)

#Condition number

print('Condition number of matrix A:')
mu = matrix_norm(A_clear) * matrix_norm(A_inv)
print('mu(A) = ||A|| * ||A^(-1)|| =', mu)

