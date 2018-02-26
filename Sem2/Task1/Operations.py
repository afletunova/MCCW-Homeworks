#Exact methods for solving systems of linear algebraic equations. The Gauss method

from decimal import *
import  numpy as np

def matrix_summ(A, B):
    (n, l) = A.shape
    (k, m) = B.shape

    if (n != k) or (l != m):
        return -1
    else:
        C = np.zeros((n, l))
        for i in range(n):
            for j in range(l):
                C[i][j] = A[i][j] + B[i][j]
        return C


def matrix_comp(A, B):
    (n, l) = A.shape
    (k, m) = B.shape

    if l != k:
        return -1
    else:
        C = np.zeros((n, m))
        for i in range(n):
            for j in range(m):
                for r in range(l):
                    C[i][j] += A[i][r] * B[r][j]
        return C


def matrix_norm(A):
    (n,m) = A.shape

    if n != m:
        return -1
    else:
        summs = np.zeros(n)
        for i in range(n):
            for j in range(n):
                summs[i] += abs(A[i][j])
        return summs.max()


def matrix_tran(A):
    (n,m) = A.shape

    C = np.zeros((m,n))

    for i in range(n):
        for j in range(m):
            C[j][i] = A[i][j]

    return C


def matrix_vector_comp(A, b):
    (n,m) = A.shape
    k = b.shape[0]

    if m != k:
        return -1
    else:
        x = np.zeros(n)
        for i in range(n):
            for j in range(m):
                x[i] += A[i][j] * b[j]
        return x

def vector_matrix_comp(b, A):
    (n,m) = A.shape
    k = b.shape[0]

    if n != k:
        return -1
    else:
        x = np.zeros(m)
        for j in range(m):
            for i in range(n):
                x[j] += b[i] * A[i][j]
        return x


def scalar_matrix_multi(alpha, A):
    (n, m) = A.shape

    C = np.zeros((n,m))
    for i in range(n):
        for j in range(m):
            C[i][j] = alpha * A[i][j]

    return C


def vectors_summ(v, u):
    n = v.shape[0]
    m = u.shape[0]

    if n != m:
        return -1
    else:
        w = np.zeros(n)
        for i in range(n):
            w[i] = v[i] + u[i]

        return w


def vectors_multi(v, u):
    n = v.shape[0]
    m = u.shape[0]

    if n != m:
        return -1
    else:
        w = np.zeros(n)
        for i in range(n):
            w[i] = v[i] * u[i]

        return w


def scalar_vector_multi(alpha, v):
    n = v.shape[0]

    w = np.zeros(n)
    for i in range(n):
        w[i] = alpha * v[i]
    return w


def vector_norm(v):
    n = v.shape[0]

    w = np.zeros(n)
    for i in range(n):
        w[i] = abs(v[i])

    return w.max()


def pretty_print(A):
    np.set_printoptions(formatter={'float': '{: 0.8f}'.format})
    print(A)
