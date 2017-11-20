import math
import scipy.integrate as ig

#1
A = 0
B = 1
m = 10

w = lambda x: 1
f = lambda x: math.exp(-x) - x ** 2 / 2

#2
J = ig.quad(f, A, B)[0]

print('Integral of w(x)f(x) over [A,B]:', J)

#3 and 4
h = (B - A) / m

points = []
for k in range(0, m + 1):
    points.append(A + k * h)

com_qf_m = 0
for i in range(0, m):
    com_qf_m += f((points[i] + points[i + 1]) / 2)

com_qf_m *= h
print('Composite Quadratic formula of medium rectangles:', com_qf_m)
print('|J - J(h)| =', abs(J - com_qf_m))

com_qf_tr = h / 2 * (f(A) + f(B))
sum = 0
for i in range(1, m):
    sum += f(points[i])

com_qf_tr += h * sum
print('Composite Quadratic formula of trapezium:', com_qf_tr)
print('|J - J(h)| =', abs(J - com_qf_tr))

com_qf_simp = 0
for i in range(0, m):
    com_qf_simp += f(points[i]) + f(points[i + 1]) + 4 * f((points[i] + points[i + 1]) / 2)

com_qf_simp *= h / 6
print('Composite Quadratic Simpson formula:', com_qf_simp)
print('|J - J(h)| =', abs(J - com_qf_simp))
