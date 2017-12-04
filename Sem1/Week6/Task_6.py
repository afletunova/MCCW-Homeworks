import math
import scipy.integrate as ig


def q_formula(points, m, f):
    summ = 0
    for i in range(0, m):
        summ += f(points[i]) + f(points[i + 1]) + 4 *f((points[i] + points[i + 1]) / 2)
    summ *= h / 6
    return summ


def solver(f, h, epsilon):
    x_1 = a
    x_2 = x_1 + h
    intervals = []
    roots = []

    while x_2 <= b:
        if w_2(x_1) * w_2(x_2) <= 0:
            intervals.append((x_1, x_2))

        x_1 = x_2
        x_2 = x_1 + h

    c = 0
    for (x, y) in intervals:
        c = (x + y) / 2
        while y - x > 2 * epsilon:
            if f(x) * f(c) < 0:
                y = c
            else:
                x = c
            c = (x + y) / 2
        roots.append(c)
    return roots


def system_solver(a_coef, b_coef):
    delta = a_coef[0] * a_coef[3] - a_coef[1] * a_coef[2]
    delta_1 = b_coef[0] * a_coef[3] - b_coef[1] * a_coef[1]
    delta_2 = b_coef[1] * a_coef[0] - b_coef[0] * a_coef[2]
    c = delta_1 / delta
    d = delta_2 / delta
    return c, d

#1
a = 0
b = 1
N = 2
m = 100
epsilon = 10 ** (-7)

#f = lambda x: math.sin(x)
f = lambda  x: x ** 2 + 7 * x


w = lambda x: abs(x - 1 / 2)

g = lambda x: w(x) * f(x)

J = ig.quad(g, a, b)[0]

h = (b - a) / m

points = []
for k in range (0, m + 1):
    points.append(a + k * h)

#2
t_values = []
t_values.append(-math.sqrt(3) / 3)
t_values.append(math.sqrt(3) / 3)

sum_gauss = 0
for k in range(0, m):
    x_values = []
    for j in range(0, 2):
        x_values.append(h / 2 * t_values[j] + a + k * h + h / 2)
    sum_gauss += (h / 2 * (g(x_values[0]) + g(x_values[1])))

print('Integral w(x)f(x) by [A,B] =', sum_gauss)
print('|J - J_f| =', abs(J - sum_gauss))

answer = 0
while answer != 1:
    print('Enter 1 to continue:')
    answer = int(input())
#3
mu = []
for i in range(0, 4):
    p = lambda x: w(x) * x ** i
    mu.append(q_formula(points, m, p))
    print('Mu_', i, ' = ', mu[i], sep='')

a_coef = [mu[1], mu[0], mu[2], mu[1]]
b_coef = [-mu[2], -mu[3]]
c, d = system_solver(a_coef, b_coef)

w_2 = lambda  x: x ** 2 + c * x + d
print('w_2(x) = x^2 +', round(c, 3), '* x +', round(d, 3))

roots = solver(w_2, h, epsilon)
print('w_2 roots:', roots[0], ',', roots[1])

a_coef = [1, 1, roots[0], roots[1]]
b_coef = [mu[0], mu[1]]

A_1, A_2 = system_solver(a_coef, b_coef)
print('A_1 =', A_1, 'A_2 =', A_2)

print('|mu_3 - A_1 * (x_1) ^ 3 - A_2 * (x_2) ^ 3| =', abs(mu[3] - A_1 * roots[0] ** 3 - A_2 * roots[1] ** 3))

J_f = A_1 * f(roots[0]) + A_2 * f(roots[1])

print('Integral w(x)f(x) by [A,B] =', J_f)
print('|J - J_f| =', abs(J - J_f))