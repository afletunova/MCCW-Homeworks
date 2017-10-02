import math


def f(x):
    a = 10
    b = -0.1
    return a * math.cos(x) + b * x ** 2


def f_der(x):
    a = -10
    b = -0.2
    return a * math.sin(x) + b * x


bound_A = -8
bound_B = 2
epsilon = 10 ** (-7)

print("A = ", bound_A, "; B = ", bound_B)
print("f(x) = 10 * cos(x) - 0.1 * x ^ 2")
print("e = ", epsilon)

h = (bound_B - bound_A) / 100


# The process of separation of roots
x_1 = bound_A
x_2 = x_1 + h

intervals = []
counter = 0

while x_2 <= bound_B:
    if f(x_1) * f(x_2) <= 0:
        intervals.append((x_1, x_2))
        print("x_1 = ", x_1, "; x_2 = ", x_2)
        counter += 1

    x_1 = x_2
    x_2 = x_1 + h
print("Number of roots: ", counter)

# Bisection method
print("Bisection method")
c = 0
for (x, y) in intervals:
    counter = 1
    c = (x + y) / 2
    print("Initial Approximation: ", c)
    while y - x > 2 * epsilon:
        if f(x) * f(c) < 0:
            y = c
        else:
            x = c
        c = (x + y) / 2
        counter += 1
    print("Step: ", counter, "; Root = ", c)
    print("|b_N - a_N| = ", math.fabs(y - x))
    print("|f(X_N) - 0| = ", math.fabs(f(c) - 0))


# Newton's method
print("Newton's method")
for (x, y) in intervals:
    a = (x + y) / 2
    b = a - f(a) / f_der(a)
    print("Initial Approximation: ", a)
    counter = 1
    while math.fabs(b - a) > epsilon:
        a = b
        b = a - f(a) / f_der(a)
        counter += 1
    print("Step: ", counter, "; Root = ", b)
    print("|X_N - X_(N-1)| = ", math.fabs(b - a))
    print("|f(X_N) - 0| = ", math.fabs(f(c) - 0))


# Modified Newton's method
print("Modified Newton's method")
for (x, y) in intervals:
    x_0 = a = (x + y) / 2
    b = a - f(a) / f_der(a)
    print("Initial Approximation: ", x_0)
    counter = 1
    while math.fabs(b - a) > epsilon:
        a = b
        b = a - f(a) / f_der(x_0)
        counter += 1
    print("Step: ", counter, "; Root = ", b)
    print("|X_N - X_(N-1)| = ", math.fabs(b - a))
    print("|f(X_N) - 0| = ", math.fabs(f(c) - 0))


# Method of secants
print("Method of secants")
for (x, y) in intervals:
    a = (x + y) / 2
    b = a - f(a) / f_der(a)
    print("Initial Approximations: ", a, b)
    counter = 1
    while math.fabs(b - a) > epsilon:
        a = b - (b - a) * f(b) / (f(b) - f(a))
        b = a + (a - b) * f(a) / (f(a) - f(b))
        counter += 1
    print("Step: ", counter, "; Root = ", b)
    print("|X_N - X_(N-1)| = ", math.fabs(b - a))
    print("|f(X_N) - 0| = ", math.fabs(f(c) - 0))

