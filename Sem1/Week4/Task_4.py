import math
from operator import itemgetter

a = 0
b = 1
m = 10

h = (b - a) / m

def f(x):
    return math.exp(-x) - x ** 2 / 2

def f_der(x):
    return -x - math.exp(-x)

def s_der(x):
    return  math.exp(-x) - 1

def P_n(y, n):
    points_values = []

    for i in range(0, m + 1):
        x = a + i * (b - a) / m
        points_values.append((x, f(x)))

    point_dist = []

    for i in range(0, m + 1):
        point_dist.append((points_values[i][0], abs(points_values[i][0] - y)))

    point_dist = sorted(point_dist, key=itemgetter(1))

    points = []
    f_points = []
    for i in range(0, n + 1):
        points.append(point_dist[i][0])

    points = sorted(points)

    for i in range(0, n + 1):
        f_points.append(f(points[i]))

    split_diff_table = [points, f_points]

    for i in range(2, n + 2):
        row = []
        for j in range(0, n + 2 - i):
            row.append((split_diff_table[i - 1][j + 1] - split_diff_table[i - 1][j]) /
                       (split_diff_table[0][j + i - 1] - split_diff_table[0][j]))
        split_diff_table.append(row)

    P_n_x = split_diff_table[1][0]
    for i in range(1, n + 1):
        tail = split_diff_table[1 + i][0]
        for j in range(0, i):
            tail *= (y - points[j])
        P_n_x += tail

    return P_n_x

points = []
f_points = []

for i in range(0, m + 1):
    x = a + i * h
    points.append(x)
    f_points.append(f(x))
    print('x = ', x, ', f(x) = ', f(x))

#First method
print('First method')

f_min = min(f(a), f(b))
total_min = min(f(a), f(b))
f_max = max(f(a), f(b))
total_max = max(f(a), f(b))

f_inv_points = []
for i in range(0, m + 1):
    x = a + i * h
    f_inv_points.append(x)
    if (f_points[i] < f_min):
        f_min = f_points[i]
    if (f_points[i] > f_max):
        f_max = f_points[i]
    print('f(x) = ', f(x), ', f^(-1)(x) = ', x)

print('Input n <=', m, 'n > 0:')
n = int(input())
while n > m or n <= 0:
    print('Input n <=', m, 'n > 0:')
    n = int(input())

print('[ c , d ] = [', f_min, ',', f_max,'] Input F:')
f_v = float(input())
while f_v < f_min or f_v > f_max:
    print('[ c , d ] = [', f_min, ',', f_max, '] Input F:')
    f_v = float(input())

choice = 0
while choice not in [1,2]:
    if f_min == total_min and f_max == total_max:
        print('Enter 1 to select the first method, 2 - second method: ')
        choice = int(input())
    else:
        choice = 2

if choice == 1:
    point_dist = []

    for i in range(0, m + 1):
        point_dist.append((f_points[i], abs(f_points[i] - f_v)))

    point_dist = sorted(point_dist, key = itemgetter(1))

    split_diff_table = [f_points, f_inv_points]

    for i in range(2, n + 2):
        row = []
        for j in range(0, n + 2 - i):
            row.append((split_diff_table[i - 1][j + 1] - split_diff_table[i - 1][j]) /
                        (split_diff_table[0][j + i - 1] - split_diff_table[0][j]))
        split_diff_table.append(row)

    for i in range(0, n):
        print(f_points[i], f_inv_points[i])


    Q_n_f = split_diff_table[1][0]
    for i in range(1, n + 1):
        tail = split_diff_table[1 + i][0]
        for j in range(0, i):
            tail *= (f_v - f_points[j])
        Q_n_f += tail

    print('Q_n(f) = x = ', Q_n_f)
    print('|f(x) - F| = ', abs(f(Q_n_f) - f_v))
elif choice == 2:
    x_1 = a
    x_2 = x_1 + h
    epsilon = 10 ** (-8)

    intervals = []

    while x_2 <= b:
        if (P_n(x_1, n) - f_v) * (P_n(x_2, n) - f_v) <= 0:
            intervals.append((x_1, x_2))
        x_1 = x_2
        x_2 = x_1 + h

    c = 0
    for (x, y) in intervals:
        c = (x + y) / 2
        while y - x > 2 * epsilon:
            if (P_n(x, n) - f_v) * (P_n(c, n) - f_v) < 0:
                y = c
            else:
                x = c
            c = (x + y) / 2
        print('Root = ', c)
    print('r_n(x) = |f(x) - F| =', abs(f(c) - f_v))

#Numerical differentiation
first_der_f = []
for i in range(0, m + 1):
    if i == 0:
        first_der_f.append((f_points[i + 1] - f_points[i]) / h)
    elif i == m:
        first_der_f.append((f_points[i] - f_points[i - 1]) / h)
    elif i == 1 or i == (m - 1):
        first_der_f.append((f_points[i + 1] - f_points[i - 1]) / (2 * h))
    elif i <= ((m + 1) / 2):
        first_der_f.append(((-3) * f_points[i] + 4 * f_points[i + 1] - f_points[i + 2]) / (2 * h))
    else:
        first_der_f.append((3 * f_points[i] - 4 * f_points[i - 1] + f_points[i - 2]) / (2 * h))

second_der_f = []
for i in range(1, m):
    second_der_f.append((f_points[i + 1] - 2 * f_points[i] + f_points[i - 1]) / h ** 2)

for i in range (0, m + 1):
    if i in [0, m]:
        print('x_', i, ' = ', points[i], ', f(x_',i, ') = ', f_points[i], 'f\'(x_', i, ') = ', first_der_f[i],
              ', |f\'(x_', i, ')_т - f\'(x_', i, ')_чд| = ', abs(f_der(points[i]) - first_der_f[i]), sep='')
    else:
        print('x_', i, ' = ', points[i], ', f(x_',i, ') = ', f_points[i], 'f\'(x_', i, ') = ', first_der_f[i],
              ', |f\'(x_', i, ')_т - f\'(x_', i, ')_чд| = ', abs(f_der(points[i]) - first_der_f[i]),
              ', f\"(x_',i, ') = ', second_der_f[i - 1], ', |f\"(x_', i, ')_т - f\"(x_', i,
              ')_чд| = ', abs(s_der(points[i]) - second_der_f[i - 1]), sep='')