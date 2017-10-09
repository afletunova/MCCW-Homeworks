import math
from operator import itemgetter

def f(x):
    return math.exp(x) - x ** 2 / 2

a = 0
b = 1
#n = 7
m = 15
#y = 0.65

print('Input n <= ', m, 'n > 0:')
n = int(input())
while n > m or n <= 0:
    print('Input n <= ', m, 'n > 0:')
    n = int(input())

print('[ a , b ] = [', a, ',', b,'] Input x:')
y = float(input())

points_values = []

for i in range(0, m + 1):
    x = a + i * (b - a) / m
    points_values.append((x, f(x)))
    print('x = ', x, ', f(x) = ', f(x))

point_dist = []

for i in range(0, m + 1):
    point_dist.append((points_values[i][0], abs(points_values[i][0] - y)))

point_dist = sorted(point_dist, key = itemgetter(1))

points = []
f_points = []
for i in range(0, n + 1):
    points.append(point_dist[i][0])

points = sorted(points)

for i in range(0, n + 1):
    f_points.append(f(points[i]))

#Newton conception
print('Newton: ')
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

print('P_n_x =', P_n_x)

print('ef_n(x) = ', abs(f(y) - P_n_x))

#Lanrange conception
print('Lagrange: ')
LP_n_x = 0
for i in range(0, n + 1):
    numerator = 1
    for j in range(0, n + 1):
        if i != j:
            numerator *= (y - points[j])
    denominator = 1
    for j in range(0, n + 1):
        if i != j:
            denominator *= points[i] - points[j]
    LP_n_x += ((numerator / denominator) * f_points[i])

print('P_n_x =', LP_n_x)

print('ef_n(x) = ', abs(f(y) - LP_n_x))
