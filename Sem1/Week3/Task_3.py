import math

def f(x):
    return 1 - math.exp(-x) + x ** 2


a = 0
b = 1
m = 10

h = (b - a) / m

points = []
f_points = []

for i in range(0, m + 1):
    x = a + i * h
    points.append(x)
    f_points.append(f(x))
    print('x = ', x, ', f(x) = ', f(x))

n = 0
print('Input n <=', m, 'n > 0:')
n = int(input())
while n > m or n <= 0:
    print('Input n <=', m, 'n > 0:')
    n = int(input())

fin_diff = [points, f_points]
for i in range(2, n + 2):
    column = []
    for j in range(0, m + 2 - i):
        column.append(fin_diff[i - 1][j + 1] - fin_diff[i - 1][j])
    fin_diff.append(column)

x = a - b
a_1 = a + ((n + 1) // 2) * h
b_1 = b - ((n + 1) // 2) * h
print ('x: ', x, 'a_1: ', a_1, 'b_1: ', b_1)
while not ((a <= x <= a + h) or (b - h <= x <= b) or (a_1 <= x <= b_1)):
    print('Input x in', '[', a, ',', a + h, '] or [', b - h, ',', b, '] or [', a + ((n + 1) // 2) * h, ',', b - ((n + 1) // 2) * h, ']')
    x = float(input())

if a <= x <= a + h:
    P_n_x = fin_diff[1][0]
    t = (x - a) / h
    for i in range (0, n):
        tail = fin_diff[i + 2][0] / math.factorial(i + 1) * t
        for j in range(1, i + 1):
            tail *= (t - j)
        P_n_x += tail
elif b - h <= x <= b:
    P_n_x = fin_diff[1][m]
    t = (x - b) / h
    for i in range(0, n):
        tail = fin_diff[i + 2][m - (i + 1)] / math.factorial(i + 1) * t
        for j in range(1, i + 1):
            tail *= (t + j)
        P_n_x += tail
else:
    left_n = int((x - a) // h)
    P_n_x = fin_diff[1][left_n]
    t = x - (fin_diff[0][left_n]) / h
    for i in range(0, n):
        tail = fin_diff[i + 2][left_n - (i + 1) // 2] / math.factorial(i + 1)
        for j in range(1, i + 1):
            tail *=(t + ((-1) ** j * (j + 1)) // 2)
        P_n_x += tail

print('P_n(x) =', P_n_x)
print('ef_n(x) = ', abs(f(x) - P_n_x))