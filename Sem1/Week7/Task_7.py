import math

def f(x, y):
    return -y + math.cos(x)

y = lambda x: 1 / 2 * (math.exp(-x) + math.sin(x) + math.cos(x))

h = 0.1
N = 10
x_0 = 0
y_0 = 1

derivatives = []
derivatives.append(-f(x_0, y_0) - math.sin(x_0))
derivatives.append(-derivatives[x_0] - math.cos(x_0))
derivatives.append(math.sin(x_0) - derivatives[1])
derivatives.append(math.cos(x_0) - derivatives[2])
derivatives.append(- derivatives[3] - math.sin(x_0))
derivatives.append(- derivatives[4] - math.cos(x_0))
derivatives.append(- derivatives[5] + math.sin(x_0))
derivatives.append(- derivatives[6] + math.cos(x_0))
derivatives.append(- derivatives[7] - math.sin(x_0))
derivatives.append(- derivatives[8] - math.cos(x_0))
derivatives.append(- derivatives[9] + math.sin(x_0))
derivatives.append(- derivatives[10] + math.cos(x_0))

# 1
points = []
y_x = []
y_N_x = []

for k in range(-2, N + 1):
    points.append(x_0 + k * h)
    y_x.append(y(x_0 + k * h))

print('{0:20} | {1:20}'.format('x', 'y(x)'))
for t, p in zip(points, y_x):
    print('{0:<20.15f} | {1:<20.15f}'.format(t, p))

answer = 'yt'
point = 0

while answer != 'n':
    if answer == 'n':
        break
    while answer not in ('n', 'y'):
        print('If you ready to start, enter y. Else - n to exit: ')
        answer = input()
    choice = -1
    while choice not in (0, 1, 2, 3, 4, 5):
        print('Choose method: ')
        print('1 - The Taylor series expansion method;')
        if point == 1:
            print('2 - The Adams Method;')
        print('3 - The Runge-Kutta Method;')
        print('4 - The Eulers Method;')
        print('5 - The Improved Euler Method')
        print('6 - The Euler-Cauchy method')
        print('   or enter 0 to exit:')
        choice = int(input())

        if choice == 0:
            break
        elif choice == 1:
        #2
            for i in range(0, 5):
                y_N = y_0 + f(x_0, y_0) * (points[i] - x_0)
                counter = 0
                values = []
                for j in range(0, 9):
                    if derivatives[j] / math.factorial(j + 2) * (points[i] - x_0) ** (j + 2) != 0:
                        counter += 1
                    y_N += derivatives[j] / math.factorial(j + 2) * (points[i] - x_0) ** (j + 2)
                values.append(y_N)
                values.append(counter)
                y_N_x.append(values)

            print('{0:20} | {1:20}'.format('x', 'y(x)'))
            for t,p in zip(points, y_N_x):
                print('{0:<20.15f} | {1:<20.15f}'.format(t, p[0]))

        #3
            print('{0:20} | {1:20}'.format('x', '|y(x)-y_N(x)|'))
            for i, j, k in zip(points, y_x, y_N_x):
                print('{0:<20.15f} | {1:<20.15f}'.format(i, abs(k[0] - j)))

            point = 1

        elif choice == 2:
            fin_diff = [[0 for x in range(7)] for y in range(N + 3)]
            for i in range(0, N + 1):
                fin_diff[i][0] = points[i]
            for i in range(0, 5):
                fin_diff[i][1] = y_N_x[i][0]
                fin_diff[i][2] = h * f(points[i], y_N_x[i][0])
            for k in range(5, N + 3):
                for i in range(3, 7):
                    for j in range(0, k - i + 3):
                        fin_diff[j][i] = fin_diff[j + 1][i - 1] - fin_diff[j][i - 1]
                ta_1 = fin_diff[k - 1][1]
                ta_2 = fin_diff[k - 1][2]
                ta_3 = 1 / 2 * fin_diff[k - 2][3]
                ta_4 = 5 / 12 * fin_diff[k - 3][4]
                ta_5 = 3 / 8 * fin_diff[k - 4][5]
                ta_6 = 251 / 720 * fin_diff[k - 5][6]
                fin_diff[k][1] = ta_1 + ta_2 + ta_3 + ta_4 + ta_5 + ta_6
                fin_diff[k][2] = h * f(fin_diff[k][0], fin_diff[k][1])

            print('{0:20} | {1:20}'.format('x', 'y(x)'))
            for i in range(0, N + 1):
                print('{0:<20.15f} | {1:<20.15f}'.format(fin_diff[i][0], fin_diff[i][1]))

            print('|y_N(x) - y(x)| =', abs(fin_diff[N][1] - y(points[N])))


        elif choice == 3:
        #5
            points = []
            points.append(x_0)
            for i in range(1, N + 1):
                points.append(x_0 + h * i)

            y_points = []
            y_points.append(y_0)
            for k in range(1, N + 1):
                k_1 = h * f(points[k - 1], y_points[k - 1])
                k_2 = h * f(points[k - 1] + h / 2, y_points[k - 1] + k_1 / 2)
                k_3 = h * f(points[k - 1] + h / 2, y_points[k - 1] + k_2 / 2)
                k_4 = h * f(points[k - 1] + h, y_points[k - 1] + k_3)
                y_points.append(y_points[k - 1] + 1 / 6 *(k_1 + 2 * k_2 + 2 * k_3 + k_4))

            print('{0:20} | {1:20}'.format('x', 'y_N(x)'))
            for i, j in zip(points, y_points):
                print('{0:<20.15f} | {1:<20.15f}'.format(i, j))

            print('|y_N(x) - y(x)| =', abs(y_points[N] - y(points[N])))

        #6
        elif choice == 4:

            points = []
            points.append(x_0)
            for i in range(1, N + 1):
                points.append(x_0 + h * i)

            y_points = []
            y_points.append(y_0)
            for k in range(1, N + 1):
                y_points.append(y_points[k - 1] + h * f(points[k - 1], y_points[k - 1]))

            print('{0:20} | {1:20}'.format('x', 'y_N(x)'))
            for i, j in zip(points, y_points):
                print('{0:<20.15f} | {1:<20.15f}'.format(i, j))

            print('|y_N(x) - y(x)| =', abs(y_points[N] - y(points[N])))

        elif choice == 5:
            points = []
            points.append(x_0)
            for i in range(1, N + 1):
                points.append(x_0 + h * i)

            y_points = []
            y_points.append(y_0)
            for k in range(1, N + 1):
                y_k_12 = y_points[k - 1] + h / 2 * f(points[k - 1], y_points[k - 1])
                y_points.append(y_points[k - 1] + h * f(points[k - 1] + h / 2, y_k_12))

            print('{0:20} | {1:20}'.format('x', 'y_N(x)'))
            for i, j in zip(points, y_points):
                print('{0:<20.15f} | {1:<20.15f}'.format(i, j))

            print('|y_N(x) - y(x)| =', abs(y_points[N] - y(points[N])))

        elif choice == 6:
            points = []
            points.append(x_0)
            for i in range(1, N + 1):
                points.append(x_0 + h * i)

            y_points = []
            y_points.append(y_0)
            for k in range(1, N + 1):
                y_k_1 = y_points[k - 1] + h * f(points[k - 1], y_points[k - 1])
                y_points.append(y_points[k - 1] + h / 2 * (f(points[k - 1] , y_points[k - 1]) + f(points[k], y_k_1)))

            print('{0:20} | {1:20}'.format('x', 'y_N(x)'))
            for i, j in zip(points, y_points):
                print('{0:<20.15f} | {1:<20.15f}'.format(i, j))

            print('|y_N(x) - y(x)| =', abs(y_points[N] - y(points[N])))

    if choice == 0:
        break
