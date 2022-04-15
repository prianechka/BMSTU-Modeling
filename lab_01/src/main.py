import matplotlib.pyplot as plt
import numpy as np

EPS = 10**(-2)

def f(x, u):
    return x**2 + u**2

def f_1(x):
    return x**3 / 3

def f_2(x, f1):
    return f1 + x**7 / 63

def f_3(x, f2):
    return f2 + 2*x**11 / 2079 + x**15 / 59535

def f_4(x, f3):
    return f3 + 2*x**15 / 93555 + 2*x**19 / 3393495 + \
        + 2*x**19 / 2488563 + 2*x**23 / 86266215 + 2*x**23 / 99411543 + \
            + 2*x**27 / 3341878155 + x**31 / 109876902975

def pikkar_method(x_min, x_max, h):
    x = x_min
    result_1 = []
    result_2 = []
    result_3 = []
    result_4 = []
    while x <= x_max:
        f1 = f_1(x)
        f2 = f_2(x, f1)
        f3 = f_3(x, f2)
        f4 = f_4(x, f3)
        result_1.append(f1)
        result_2.append(f2)
        result_3.append(f3)
        result_4.append(f4)
        x += h
    return result_1, result_2, result_3, result_4

# Alpha = 1
def runge_kutt_1(x_min, x_max, h, y0):
    x = x_min
    y = y0
    result_y = []
    result_y.append(y)
    while x <= x_max:
        try:
            y_temp = y + h/2 * f(x, y)
            y_pr_temp = f(x + h/2, y_temp)
            y += h * y_pr_temp
            x += h
            result_y.append(y)
        except:
            while x <= x_max:
                result_y.append('-')
                x += h
            return result_y
    
    return result_y

# Alpha = 0.5
def runge_kutt_2(x_min, x_max, h, y0):
    x = x_min
    y = y0
    result_y = []
    result_y.append(y)
    while x <= x_max:
        try:
            y_temp = y + h/2 * f(x, y)
            y_pr_temp = f(x + h, y_temp)
            y_pr_temp_2 = 0.5 * (f(x, y) + y_pr_temp)
            y += h * y_pr_temp_2
            x += h
            result_y.append(y)
        except:
            while x <= x_max:
                result_y.append('-')
                x += h
            return result_y
    
    return result_y

def euler(x_min, x_max, h, y0):
    x = x_min
    y = y0
    result = []
    result.append(y)
    while x <= x_max:
        try:
            y += h * f(x, y)
            x += h
            result.append(y)
        except:
            while x <= x_max:
                result.append('-')
                x += h
            return result

    return result

def print_result(x):
    print("|{:^14.8f}".format(x), end = "")

def print_x(x):
    print("|{:^14.2f}".format(x), end = "")

def print_str():
    print("|" + "—" * 119 + "|")

def print_header():
    print("|{:^14}|{:^14}|{:^14}|{:^14}|{:^14}|{:^14}|{:^14}|{:^14}|".format("X", \
        "Pikkar_1", "Pikkar_2", "Pikkar_3", "Pikkar_4", "Runge(a=1)", \
        "Runge(a=0.5)", "Euler"))

def find_accuracy(x_min, x_max, h, y):
    res5 = runge_kutt_1(x_min, x_max, h, y)
    res6 = runge_kutt_2(x_min, x_max, h, y)
    res7 = euler(x_min, x_max, h, y)

    x = 0
    N = len(res5)
    for i in range(0, N):
        if (abs(res5[i] - res6[i]) <= EPS) and \
            (abs(res6[i] - res7[i]) <= EPS) and \
                (abs(res5[i] - res7[i]) <= EPS):
            x +=  h
        else:
            print()
            print("Xmax = ", x)
            break
    return x

def draw_graph(x_max, h):
    X = np.arange(0, x_max + h, h)
    y0 = np.array(runge_kutt_1(0, x_max, h, 0))
    y1 = y0
    X1 = X * (-1)

    plt.figure()
    plt.xlabel("Ось X")
    plt.ylabel("Значение функции")
    plt.title("График функции")
    plt.plot(X, y0, color = "blue")
    plt.plot(X1, y1, color = "blue")
    plt.show()

def find_pikkar_accuracy(x_min, x_max, h, y):
    res1, res2, res3, res4 = pikkar_method(x_min, x_max, h)
    x = x_min
    N = len(res1)
    for i in range(0, N):
        if (abs(res1[i] - res2[i]) <= EPS):
            x +=  h
        else:
            print("Правая граница точности для Пиккара 1-го порядка: ", x)
            break

    x = x_min
    for i in range(0, N):
        if (abs(res2[i] - res3[i]) <= EPS):
            x += h
        else:
            print("Правая граница точности для Пиккара 2-го порядка: ", x)
            break

    x = x_min
    for i in range(0, N):
        if (abs(res3[i] - res4[i]) <= EPS):
            x +=  h
        else:
            print("Правая граница точности для Пиккара 3-го порядка: ", x)
            break
    
    

x_min = 0
x_max = 2.02
h = 10**(-5)
y = 0

find_pikkar_accuracy(x_min, x_max, h, y)
#draw_graph(x_max, h)

'''

# find_accuracy(x_min, x_max, h, y)


res1, res2, res3, res4 = pikkar_method(x_min, x_max, h)
res5 = runge_kutt_1(x_min, x_max, h, y)
res6 = runge_kutt_2(x_min, x_max, h, y)
res7 = euler(x_min, x_max, h, y)

N = len(res1)

print_str()
print_header()
print_str()

x = 0
for i in range(0, N, 5000):
    print_x(x)
    print_result(res1[i])
    print_result(res2[i])
    print_result(res3[i])
    print_result(res4[i])
    print_result(res5[i])
    print_result(res6[i])
    print_result(res7[i])
    x += 5000 * h
    print("|")
print_str()

'''