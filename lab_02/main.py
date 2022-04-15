import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from prettytable import PrettyTable

EPS = 1e-4
c = 30_000_000_000
R = 0.35


Tw = 2000
T0 = 10000
k0 = 0.08
p = 4
m = 0.786

def T(z):
    return (Tw - T0) * (z ** p) + T0

def k(z):
    t = T(z)
    return k0 * ((t/300)**2)

def Up(z):
    return (3.084 * 10**(-4)) / (np.exp(4.799 * 10000 / T(z)) - 1)

# Поправка F / Z - это из-за свойства производной
def F_z(Z, F, u):
    if Z == 0:
        return ((c * R)/2) * k(Z) * (Up(Z) - u)
    else:
        return c * R * k(Z) * (Up(Z) - u) - (F / Z)

def U_z(Z, F):
    return -(3 * R * F * k(Z)) / c

# Метод Рунге-Кутта 4 степени точности
def runge_kutt_4(z0, f0, u0, h, z_max):
    z_res = []
    u_res = []
    f_res = []

    z = z0
    u = u0
    f = f0

    z_res.append(z)
    u_res.append(u)
    f_res.append(f)

    while z < z_max:
        q_1 = h * U_z(z, f)
        k_1 = h * F_z(z, f, u)

        q_2 = h * U_z(z + h / 2, f + k_1 / 2)
        k_2 = h * F_z(z + h / 2, f + k_1 / 2, u + q_1 / 2)

        q_3 = h * U_z(z + h / 2, f + k_2 / 2)
        k_3 = h * F_z(z + h / 2, f + k_2 / 2, u + q_2 / 2)

        q_4 = h * U_z(z + h, f + k_3)
        k_4 = h * F_z(z + h, f + k_3, u + q_3)

        u += ((q_1 + 2 * q_2 + 2 * q_3 + q_4) / 6)
        f += ((k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6)
        z += h

        z_res.append(z)
        u_res.append(u)
        f_res.append(f)
    return z_res, u_res, f_res


# Подбор методом стрельбы начального условия для функции u(z)
# Для начала ищем интересующий интервал функции Psi, а затем итерационно будем искать решения Psi(Xi) == 0

# Объявляем функцию Psi, которую используем для поиска решения для кси
def Psi(F, m, c, u):
    return F - m * c * u / 2

# В этой части мы ищем, есть ли в принципе решение у уравнения
def findPossibleSolutions(xi, xi_step, h, z0, f0):
    hasSolution = False
    _, u_res, f_res = runge_kutt_4(z0, f0, xi * Up(z0), h, 1)
    psi_last = Psi(f_res[-1], m, c, u_res[-1])

    while xi <= 1:
        _, u_res, f_res = runge_kutt_4(z0, f0, xi * Up(z0), h, 1)

        tmpPsi = Psi(f_res[-1], m, c, u_res[-1])

        if (psi_last > 0 and tmpPsi < 0) or (psi_last < 0 and tmpPsi > 0):
            hasSolution = True
            break
        
        psi_last = tmpPsi
        xi += xi_step
    
    return hasSolution, xi

# Вторая часть итерационного метода
def FindXtrue(xi, xi_step, z0, f0, h):
    xi_1 = xi - xi_step
    xi_2 = xi
    xi_true = (xi_1 + xi_2)/2

    while abs(xi_1 - xi_2)/xi_true > EPS:
        _, u_res_1, f_res_1 = runge_kutt_4(z0, f0, xi_1 * Up(z0), h, 1)
        _, u_res_t, f_res_t = runge_kutt_4(z0, f0, xi_true * Up(z0), h, 1)
        _, u_res_2, f_res_2 = runge_kutt_4(z0, f0, xi_2 * Up(z0), h, 1)

        psi_1 = Psi(f_res_1[-1], m, c, u_res_1[-1])
        psi_t = Psi(f_res_t[-1], m, c, u_res_t[-1])
        psi_2 = Psi(f_res_2[-1], m, c, u_res_2[-1])

        if psi_1 > 0 and psi_t < 0:
            xi_2 = xi_true
        elif psi_t > 0 and psi_2 < 0:
            xi_1 = xi_true
        elif psi_1 < 0 and psi_t > 0:
            xi_2 = xi_true
        elif psi_t < 0 and psi_2 > 0:
            xi_1 = xi_true

        xi_true = (xi_1 + xi_2)/2
    
    return xi_true

def Solve(xi, xi_step, h, z0, f0):

    hasSolution, xi = findPossibleSolutions(xi, xi_step, h, z0, f0)

    # Если решения нет, то покидаем программу
    if hasSolution == False:
        return None, None, None, None
    else:
        xi_true = FindXtrue(xi, xi_step, z0, f0, h)
        # Запускаем результат на подготовленных данных
        z_res, u_res, f_res = runge_kutt_4(z0, f0, xi_true * Up(z0), h, 1)
        return z_res, u_res, f_res, xi_true

xi = 0.01
xi_step = 0.001
h = 0.01
z0 = 0
f0 = 0

z_res, u_res, f_res, xi_true = Solve(xi, xi_step, h, z0, f0)
if (z_res is None):
    exit()
u_p_res = []
for z in z_res:
    u_p_res.append(Up(z))

algos = []

cols = ['z', 'u_p', 'u', 'F']
table = PrettyTable(cols)
for z, u_p_r, u_r, f_r in zip(z_res, u_p_res, u_res, f_res):
    table.add_row([round(z, 14), round(u_p_r, 14), round(u_r, 14), round(f_r, 14)])

print('h: ' + str(h))
print('xi step: ' + str(xi_step))
print('xi: ' + str(xi_true))
print(table)

with open("/home/prianechka/Education/BMSTU/Modeling/lab_02/result.txt", 'w') as f:
    f.write('h: ' + str(h) + '\n')
    f.write('xi step: ' + str(xi_step) + '\n')
    f.write('xi: ' + str(xi_true) + '\n')
    f.write(str(table))
    f.close()

fig, axs = plt.subplots(nrows=1, ncols=3)
sns.lineplot(x=z_res, y=u_p_res, hue=['Up(z)'] * len(u_p_res), ax=axs[0]).set(\
         title = "Зависимость Up(Z) от Z", ylabel = "Up(Z)")
sns.lineplot(x=z_res, y=u_res, hue=['u(z)'] * len(u_res), ax=axs[1]).set(\
         title = "Зависимость U(z) от Z", xlabel = "Z", ylabel = "U(Z)")
sns.lineplot(x=z_res, y=f_res, hue=['F(z)'] * len(f_res), ax=axs[2]).set(\
         title = "Зависимость F(Z) от Z", xlabel = "Z", ylabel = "F(Z)")
plt.show()


'''
def FindRelativesT0(xi, xi_step, h, z0, f0):
    arrayT = []
    result = []
    global T0
    tmpT = T0
    T0 = 2000
    Tmax = 50000
    while T0 < Tmax:
        arrayT.append(T0)
        _, u_res, _, _ = Solve(xi, xi_step, h, z0, f0)
        if (u_res is None):
            result.append(None)
        else:
            result.append(u_res[0])
        T0 += 1000
    T0 = tmpT

    return arrayT, result

def FindRelativesTw(xi, xi_step, h, z0, f0):
    arrayT = []
    result = []
    global Tw
    tmpT = Tw
    Tw = 200
    Tmax = 5000
    while Tw < Tmax:
        arrayT.append(Tw)
        _, u_res, _, _ = Solve(xi, xi_step, h, z0, f0)
        if (u_res is None):
            result.append(None)
        else:
            result.append(u_res[0])
        Tw += 100
    Tw = tmpT

    return arrayT, result

def FindRelativesP(xi, xi_step, h, z0, f0):
    arrayT = []
    result = []
    global p
    tmpT = p
    p = 1
    Pmax = 15
    while p < Pmax:
        arrayT.append(p)
        _, u_res, _, _ = Solve(xi, xi_step, h, z0, f0)
        if (u_res is None):
            result.append(None)
        else:
            result.append(u_res[0])
        p += 1
    p = tmpT

    return arrayT, result

def FindRelativesK0(xi, xi_step, h, z0, f0):
    arrayT = []
    result = []
    global k0
    tmpT = k0
    k0 = 0.0001
    kmax = 0.002
    while k0 < kmax:
        arrayT.append(k0)
        _, u_res, _, _ = Solve(xi, xi_step, h, z0, f0)
        if (u_res is None):
            result.append(None)
        else:
            result.append(u_res[0])
        k0 += 0.0001
    k0 = tmpT

    return arrayT, result

def FindRelativesR(xi, xi_step, h, z0, f0):
    arrayT = []
    result = []
    global R
    tmpT = R
    R = 0.01
    Rmax = 1
    while R < Rmax:
        arrayT.append(R)
        _, u_res, _, _ = Solve(xi, xi_step, h, z0, f0)
        if (u_res is None):
            result.append(None)
        else:
            result.append(u_res[0])
        R += 0.01
    R = tmpT

    return arrayT, result

#arrayT, result = FindRelativesT0(xi, xi_step, h, z0, f0) 
#cols = ['T0', 'U(0)']
#table = PrettyTable(cols)
#for t, res in zip(arrayT, result):
#    table.add_row([t, res])
#print(table)


#plt.title("Зависимость U(0) от T(0)")
#plt.plot(arrayT, result)
#plt.xlabel("T0")
#plt.ylabel("U(0)")
#plt.show()

#arrayTw, resultW = FindRelativesTw(xi, xi_step, h, z0, f0) 
#plt.title("Зависимость U(0) от Tw")
#plt.plot(arrayTw, resultW)
#plt.xlabel("Tw")
#plt.ylabel("U(0)")
#plt.show()

#arrayP, resultP = FindRelativesP(xi, xi_step, h, z0, f0) 
#plt.title("Зависимость U(0) от P")
#plt.plot(arrayP, resultP)
#plt.xlabel("P")
#plt.ylabel("U(0)")
#plt.show()

#arrayP, resultP = FindRelativesK0(xi, xi_step, h, z0, f0) 
#plt.title("Зависимость U(0) от k0")
#plt.plot(arrayP, resultP)
#plt.xlabel("k0")
#plt.ylabel("U(0)")
#plt.show()

#arrayP, resultP = FindRelativesR(xi, xi_step, h, z0, f0) 
#plt.title("Зависимость U(0) от R")
#plt.plot(arrayP, resultP)
#plt.xlabel("R")
#plt.ylabel("U(0)")
#plt.show()
'''