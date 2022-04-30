import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable

EPS_STEP = 1e-4
EPS = 1e-4
c = 30_000_000_000
R = 0.35


Tw = 2000
T0 = 10000
k0 = 0.008
p = 4
m = 0.786
Zmin = 0
Zmax = 1

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

def divF(Z, U):
    return c * k(Z) * (Up(Z) - U)

def Kn(Z):
    return c / (3 * R * k(Z))

def Fn(Z):
    return c * k(Z) * Up(Z)

def pn(Z):
    return c * k(Z)

def Khalf(Z):
    return (Kn(Z + EPS_STEP / 2) + Kn(Z - EPS_STEP / 2)) / 2

def h2(Z):
    return Z * EPS_STEP ** 2 * R

def A(Z):
    return (Z - EPS_STEP / 2) * (Khalf(Z - EPS_STEP / 2))

def C(Z):
    return (Z + EPS_STEP / 2) * Khalf(Z + EPS_STEP / 2)

def B(Z):
    return A(Z) + C(Z) + pn(Z) * h2(Z)

def D(Z):
    return Fn(Z) * h2(Z)

def LeftBoundCond(Z0, F0, h):
    K0 = -Khalf(Z0 + h / 2) * (Z0 + h / 2) + c * R * h * h / 8 * k(Z0 + h / 2) * (Z0 + h / 2)
    M0 = -K0
    P0 = c * R * h * h / 4 * k(Z0 + h / 2) * Up(Z0 + h / 2) * (Z0 + h / 2)
    return K0, M0, P0

def RightBoundCond(Z, h):
    KN = Khalf(Z - h / 2) * (Z - h / 2) + m * c * Z * h / 2 + c * R * h * h / 8 * k(Z - h / 2) * (Z - h / 2) + R * c * h * h * k(Z) / 4
    MN = -Khalf(Z - h / 2) * (Z - h / 2) + c * R * h * h / 8 * k(Z - h / 2) * (Z - h / 2)
    PN = c * R * h * h / 4 * (k(Z - h / 2) * Up(Z - h / 2) * (Z - h / 2) + k(Z) * Up(Z))
    return KN, MN, PN

def Sweep():
    # Прямой ход
    h = EPS_STEP
    K0, M0, P0 = LeftBoundCond(0, 0, EPS_STEP)
    KN, MN, PN = RightBoundCond(1, EPS_STEP)    
    eps = [0, -K0 / M0]
    eta = [0, P0 / M0]

    x = h
    i = 1
    while x < Zmax:
        eps.append(C(x) / (B(x) - A(x) * eps[i]))
        eta.append((A(x) * eta[i] + D(x)) / (B(x) - A(x) * eps[i]))
        i += 1
        x += h
    
    # Обратный ход
    n = i
    u = [0] * (n)
    u[n - 1] = (PN - MN * eta[i]) / (KN + MN * eps[n])

    for i in range(n - 2, -1, -1):
        u[i] = eps[i + 1] * u[i + 1] + eta[i + 1]

    return u

def GetCenter(y, Z, h):
    res = [(-3 * y[0] + 4 * y[1] - y[2]) / 2 / h]
    for i in range(1, len(y) - 1):
        r = (y[i + 1] - y[i - 1]) / 2 / h
        res.append(r)
    res.append((3 * y[-1] - 4 * y[-2] + y[-3]) / 2 / h)
    return res

def F_res_deriv(u, Z):
    f = [0]
    u_res = GetCenter(u, Z, EPS_STEP)
    for i in range(1, len(Z)):
        r = -c / 3 / R / k(Z[i]) * u_res[i]
        f.append(r)
    return f

def F_res_integ(Z, un, un1, f):
    if Z == 1:
        return m * c * un / 2
    return Khalf(Z - EPS_STEP / 2) * (un - un1) / EPS_STEP

# График
name = ['U(z)', 'F(z)']
u_res = Sweep()
z_res = [i for i in np.arange(0, 1 + EPS_STEP, EPS_STEP)]

f_res = [0] * len(z_res)
up_res = [0] * len(z_res)
_divF = [0] * len(z_res)

f3_res = F_res_deriv(u_res, z_res)

for i in range(0, len(z_res) - 1):
    up_res[i] = Up(z_res[i])
    _divF[i] = divF(z_res[i], u_res[i])

for i in range(1, len(z_res)):
    f_res[i] = F_res_integ(z_res[i], u_res[i - 1], u_res[i], f_res[i - 1])

for i in range(len(z_res)):
    z_res[i] = round(z_res[i], 3)

for i in range(len(f_res)):
    f_res[i] = round(f_res[i], 10)

for i in range(len(f3_res)):
    f3_res[i] = round(f3_res[i], 10)

for i in range(len(u_res)):
    u_res[i] = round(u_res[i], 10)

for i in range(len(_divF)):
    _divF[i] = round(_divF[i], 3)

tb = PrettyTable()
tb.add_column("Z", z_res)
tb.add_column("F", f_res)
tb.add_column("F deriv", f3_res)
tb.add_column("U", u_res)
tb.add_column("divF", _divF)

with open('result.txt', 'w') as f:
    f.write(str(tb))

plt.subplot(2, 2, 1)
plt.plot(z_res, u_res, 'r', label='u')
plt.plot(z_res, up_res, 'g', label='u_p')
plt.legend()
plt.title(name[0])

plt.subplot(2, 2, 3)
plt.plot(z_res, f_res, 'g')
plt.title(name[1])

plt.subplot(2, 2, 4)
plt.plot(z_res, _divF, 'b')
plt.title("divF")

plt.subplot(2, 2, 2)
plt.plot(z_res, f3_res, 'g')
plt.title("F(z) deriv")

plt.show()

