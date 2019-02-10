import mpmath as mp
from scipy.integrate import dblquad
import numpy as np
import cmath
import math
import pylab

mev = 0.511
d = 3.6
df = 0.0
z1 = 3.4e1
z2 = 1
a1 = 7.8e1
a2 = 1

mn = 1.835E3
gv = 2.9899e-12
aplha_e = 1/137.03

# substitute
ksi_b = 1
plank = 1
e = 1

m = mn*a1*a2/(a1+a2)

d = d/mev
df = df/mev


# global
lambda_i = 0

def Phi(E):
    return 1/60 * (E**2 - 1)**0.5 * (2*E**4 - 9*E**2 - 8) \
        + 1/4 * E * np.log(E + (E**2 - 1)**0.5)


def Get_k_f(ef):
    return math.sqrt(2*m*ef)


def Lambda_f(ef):
    return (z1+1)*z2*aplha_e*np.sqrt(m/2/ef)


def Lambda_i(ei):
    return z1*z2*aplha_e*np.sqrt(m/2/ei)


def Get_k_i(ei):
    return np.sqrt(2*m*ei)


def Inner_func(x, eps_f):
    k_f = Get_k_f(eps_f)
    lambda_f = Lambda_f(eps_f)
    return abs(mp.hyp2f1(-lambda_i*1j, -lambda_f*1j, 1, x)**2) / (1 - x)**2


def E_f(eps_f, eps_i):
    return eps_i - eps_f - d - df


def Outer_func(eps_f, eps_i):
    e_f = E_f(eps_f, eps_i)
    k_f = Get_k_f(eps_f)
    result = Phi(e_f) / (np.exp(2*cmath.pi*float(Lambda_f(eps_f)))
                         * k_f*(k_i - k_f)**4 * (k_i+k_f)**2)
    return result


def X_0(eps_f):
    k_f = Get_k_f(eps_f)
    result = -4*k_i*k_f / (k_i - k_f)**2
    return result


def ResultFunc(x, eps_f, args):
    eps_i = args
    result = Outer_func(eps_f, eps_i)*Inner_func(x, eps_f)
    return float(result.real)


def sig(eps_i):
    global lambda_i
    global k_i
    k_i = Get_k_i(eps_i)
    lambda_i = Lambda_i(eps_i)
    top_limit = eps_i - d - df
    integral = dblquad(ResultFunc, 0, top_limit, X_0, lambda eps_f: 0,
                       args=[eps_i])
    mltiple = (4 * 2**0.5 / cmath.pi) * \
        ((gv**2 * aplha_e**4 * z1 * (z1 + 1) * z2**4 * m**(9/2))) \
        / (eps_i**1.5 * (1 - np.exp(-2 * cmath.pi * float(lambda_i)))) * ksi_b
    return mltiple * integral[0]

x_array = np.linspace(9.04500978473581, 20, 2)
y_array = []
for i in range(0, len(x_array)):
        x = float(x_array[i])
        y_array.append(sig(x))
        pass

pylab.plot(x_array, y_array)
pylab.yscale('log')
pylab.show ()


