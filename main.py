import mpmath as mp
from scipy.integrate import dblquad
import scipy.integrate as integrate
import numpy as np
import cmath
import math
import pylab
from multiprocessing import Pool
import time


class Sig_Calculate:
    mev = 0.511
    d = 3.6
    df = 0.0
    z1 = 48
    z2 = 1
    a1 = 112
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
    k_i = 0

    def __init__(self, a1, z1, d, df, beta):
        self.a1 = a1; self.z1 = z1; self.d = d/self.mev; self.df = df/self.mev
        self.m = self.mn*a1*self.a2/(a1+self.a2)
        self.ksi_b = 6250 / 10**beta

    def Phi(self, E):
        return 1/60 * (E**2 - 1)**0.5 * (2*E**4 - 9*E**2 - 8) \
            + 1/4 * E * np.log(E + (E**2 - 1)**0.5)

    def Get_k_f(self, ef):
        return math.sqrt(2*self.m*ef)

    def Lambda_f(self, ef):
        return (self.z1+1)*self.z2*self.aplha_e*np.sqrt(self.m/2/ef)

    def Lambda_i(self, ei):
        return self.z1*self.z2*self.aplha_e*np.sqrt(self.m/2/ei)

    def Get_k_i(self, ei):
        return np.sqrt(2*self.m*ei)

    def Inner_func(self, x, eps_f):
        k_f = self.Get_k_f(eps_f)
        lambda_f = self.Lambda_f(eps_f)
        return abs(mp.hyp2f1(-self.lambda_i*1j, -lambda_f*1j, 1, x))**2 / (1 - x)**2

    def E_f(self, eps_f, eps_i):
        return eps_i - eps_f - self.d - self.df

    def Outer_func(self, eps_f, eps_i):
        e_f = self.E_f(eps_f, eps_i)
        k_f = self.Get_k_f(eps_f)
        result = self.Phi(e_f) / (np.exp(2*cmath.pi*float(self.Lambda_f(eps_f)) - 1)
                                  * k_f*(self.k_i - k_f)**4 * (self.k_i+k_f)**2)
        return result

    def X_0(self, eps_f):
        k_f = self.Get_k_f(eps_f)
        result = -4*self.k_i*k_f / (self.k_i - k_f)**2
        return result

    def ResultFunc(self, x, eps_f, args):
        eps_i = args
        result = self.Outer_func(eps_f, eps_i)*self.Inner_func(x, eps_f)
        return float(result.real)

    def sig(self, eps_i):
        self.k_i = self.Get_k_i(eps_i)
        self.lambda_i = self.Lambda_i(eps_i)
        top_limit = eps_i - self.d - self.df
        integral = dblquad(self.ResultFunc, 0, top_limit, self.X_0, lambda eps_f: 0,
                           args=[eps_i])
        mltiple = (4 * 2**0.5 / cmath.pi) * \
            ((self.gv**2 * self.aplha_e**4 * self.z1 * (self.z1 + 1) * self.z2**4 * self.m**(9/2))) \
            / (eps_i**1.5 * (1 - np.exp(-2 * cmath.pi * float(self.lambda_i)))) * self.ksi_b
        return mltiple * integral[0]

    def fun3(self, eps_i, tt):
        return math.exp(-eps_i/tt)*eps_i*self.sig(eps_i)

    def sigv(self, tt):        
        tt = tt * 1e-10/1.16/self.mev
        res = integrate.quad(self.fun3, self.d+self.df+1, np.inf, args=tt)
        return math.sqrt(8/math.pi/self.m/tt/tt/tt)*res[0] * 0.19448 * 44.7e-12

def thread_work(input_array):
    tt = input_array[0]
    a1 = input_array[1]
    z1 = input_array[2]
    d = input_array[3]
    df = input_array[4]
    beta = input_array[5]
    calc = Sig_Calculate(a1, z1, d, df, beta)
    return calc.sigv(tt)

def build_full_range(a1, z1, d, df, beta):
    t_array = np.linspace(1e8, 1e10, 24)
    input_array = []
    for i in range(0, len(t_array)):
        input_array.append([t_array[i], a1, z1, d, df, beta])
    p = Pool()
    #y_array = [thread_work(input_array[0])]
    start_time = time.time()
    f = open('out-main.txt', 'w')
    y_array = p.map(thread_work, input_array)
    for i in range(0, len(y_array)):
        print(str(t_array[i]) + ' ' + str(y_array[i]), file=f)
    print("--- %s seconds ---" % (time.time() - start_time))
    return y_array


if __name__ == "__main__":
    result = build_full_range(108, 47, 3.6, 0, 4.8)
    print(result)