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
        return abs(mp.hyp2f1(-self.lambda_i*1j, -lambda_f*1j, 1, x)**2) / (1 - x)**2

    def E_f(self, eps_f, eps_i):
        return eps_i - eps_f - self.d - self.df

    def Outer_func(self, eps_f, eps_i):
        e_f = self.E_f(eps_f, eps_i)
        k_f = self.Get_k_f(eps_f)
        result = self.Phi(e_f) / (np.exp(2*cmath.pi*float(self.Lambda_f(eps_f)))
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
        res = integrate.quad(self.fun3, self.d+self.df+1, np.inf, args=tt)
        return math.sqrt(8/math.pi/self.m/tt/tt/tt)*res[0] #*0.09747


def work(value):
    if isinstance(value, float):
        calc = Sig_Calculate()
        return calc.sig(float(value))
    else:
        raise "error"

def work2(value):
    if isinstance(value, float):
        calc = Sig_Calculate()
        return calc.sigv(float(value)) * 44.7e-12
    else:
        raise "error"


def compare_fun3_for_test():
    x_array = np.linspace(1e8, 1e10, 6)
    x_array = x_array * 1e-10/1.16/0.511
    p = Pool()
    start_time = time.time()
    f = open('out-fuc-3-true.txt', 'w')
    y_array = p.map(work2, x_array)
    for i in range(0, len(y_array)):
        print(str(x_array[i]) + ' ' + str(y_array[i]), file=f)
    print("--- %s seconds ---" % (time.time() - start_time))

def compare_files():
    x_array = np.linspace(9.04500978473581, 20, 100)
    y_array_true = []
    y_array = []
    with open('out-fuc-3-true.txt') as inf:
        for line in inf:
            parts = line.split()  # split line into parts
            if len(parts) > 1:   # if at least 2 parts/columns
                y_array_true.append(float(parts[1]))   # print column 2
    with open('out-fuc-3.txt') as inf:
        for line in inf:
            parts = line.split()  # split line into parts
            if len(parts) > 1:   # if at least 2 parts/columns
                y_array.append(float(parts[1]))   # print column 2
    pylab.plot(x_array, y_array, x_array, y_array_true)
    pylab.yscale('log')
    pylab.show()


if __name__ == "__main__":
    compare_fun3_for_test()
    #compare_files()
    # x_array = np.linspace(9.04500978473581, 20, 100)
    # p = Pool()
    # start_time = time.time()
    # y_array = p.map(work, x_array)
    # print("--- %s seconds ---" % (time.time() - start_time))
    # pylab.plot(x_array, y_array)
    # pylab.yscale('log')
    # pylab.show()
