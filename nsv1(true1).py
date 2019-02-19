import mpmath as mp
import math
import scipy.integrate as integrate
import numpy as np
import pylab
import time
from multiprocessing import Pool

np.seterr(all='ignore')

mev = 0.511
n = 1
d = 3.6
df = 0.0
tt = 1.0e8
z1 = 48
z2 = 1
a1 = 112
a2 = 1

mn = 1.835E3
gv = 2.9899e-12
aa = 1/137.03

m = mn*a1*a2/(a1+a2)

d = d/mev
df = df/mev
tt = tt*1e-10/1.16/mev


def kf(ef):
    return math.sqrt(2*m*ef)


def lf(ef):
    return (z1+1)*z2*aa*math.sqrt(m/2/ef)


def li(ei):
    return z1*z2*aa*math.sqrt(m/2/ei)


def ki(ei):
    return math.sqrt(2*m*ei)


def fun1(t, ei, ef):
    kf = np.sqrt(2*m*ef)
    lf = (z1+1)*z2*aa*np.sqrt(m/2/ef)
    ki = np.sqrt(2*m*ei)
    li = z1*z2*aa*np.sqrt(m/2/ei)
    x = 2*(1-t)/(ki/kf+kf/ki-2*t)
    a = 1 + li*1j
    b = -lf*1j
    c = 1
    res = mp.hyp2f1(a, b, c, x)
    #checkRes = (1-x)**(c-a-b) * mp.hyp2f1(c-a, c-b, c, x)
    #diff = np.sqrt((res.imag - checkRes.imag)**2 +
    #               (res.real - checkRes.real)**2)
    # if (diff > 1e-9):
    #    print("Error: " + str(diff))
    return (res.real*res.real+res.imag*res.imag)/(ki*ki+kf*kf-2*ki*kf*t)/(ki*ki+kf*kf-2*ki*kf*t)


def sig1(ei, ef):
    return integrate.quad(lambda x, args: fun1(x, args[0], args[1]), -1, 1, args=[ei, ef])[0]


def fun2(ef, ei):
    kf = np.sqrt(2*m*ef)
    lf = (z1+1)*z2*aa*np.sqrt(m/2/ef)
    ki = np.sqrt(2*m*ei)
    e = ei-d-df-1
    c = np.sqrt(e-ef)*(e-ef)*(e-ef)*(e-ef)/(ki*ki-kf*kf)/(ki*ki-kf*kf)
    sig1res = sig1(ei, ef)
    return c/(np.exp(2*np.pi*lf)-1)*sig1res


def fun3(ei, tt):
    li = z1*z2*aa*np.sqrt(m/2/ei)
    res = integrate.quad(lambda x, args: fun2(x, args), 0.0, ei-d-df-1, args=ei)[0]
    ww = ei*np.exp(-ei/tt)
    return ww*z2*z2*256*np.sqrt(2)*aa*aa*aa*aa*gv*gv*m*m*m*m*m*z1*(z1+1) * z2*z2/(105*math.pi*ei)/(1-np.exp(-2*math.pi*li))*res

def fun3_for_test(ei):
    li = z1*z2*aa*np.sqrt(m/2/ei)
    return integrate.quad(lambda x, args: fun2(x, args), 0.0, ei-d-df-1, args=ei)[0]

def nsv(tt):
    res = integrate.quad(lambda x, args: fun3(x, args), d+df+1, np.inf, args=tt)
    return math.sqrt(8/math.pi/m/tt/tt/tt)*res[0]*0.19448


def nsv_norm(tt):
    return nsv(tt)*44.722E-12

def work2(value):
    if isinstance(value, float):
        return nsv_norm(float(value)) * 44.7e-12
    else:
        raise "error"

def compare_fun3_for_test():
    x_array = np.linspace(1e8, 1e10, 6)
    x_array = x_array * 1e-10/1.16/0.511
    p = Pool()
    start_time = time.time()
    f = open('out-fuc-3.txt', 'w')
    y_array = p.map(work2, x_array)
    for i in range(0, len(y_array)):
         print(str(x_array[i]) + ' ' + str(y_array[i]), file=f)
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == "__main__":
    # Compare with original version
    # verctorized_nsv_norm = np.vectorize(nsv_norm)
    # print('Original:' + str(verctorized_nsv_norm(x_tt)))    

    compare_fun3_for_test()

    # p = Pool()
    # start_time = time.time()
    # mp_solutions = p.map(nsv_norm, x_tt)
    # f = open('out.txt', 'w')
    # print("--- %s seconds ---" % (time.time() - start_time), file=f)
    # print(str(mp_solutions))    
    # for i in range(0, len(mp_solutions)):
    #     print(str(x[i]) + ' ' + str(mp_solutions[i]), file=f)
