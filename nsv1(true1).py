import mpmath as mp
import math
import scipy.integrate as integrate
import numpy as np

np.seterr(all='raise')

mev = 0.511
n = 1
d = 3.6
df = 0.0
tt = 1.0e8
z1 = 3.4e1
z2 = 1
a1 = 7.8e1
a2 = 1

mn = 1.835E3
gv = 2.9899e-12
aa = 1/137.03

send_ef = 0
send_ei = 0

m = mn*a1*a2/(a1+a2)
print(a1*a2/(a1+a2))

d = d/mev
df = df/mev
tt = tt*1e-10/1.16/mev


def kf(ef):
    return math.sqrt(2*m*ef)


def lf(ef):
    return (z1+1)*z2*aa*math.sqrt(m/2/ef)


def li(ei):
    return z1*z2*aa*math.sqrt(m/2/ei)


def fun3(ei):
    send_ei = ei
    return exp(-ei/tt)*ei*sigb(ei)


def ki(ei):
    return math.sqrt(2*m*ei)


def fun1(t):
    ef = send_ef
    ei = send_ei
    kf = np.sqrt(2*m*ef)
    lf = (z1+1)*z2*aa*np.sqrt(m/2/ef)
    ki = np.sqrt(2*m*ei)
    li = z1*z2*aa*np.sqrt(m/2/ei)
    x = 2*(1-t)/(ki/kf+kf/ki-2*t)
    a = 1 + li*1j
    b = -lf*1j
    c = 1
    rez = mp.hyp2f1(a, b, c, x)
    return 1e-300*(rez.real*rez.real+rez.imag*rez.imag)/(ki*ki+kf*kf-2*ki*kf*t)/(ki*ki+kf*kf-2*ki*kf*t)


def sig1(ef):
    global send_ef
    send_ef = ef
    return integrate.quad(lambda x: fun1(x), -1, 1)[0]


def fun2(ef):
    ei = send_ei
    kf = np.sqrt(2*m*ef)
    lf = (z1+1)*z2*aa*np.sqrt(m/2/ef)
    ki = np.sqrt(2*m*ei)
    e = ei-d-df-1
    c = 1e300*np.sqrt(e-ef)*(e-ef)*(e-ef)*(e-ef)/(ki*ki-kf*kf)/(ki*ki-kf*kf)
    return c/(np.exp(2*np.pi*lf)-1)*sig1(ef)


def sigb(ei):
    res = integrate.quad(lambda x: fun2(x), 0.0, ei-d-df-1)[0]
    return z2*z2*256*math.sqrt(2)*aa*aa*aa*aa*gv*gv*m*m*m*m*m*z1*(z1+1)*z2*z2/(105*math.pi*ei)/(1-math.exp(-2*math.pi*li(ei)))*res


def fun3(ei):
    global send_ei
    send_ei = ei
    li = z1*z2*aa*np.sqrt(m/2/ei)
    res = integrate.quad(lambda x: fun2(x), 0.0, ei-d-df-1)[0]
    ww = ei*np.exp(-ei/tt)
    return ww*z2*z2*256*np.sqrt(2)*aa*aa*aa*aa*gv*gv*m*m*m*m*m*z1*(z1+1)* z2*z2/(105*math.pi*ei)/(1-np.exp(-2*math.pi*li))*res


def nsv(d):
    res = integrate.quad(lambda x: fun3(x), d+df+1, 100)
    return math.sqrt(8/math.pi/m/tt/tt/tt)*res[0]*0.19448


print(nsv(d)*44.722E-12)
exit
