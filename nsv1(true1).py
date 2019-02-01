import mpmath as mp
import math
import scipy.integrate as integrate
import numpy as np
import pylab
import multiprocessing

np.seterr(all='print')

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
    res = mp.hyp2f1(a, b, c, x)
    checkRes = (1-x)**(c-a-b) * mp.hyp2f1(c-a, c-b, c, x)
    diff = np.sqrt((res.imag - checkRes.imag)**2 +
                   (res.real - checkRes.real)**2)
    if (diff > 1e-9):
        print("Error: " + str(diff))
    return (res.real*res.real+res.imag*res.imag)/(ki*ki+kf*kf-2*ki*kf*t)/(ki*ki+kf*kf-2*ki*kf*t)


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
    c = np.sqrt(e-ef)*(e-ef)*(e-ef)*(e-ef)/(ki*ki-kf*kf)/(ki*ki-kf*kf)
    sig1res = sig1(ef)
    return c/(np.exp(2*np.pi*lf)-1)*sig1res


def fun3(ei):
    global send_ei
    send_ei = ei
    li = z1*z2*aa*np.sqrt(m/2/ei)
    res = integrate.quad(lambda x: fun2(x), 0.0, ei-d-df-1)[0]
    ww = ei*np.exp(-ei/tt)
    return ww*z2*z2*256*np.sqrt(2)*aa*aa*aa*aa*gv*gv*m*m*m*m*m*z1*(z1+1) * z2*z2/(105*math.pi*ei)/(1-np.exp(-2*math.pi*li))*res


def nsv(d):
    res = integrate.quad(lambda x: fun3(x), d+df+1, np.inf)
    return math.sqrt(8/math.pi/m/tt/tt/tt)*res[0]*0.19448

def nsv_norm(d):
    return nsv(d)*44.722E-12

t_array = np.linspace(1e8, 1e10, 24)
#vectorisedFunc3 = np.vectorize(nsv_norm)
#t_result = np.linspace(1e8, 1e10, 24)
tt = 1e8
print(nsv_norm(d))
f1 = open('./nsv1(true1).out', 'w')
#for i in range(0, len(t_array)):
    #tt = t_array[i]
    #res = nsv_norm(d)
    #print(str(t_array[i]) + ' ' + str(res), file=f1)

# print(nsv(d)*44.722E-12)
# if __name__ == "__main__":
#     vectorisedFunc3 = np.vectorize(fun3)

#     xz = np.linspace(d+df+1+0.2, 10, 5)

#     yz = vectorisedFunc3(xz)
#     #p = multiprocessing.Pool(11)
#     #mp_solutions = p.map(fun3, x)

#     pylab.autoscale(enable=True, axis='both', tight=True)
#     pylab.plot(xz,yz)
#     pylab.show()
