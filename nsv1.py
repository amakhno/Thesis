import mpmath as mp
import math
import scipy.integrate as integrate

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
send_ef = 0

m = mn*a1*a2/(a1+a2)
print(a1*a2/(a1+a2))


def kf(ef):
    return math.sqrt(2*m*ef)


def lf(ef):
    return (z1+1)*z2*aa*math.sqrt(m/2/ef)


def li(ei):
    return z1*z2*aa*math.sqrt(m/2/ei)


def sigv(d):
    integral(d+df+1, 100, 3, igt)
    return math.sqrt(8/pi/m/tt/tt/tt)*igt*0.09747


def fun3(ei):
    send_ei = ei
    return exp(-ei/tt)*ei*sigb(ei)


def ki(ei):
    return math.sqrt(2*m*ei)


def fun1(t):
    a = 1 + 2j
    b = 2 + 3j
    c = 3 + 2j
    x = 0.2
    rez = mp.hyp2f1(a, b, c, x)
    print(rez)
    ef = send_ef
    ei = send_ei
    kki = math.sqrt(2*m*ei)
    kkf = math.sqrt(2*m*ef)
    print(kki, kkf, t)
    x = 2*(1-t)/(ki(ei)/kf(ef)+kf(ef)/ki(ei)-2*t)
    llf = (z1+1)*z2*aa*math.sqrt(m/2/ef)
    lli = z1*z2*aa*math.sqrt(m/2/ei)
    a = 1 + 2j
    b = 3j
    c = 3 + 2j
    x = 0.2
    print(rez)
    per = kki*kki+kkf*kkf
    return 1e-300*(rez.real*rez.real+rez.imag*rez.imag)/(per-2*kki*kkf*t)/(per-2*kki*kkf*t)


def sig1(ef):
    send_ef = ef
    return integrate.quad(lambda x: fun1(x), -1, 1)[0]


def fun2(ef):
    ei = send_ei
    e = ei-d-df-1
    c = 1e300 * math.sqrt(e-ef)*(e-ef)*(e-ef)*(e-ef)/(ki(ei)
                                                      * ki(ei)-kf(ef)*kf(ef))/(ki(ei)*ki(ei)-kf(ef)*kf(ef))
    return c/(math.exp(2*math.pi*lf(ef))-1)*sig1(ef)


def sigb(ei):
    res = integrate.quad(lambda x: fun2(x), 0.0, ei-d-df-1)[0]
    return z2*z2*256*math.sqrt(2)*aa*aa*aa*aa*gv*gv*m*m*m*m*m*z1*(z1+1)*z2*z2/(105*math.pi*ei)/(1-exp(-2*pi*li(ei)))*res


def fun3(ei):
    send_ei = ei
    return math.exp(-ei/tt)*ei*sigb(ei)


def sigv(d):
    res = integrate.quad(lambda x: fun3(x), d+df+1, 100)
    return math.sqrt(8/pi/m/tt/tt/tt)*igt*0.09747


print(sigv(d)*44.722E-12)
exit
