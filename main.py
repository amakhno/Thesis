from scipy import special
from scipy.integrate import dblquad
import numpy as np
import math

#Function
def ReadVariable (input_file, var_name):
    line = input_file.readline()
    if var_name != line[0: len(var_name)]:
        return 0
    value = line[len(var_name) + 3:]
    return float(value)

input_file = open("input.data", "r") 
temp_value = ReadVariable(input_file, "tempValue")
b = ReadVariable(input_file, "b")
c = ReadVariable(input_file, "c")
input_file.close()
g_v = 0.1
aplha_e = 0.1
Z = 1
Z1 = 3
mu = 0.1
plank = 1.054572e-27
ksi_b = 1
e = 0.5
lambda_i = 0.1
lambda_f = 0.3
k_i = 0.2
eps_i = (plank**2 * k_i**2)/(2*mu)
multiple = (4*2**0.5/math.pi) * ( (g_v**2 * aplha_e**4 * Z * (Z + 1) * Z1**4 * mu**(9/2)) ) \
    / ( eps_i**1.5 * (1 - math.exp(-2*math.pi * lambda_i))) * ksi_b

def Lambda_i (k_i) :
    return Z*Z1*e**2*mu / (plank**2*k_i)

def Lambda_f (k_f) :
    return (Z+1)*Z1*e**2*mu/ (plank**2*k_f)

def Get_k_i (eps_i) :
    return math.sqrt(2*mu*eps_i) / plank**2

def Get_k_f (eps_f) :
    return math.sqrt(2*mu*eps_f) / plank**2

def Phi (value):
    return 0.2

def ten (x, y) : 
    return x*y

inner_func = lambda x: abs(special.hyp1f2(-lambda_i*1j, -lambda_f*1j, 1, x)**2)/(1 - x)**2
"outer_func = lambda eps_f:"
integral = dblquad(ten, 0, 0.5, lambda x: 0, lambda x: 1-2*x)
print (integral)