import numpy as np

with open('out.txt') as f:
    content = f.readlines()

content = [x.strip() for x in content]

t_array = []
sig_array = []

for index in range(1, len(content)):
    splitted = content[index].rsplit(' ')
    t_array.append(float(splitted[0]))
    sig_array.append(float(splitted[1]))


def build_a(t_arr):
    if isinstance(t_arr, list):
        if len(t_arr) >= 7:
            a = []
            for index in range(len(t_arr)):
                t = t_arr[index]
                a.append([1, t, t**2, t**3, t**4, t**5, t**6])
            return a
        else:
            pass
    else:
        pass
a = build_a(t_array)
b = sig_array
f=open('solve.txt', 'w')
print(np.linalg.solve(a[0:7], b[0:7]), file=f)

x, _, _, _ =np.linalg.lstsq(a, b)
print(x, file=f)