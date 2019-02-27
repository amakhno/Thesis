import numpy as np

with open('out-fuc-3-true.txt') as f:
    content = f.readlines()

content = [x.strip() for x in content]

t_array = np.linspace(1e8, 1e10, 24)
sig_array = []

for index in range(0, len(content)):
    splitted = content[index].rsplit(' ')
    # t_array.append(float(splitted[0]))
    val = float(splitted[1])
    if (val > 0):
        sig_array.append(np.log(val))
    else:
        sig_array.append(np.log(np.abs(val)))
    


def build_a(t_arr):
    if isinstance(t_arr, np.ndarray):
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
f=open('solved.txt', 'w')
print(np.linalg.solve(a[0:7], b[0:7]), file=f)

x, _, _, _ =np.linalg.lstsq(a, b, rcond=1)
print(x, file=f)
f.close()
