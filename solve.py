import numpy as np
import pylab

with open('out-main.txt') as f:
    content = f.readlines()

content = [x.strip() for x in content]

t_array = []
sig_array = []

for index in range(0, len(content)):
    splitted = content[index].rsplit(' ')
    t_array.append(float(splitted[0])/1e9)
    val = float(splitted[1]) * 1e25
    if (val > 0):
        sig_array.append(np.log(val))
    else:
        sig_array.append(np.log(np.abs(val)))
    

def build_a(t_arr):
    if isinstance(t_arr, list):
        if len(t_arr) >= 7:
            a = []
            for index in range(len(t_arr)):
                t = t_arr[index]
                a.append([1, t**(-1), t**(-1/3), t**(1/3), t, t**(5/3), np.log(t)])
            return a
        else:
            pass
    else:
        pass

a = build_a(t_array)
b = sig_array
f=open('solved.txt', 'w')

x = []; y = []

for index in range(0, 24-7):
    x.append(t_array[index])
    y.append(np.linalg.solve(a[index:index + 7], b[index:index + 7]))

def plot_solution(t_array, solution, block):
    if (block != True):
        block = False
    a1_array = []; b1_array = []
    for t_index in range(0, len(t_array)):
        current_t = t_array[t_index]
        current_a1 = 0
        for a_index in range(0, len(solution)):
            current_a1 += solution[a_index] * a[t_index][a_index]
        a1_array.append(current_a1)
        b1_array.append(b[t_index])
    pylab.figure()
    line1, = pylab.plot(t_array, a1_array, label="A", linestyle="--")
    line2, = pylab.plot(t_array, b1_array, label="B")    
    pylab.legend(handles=[line1, line2], loc=1)
    pylab.show(block=block)

print("lstsq1")
x_lstsq, _, _, _ =np.linalg.lstsq(a, b, rcond=1)
plot_solution(t_array, x_lstsq, True)
print("lstsq2")
x_lstsq, _, _, _ =np.linalg.lstsq(a[2:20], b[2:20], rcond=1)
plot_solution(t_array, x_lstsq, True)
print(x_lstsq, file=f)
f.close()

#pylab.plot(x, y)
#pylab.xticks(x, x, rotation='vertical', )
#pylab.show()

# print(np.linalg.solve(a[0:7], b[0:7]), file=f)

# f.close()
