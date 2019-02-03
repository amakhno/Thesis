with open('out.txt') as f:
    content = f.readlines()

content = [x.strip() for x in content] 

t_array = []
sig_array = []

for index in range(0, len(content)):
    splitted = content[index].rsplit(' ')
    t_array.append(splitted[0])
    sig_array.append(splitted[1])

