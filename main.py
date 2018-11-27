from scipy import special
import numpy as np

#Function
def ReadVariable (input_file, var_name):
    line = input_file.readline()
    if var_name != line[0: len(var_name)]:
        return 0
    value = line[len(var_name) + 3:]
    return float(value)

input_file = open("input.data", "r") 
tempValue = ReadVariable(input_file, "tempValue")
item = special.hyp1f1(1, tempValue, 1)
print (item)