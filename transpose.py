from functools import *
import fileinput
import math
data = []
for line in fileinput.input():
    data.append(list(map(lambda x: float(x), line.split())))
for i in range(len(data[0])):
    print(reduce(lambda x,y: x + "," + str(y[i]), data, "")[1:] + "," + str(math.sqrt(abs(data[1][i]))))
    



          
