import matplotlib.pyplot as plot

import fileinput

data = []
for line in fileinput.input():
    data.append(list(map(lambda x: float(x), line.split())))
for i in range(len(data)//2):
    
    plot.plot(data[2*i], data[2*i+1])
plot.show()
