import matplotlib.pyplot as plot
import numpy as np
import math

N = int(input().split()[0])
print("Plotting density function for {0} particles...".format(N))

class Particle():
    def __init__(self):
        self.mass = 0
        self.pos = [0,0,0]
        self.vel = [0,0,0]
        self.softening = 0
        self.potential = 0

particles = [Particle() for i in range(N)]

xcoords = []
ycoords = []
zcoords = []
for p in particles:
    p.mass = float(input())
for p in particles:
    p.pos[0] = float(input())
    xcoords.append(p.pos[0])
for p in particles:
    p.pos[1] = float(input())
    ycoords.append(p.pos[1])
    
for p in particles:
    p.pos[2] = float(input())
    zcoords.append(p.pos[2])
for p in particles:
    p.vel[0] = float(input())
for p in particles:
    p.vel[1] = float(input())
for p in particles:
    p.vel[2] = float(input())
for p in particles:
    p.softening = float(input())
for p in particles:
    p.potential = float(input())

distances = []
totalMass = 0
bucketSize = 1.0
a = 10
ind = np.arange(1000)
buckets= [0 for i in range(1000)]
for i, p in enumerate(particles):
    d =(np.linalg.norm(p.pos))
    buckets[int(d/bucketSize)]+=1
    totalMass += p.mass

errors = [math.sqrt(buckets[i]) for i in range(len(buckets))]
analytical= [totalMass / 2 / 3.14159265 * a / (r+0.5)/bucketSize / (bucketSize * (r+0.5) + a)**3 for r in range(len(buckets))]
#print("Particle {0} has a maximum dist of {1} at pos {2}".format(maxP, maxDist, particles[maxP].pos))

plot.gca().set(title='Distance distribution', ylabel='rho')
plot.bar(ind, buckets)
plot.bar(ind, errors)
plot.plot(ind, analytical)
plot.plot(ind, analytical)
plot.figure()
plot.scatter(xcoords,ycoords)
plot.scatter(xcoords,zcoords)

plot.show()


