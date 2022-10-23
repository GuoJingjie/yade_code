#author: bruno.chareyre@grenoble-inp.fr

# demonstrates using standalone bridge solver and some timing

# output ont intel i7:
    # triangulate data: 0.42063212394714355 sec
    # 100000 solutions with smart locate: 0.10309433937072754 sec
    # 100000 solutions with no cache locate: 0.2700941562652588 sec



import time
 
l=CapillarityEngine()
l.liquidTension = 1
Nsolve = 100000

import numpy as np

t1 = time.time()
i = l.solveStandalone(1,1,1,0.1) # dummy solve to trigger triangulation
dt = (time.time() - t1)
print("triangulate data:",dt,"sec")

t1 = time.time()
for pp in np.linspace(0,10,Nsolve):
    i = l.solveStandalone(1,1,pp,0.1,i)
dt = (time.time() - t1)
print(Nsolve,"solutions with smart locate:",dt,"sec")

t1 = time.time()
for pp in np.linspace(0,10,Nsolve):
    i = l.solveStandalone(1,1,pp,0.1)
dt = (time.time() - t1)
print(Nsolve,"solutions with no cache locate:",dt,"sec")

## Plotting 
v = []
p = []
f = []
for pp in np.linspace(0,10,Nsolve):
    i = l.solveStandalone(1,1,pp,0.1,i)
    v.append(i.vMeniscus)
    p.append(pp)
    f.append(i.fCap[0])
dt = (time.time() - t1)

from matplotlib import pyplot as plt
plt.plot(p,v)
plt.show(False)
