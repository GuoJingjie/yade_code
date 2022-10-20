#author: bruno.chareyre@grenoble-inp.fr

# demonstrates using standalone bridge solver and some timing

import time
 
l=Law2_ScGeom_CapillaryPhys_Capillarity1()
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
