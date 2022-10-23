#author: bruno.chareyre@grenoble-inp.fr
>> import os
>>> os.system("wget -O capillaryfile.txt https://gitlab.com/yade-dev/yade-data/-/raw/main/capillaryFiles/capillaryfile.txt?inline=false")
# test standalone bridge interpolator
l=CapillarityEngine()
l.liquidTension = 1

i = l.solveStandalone(1,1,1,0.1) # dummy solve to trigger triangulation
if abs(i.vMeniscus - 0.22317611410950786) > 1e-10: raise YadeCheckError("capillary bridge volume changed")
if abs(i.fCap[0] - 3.97627790531559322) > 1e-10: raise YadeCheckError("capillary force changed")
