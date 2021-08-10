# -*- encoding=utf-8 -*-
from __future__ import print_function
import sys
import re
import os
import time
from yade.params import table
# see https://yade-dem.org/doc/user.html#batch-execution-on-job-based-clusters-oar

readParamsFromTable(unknownOk=True,
                    important=6,
                    unimportant='foo',
                    this=-1,
                    notInTable='notInTable'
                    )
print(O.tags['description'])
print('important', table.important)
print('unimportant', table.unimportant)
print(O.tags['params'].replace(',', '_'))
print(O.tags['defaultParams'])

# Stop script before walltime ends
m_hasWalltime = 'YADE_WALLTIME' in os.environ
m_Walltime = 0
m_WalltimeStopAt = 0.95

if m_hasWalltime:
    walltime = os.environ['YADE_WALLTIME']
    if re.search("^([0-9]+):([0-9]{2}):([0-9]{2})$", walltime):
        w = re.match("^([0-9]+):([0-9]{2}):([0-9]{2})$", walltime)
        g = w.groups()
        m_Walltime = 3600 * float(g[0]) + 60 * float(g[1]) + float(g[2])
        print("Will run for %i seconds" % m_Walltime)
    else:
        print("Wrong walltime format.")
        m_hasWalltime = False


def SaveAndQuit():
    print("Quit due to walltime expiration")
    O.stopAtIter = O.iter + 1


# time.sleep(5)
O.engines = [PyRunner(command='time.sleep(.005)', iterPeriod=1)]

if m_hasWalltime:
    O.engines += [PyRunner(command='SaveAndQuit()', realPeriod=m_Walltime * m_WalltimeStopAt)]

O.run(1000, True)
print('finished')
sys.stdout.flush()
sys.exit(0)
