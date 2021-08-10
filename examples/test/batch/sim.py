# -*- encoding=utf-8 -*-
from __future__ import print_function
import sys
import time
from yade.params import table
# see https://yade-dem.org/doc/user.html#batch-queuing-and-execution-yade-batch

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
# time.sleep(5)
O.engines = [PyRunner(command='time.sleep(.005)', iterPeriod=1)]
O.run(1000, True)
print('finished')
sys.stdout.flush()
sys.exit(0)
