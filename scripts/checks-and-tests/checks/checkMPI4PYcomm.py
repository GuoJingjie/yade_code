import time
import numpy
if 'MPI' in yade.config.features:

	try:
		from yade import mpy as mp
		mp.initialize(2)
		if mp.rank > 0:
			mp.mprint("spawned!")
	except:
		raise YadeCheckError("Error in initializing mpy ")

	try:
		numTrial = 10

		if mp.rank == 0:
			listData0 = [k for k in range(10000)]
			mp.comm.barrier()
		else:
			listData1 = [k * 3 for k in range(10000)]
			mp.comm.barrier()

		t1 = time.time()
		for count in range(numTrial):
			if mp.rank == 0:
				mp.comm.send(listData0, dest=1)
				listData1 = mp.comm.recv(source=1)
			else:
				lisData0 = mp.comm.recv(source=0)
				mp.comm.send(listData1, dest=0)
		dt1 = time.time() - t1

		arrayData = numpy.arange(10000, dtype='i')
		t1 = time.time()
		for count in range(numTrial):
			if mp.rank == 0:
				mp.comm.Send([arrayData, mp.MPI.INT], dest=1, tag=77)
				status = mp.MPI.Status()
				mp.comm.Probe(source=1, tag=77, status=status)
				size = status.Get_count(mp.MPI.INT)
				data = numpy.empty(size, dtype='i')
				mp.comm.Recv([data, mp.MPI.INT], source=1, tag=77)
			elif mp.rank == 1:
				status = mp.MPI.Status()
				mp.comm.Probe(source=mp.MPI.ANY_SOURCE, tag=77, status=status)
				size = status.Get_count(mp.MPI.INT)
				data = numpy.empty(size, dtype='i')
				mp.comm.Recv([data, mp.MPI.INT], source=0, tag=77)
				mp.comm.Send([arrayData, mp.MPI.INT], dest=0, tag=77)
		dt2 = time.time() - t1

		if mp.rank == 0:
			print("____ MPI comm times: ____")
			print("1e4 integers in a python list:", dt1 / numTrial)
			print("1e4 integers in a python array:", dt2 / numTrial)
	except:
		raise YadeCheckError("Error in some MPI4PY communications")

	try:
		if mp.rank == 0:
			mp.disconnect()
	except:
		raise YadeCheckError("Error in disconnecting mpy ")
