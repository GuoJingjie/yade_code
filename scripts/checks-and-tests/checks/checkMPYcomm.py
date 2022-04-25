if 'MPI' in yade.config.features:

	try:
		from yade import mpy as mp
		mp.initialize(3)
	except:
		raise YadeCheckError("Error in initializing mpy ")

	if mp.rank == 0:
		lenOther = mp.sendCommand(executors="all", command="len(O.bodies)", wait=True)
		if lenOther != [0, 0, 0]:
			raise YadeCheckError("wrong subdomains sizes")

		## pickable objects
		## pay attention to pointer adresses, they are different! (as expected)
		## this is moving data around between independent parts of memory
		mp.sendCommand(
		        executors=1, command="O.bodies.append(Body()); O.bodies[0].shape=Sphere(radius=0.456); comm.send(O.bodies[0],dest=0)", wait=False
		)
		bodyFrom1 = mp.comm.recv(source=1)
		if not (bodyFrom1.shape.radius == 0.456):
			raise YadeCheckError("communication of pickled body failed")

		# More bodies
		# note that 'rank' is used instead of mp.rank, the scope of the mpy module is accessed here
		added = mp.sendCommand(executors="slaves", command="ids=O.bodies.append([sphere((xx,1.5+rank,0),0.5) for xx in range(-1,2)])", wait=True)

		wallId = O.bodies.append(box(center=(0, 0, 0), extents=(2, 0, 1), fixed=True))  #master only

		added = mp.sendCommand(executors="all", command="len(O.bodies)", wait=True)
		if added != [1, 4, 3]:
			raise YadeCheckError("wrong subdomains sizes")

		# set remote attributes
		mp.sendCommand(executors="slaves", command="list(map(lambda b: setattr(b,'subdomain',rank),O.bodies))", wait=True)
		assigned = mp.sendCommand("slaves", "len([b for b in O.bodies if b.subdomain==rank])", True)

		if assigned != [4, 3]:
			raise YadeCheckError("failed setting remote attribute")

		# Kill workers pool, then start again with a new
		# any data in the old threads is lost
		mp.disconnect()
		mp.initialize(4)
		# not all workers will be displayed because mp.MAX_RANK_OUTPUT=5 by default
		newSizes = mp.sendCommand(executors="slaves", command="len(O.bodies)", wait=True)
		if newSizes != [0, 0, 0]:
			raise YadeCheckError("dirty scenes after a reconnexion")

	elif mp.rank == None:
		raise YadeCheckError("uninitialized mpi rank")

	if mp.rank == 0:
		try:
			mp.mprint("MPYcomm disconnect")
			mp.disconnect()
		except:
			raise YadeCheckError("Error in disconnecting mpy ")
