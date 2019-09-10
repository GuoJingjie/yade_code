'''
This module defines mpirun(), a parallel implementation of run() using a distributed memory approach. Message passing is done with mpi4py mainly, however some messages are also handled in c++ (with openmpi).


The distribution logic is as follows:
1. Instanciate a complete, ordinary, yade scene
2. Insert subdomains as special yade bodies. This is somehow similar to adding a clump body on the top of clump members
3. Broadcast this scene to all workers. In the initialization phase the workers will:
	- define the bounding box of their assigned bodies and return it to other workers
	- detect which assigned bodies are virtually in interaction with other domains (based on their bounding boxes) and communicate the lists to the relevant workers
	- erase the bodies which are neither assigned nor virtually interacting with the subdomain
4. Run a number of 'regular' iterations without re-running collision detection (verlet dist mechanism)
5. In each regular iteration the workers will:
	- calculate internal and cross-domains interactions
	- execute Newton on assigned bodies (modified Newton skips other domains)
	- send updated positions to other workers and partial force on floor to master

Yet to be implemented is the global update of domain bounds and new collision detection. It could be simply achieved by importing all bodies back in master process and re-running an initialization/distribution but there are certainly mmore efficient techniques to find.

#RULES:
	#- intersections[0] has 0-bodies (to which we need to send force)
	#- intersections[thisDomain] has ids of the other domains overlapping the current ones
	#- intersections[otherDomain] has ids of bodies in _current_ domain which are overlapping with other domain (for which we need to send updated pos/vel)

#HINTS:
- handle subD.intersections with care (same for mirrorIntersections). subD.intersections.append() will not reach the c++ object. subD.intersections can only be assigned (a list of list of int)

'''

import sys,os,inspect
import time
from mpi4py import MPI
import numpy as np
import yade.bisectionDecomposition as dd

sys.stderr.write=sys.stdout.write #so we see error messages from workers

this = sys.modules[__name__]

commSplit = False
worldComm = MPI.COMM_WORLD
color = 3; key =0; 
comm = worldComm.Split(color, key)
rank = comm.Get_rank() 
numThreads = comm.Get_size() 


waitingCommands=False #are workers currently interactive?
userScriptInCheckList=""	# the simulation script from which mpy.py is used

ACCUMULATE_FORCES=True #control force summation on master's body. FIXME: if false master goes out of sync since nothing is blocking rank=0 thread
VERBOSE_OUTPUT=False
NO_OUTPUT=False
MAX_RANK_OUTPUT=5 #larger ranks will be skipped in mprint
SEND_SHAPES=False #if false only bodies' states are communicated between threads, else shapes as well (to be implemented)
ERASE_REMOTE = True # True is MANDATORY. Erase bodies not interacting wit a given subdomain? else keep dead clones of all bodies in each scene
ERASE_REMOTE_MASTER = True # erase remotes on master or keep them for fast merge (updating just b.state)
OPTIMIZE_COM=True
USE_CPP_MPI=True and OPTIMIZE_COM
YADE_TIMING=True #report timing.stats()?
MERGE_SPLIT = False
MERGE_W_INTERACTIONS = True
COPY_MIRROR_BODIES_WHEN_COLLIDE = True  # True is MANDATORY 
RESET_SUBDOMAINS_WHEN_COLLIDE = False
DOMAIN_DECOMPOSITION = False
NUM_MERGES = 0
SEND_BYTEARRAYS = True
ENABLE_PFACETS = False    #PFacets need special (and expensive) tricks, if PFacets are not used skip the tricks
DISTRIBUTED_INSERT = False  #if True each worker is supposed to "O.bodies.insertAtId" its own bodies 
fibreList = []
FLUID_COUPLING = False
fluidBodies = [] 

#tags for mpi messages
_SCENE_=11
_SUBDOMAINSIZE_=12
_INTERSECTION_=13
_ID_STATE_SHAPE_=14
_FORCES_=15
_MIRROR_INTERSECTIONS_ = 16
_POS_VEL_ = 17
_BOUNDS_ = 18
_MASTER_COMMAND_ = 19
_RETURN_VALUE_ = 20

#for coloring processes outputs differently
bcolors=['\033[95m','\033[94m','\033[93m','\033[92m','\033[91m','\033[90m','\033[95m','\033[93m','\033[91m','\033[1m','\033[4m','\033[0m']

#def mprint(m): #this one will print regardless of VERBOSE_OUTPUT
	#if 1:
		#if rank>0:
			#print (bcolors.WARNING if rank==1 else bcolors.OKBLUE) +"worker"+str(rank)+": "+m+bcolors.ENDC
		#else: print bcolors.HEADER+"master:" +m+bcolors.ENDC

##for coloring processes outputs differently
#red = [255]*256 + range(255,-1,-1) + [0]*256*2 + range(0,256) + [255]*256
#green = red[256*2:] + red[0:256*2]
#blue = green[256*2:] + green[0:256*2]
#rank_color_id = int(len(red)/numThreads*rank)
#rank_color = [str(red[rank_color_id]),str(green[rank_color_id]),str(blue[rank_color_id])]

def mprint(*args): #this one will print regardless of VERBOSE_OUTPUT
	if NO_OUTPUT or rank>MAX_RANK_OUTPUT: return
	m=bcolors[min(rank,len(bcolors)-2)]
	if rank==0:
		m+='1;' #bold for master
	#m+='38;2;'+rank_color[0]+';'+rank_color[1]+';'+rank_color[2]+'m'
	if rank==0:
		m+='Master: '
	else:
		m+='Worker'+str(rank)+": "
	for a in args:
		m+=str(a)
	m+=bcolors[-1]
	print (m)

def wprint(*args):
	if not VERBOSE_OUTPUT: return
	mprint(*args)

# yadeimport.py is generated by `ln -s yade-versionNo yadeimport.py`, where yade-versionNo is the yade executable
#from yadeimport import yade,sphere,box,Sphere,Body,Subdomain,Bo1_Subdomain_Aabb,typedEngine,PFacet,GridConnection,GridNode,PyRunner,kineticEnergy
from yade import *
from yade.utils import *
from yade.wrapper import *
import yade.runtime

mit_mode = yade.runtime.opts.mit>1

def initialize():
	global comm,rank,numThreads, worldComm
	if(mit_mode):
		numThreads=yade.runtime.opts.mit
		process_count = comm.Get_size()
		if not yade.runtime.opts.mpi_mode and process_count<numThreads: #MASTER ONLY
			mprint("I will spawn ",numThreads-process_count," workers")
			if (userScriptInCheckList==""): #normal case
				comm = MPI.COMM_WORLD.Spawn(sys.yade_argv[0], args=sys.yade_argv[1:],maxprocs=numThreads-process_count).Merge()
			else: #HACK, otherwise, handle execution from checkList.py otherwise will we run checkList.py in parallel
				os.environ['OMPI_MCA_rmaps_base_oversubscribe'] = "1" #else spawn fails randomly
				comm = MPI.COMM_WORLD.Spawn(sys.yade_argv[0], args=[userScriptInCheckList],maxprocs=numThreads-process_count).Merge()
			#TODO: if process_count>numThreads, free some workers
			rank=0
		else:	#WORKERS
			comm = MPI.Comm.Get_parent().Merge()
			rank=comm.Get_rank()
			mprint("Hello, I'm worker "+str(rank))
			
		
	else:
		rank = os.getenv('OMPI_COMM_WORLD_RANK')
		numThreads=None
		if rank is not None: #mpiexec was used
			rank=int(rank)
			numThreads=int(os.getenv('OMPI_COMM_WORLD_SIZE'))
		#else monolithic simulation (no mpiexec, no mit)
	return rank,numThreads

def autoInitialize(np):
	global mit_mode
	yade.runtime.opts.mit=np
	mit_mode=True
	return initialize()

def spawnedProcessWaitCommand():
	global waitingCommands
	if waitingCommands: return
	waitingCommands = True
	sys.stderr.write=sys.stdout.write
	mprint("I'm now waiting")
	while 1:
		while not comm.Iprobe(source=0, tag=_MASTER_COMMAND_):
			time.sleep(0.001)
		command = comm.recv(source=0,tag=_MASTER_COMMAND_)
		mprint("I will now execute ",command)
		try:
			exec(command)
		except:
			mprint(sys.exc_info())
		
def sendCommand(executors,command,wait=True):
	'''
	Send a command to a worker (or list of) from master or from another worker
	'''
	start=time.time()
	if not mit_mode: mprint("sendCommand in interactive mode only"); return
	argIsList=isinstance(executors,list)
	if not argIsList: executors = [executors]
	if 0 in executors: mprint("master does not accept mpi commands"); return
	if len(executors)>=numThreads: mprint("executors > numThreads"); return
	
	if wait:
		command="resCommand="+command+";O.subD.comm.send(resCommand,dest="+str(rank)+",tag=_RETURN_VALUE_)"
	for w in executors:
		O.subD.comm.isend(command,dest=w,tag=_MASTER_COMMAND_)
	if wait:
		resCommand= [O.subD.comm.recv(source=w,tag=_RETURN_VALUE_) for w in executors]
		mprint("sendCommand returned in "+str(time.time()-start)+" s")
		return (resCommand if argIsList else res[0])
	else:
		return None


def probeRecvMessage(source, tag):
	msgStat = MPI.Status() 
	O.subD.comm.Probe(source=source, tag=tag, status=msgStat)
	if msgStat.tag == tag : print("message size recvd")
	data = bytearray(msgStat.Get_count(MPI.BYTE))  
	O.subD.comm.Recv([data, MPI.BYTE], source=source, tag=tag)
	return data


class Timing_comm():
	def __init__(self):
		self.timings={}
		
	def clear(self):
		self.timings={}
	
	def print_all(self):
		time.sleep((numThreads-rank)*0.001)
		message = "COMMUNICATION TIMINGS:\n" 
		max_string_len = len(max(self.timings.keys(),key=len))
		for k,v in sorted(self.timings.items(), key=lambda x: x[1][1], reverse=True):
			message += ("{:<"+str(max_string_len)+"}").format(k) + " " + str(v) + "\n"
		mprint(message)
	
	def enable_timing(comm_function):
		def wrapper(self,timing_name, *args, **kwargs):
			#pre-exec
			ti=time.time()
			#exec
			rvalue=comm_function(self, *args, **kwargs)
			#post-exec
			if(not timing_name in self.timings.keys()):
				self.timings[timing_name]=[0,0]
			self.timings[timing_name][0]+=1
			self.timings[timing_name][1]+=time.time()-ti
			return rvalue
		return wrapper
	
	@enable_timing
	def send(self, *args, **kwargs):
		return comm.send(*args,**kwargs)
	
	@enable_timing
	def bcast(self, *args, **kwargs):
		return comm.bcast(*args,**kwargs)
	
	@enable_timing
	def allreduce(self, *args, **kwargs):
		return comm.allreduce(*args,**kwargs)
	
	@enable_timing
	def Gather(self, *args, **kwargs):
		return comm.Gather(*args, **kwargs)
	
	@enable_timing
	def Gatherv(self, *args, **kwargs):
		return comm.Gatherv(*args, **kwargs)
	
	@enable_timing
	def Allgather(self, *args, **kwargs):
		return comm.Allgather(*args, **kwargs)

timing_comm = Timing_comm()

def receiveForces(subdomains):
	'''
	Accumulate forces from subdomains (only executed by master process), should happen after ForceResetter but before Newton and before any other force-dependent engine (e.g. StressController), could be inserted via yade's pyRunner.
	'''
	if 1: #non-blocking:
		reqForces=[]
		for sd in subdomains:#would be better as a loop on subdomains directly, but we don't have those
			
			#wprint( "master getting forces from "+str(b.subdomain)+"(id="+str(b.id)+")")		
			reqForces.append(comm.irecv(None,sd, tag=_FORCES_))
			#wprint( "master got forces from "+str(b.subdomain)+": "+str(forces))		
		for r in reqForces:
			forces=r.wait()
			for ft in forces:
				#wprint(  "adding force "+str(ft[1])+" to body "+str(ft[0]))
				O.forces.addF(ft[0],ft[1])
				O.forces.addT(ft[0],ft[2])
	else:
		for sd in subdomains:
			forces=comm.recv(source=sd, tag=_FORCES_)
			#wprint( "master got forces from "+str(sd)+": "+str(forces)+" iter="+str(O.iter)+" dt="+str(O.dt))
			for ft in forces:
				#wprint(  "adding force "+str(ft[1])+" to body "+str(ft[0]))
				O.forces.addF(ft[0],ft[1])
				O.forces.addT(ft[0],ft[2])
				
def checkNeedCollide():
	'''
	return true if collision detection needs activation in at least one SD, else false. If COPY_MIRROR_BODIES_WHEN_COLLIDE run collider when needed, and in that case return False.
	'''
	needsCollide = int(utils.typedEngine("InsertionSortCollider").isActivated())
	if(needsCollide!=0):
		wprint("triggers collider at iter "+str(O.iter))
	needsCollide = timing_comm.allreduce("checkcollider",needsCollide,op=MPI.SUM)
	if needsCollide:
		if(COPY_MIRROR_BODIES_WHEN_COLLIDE):
			updateMirrorIntersections()
			return False
		else: return True
	return False

def unboundRemoteBodies():
	'''
	Turn bounding boxes on/off depending on rank
	'''
	for b in O.bodies:# unbound the bodies assigned to workers (not interacting directly with other bodies in master scene)
		if (not (b.isSubdomain or isinstance(b.shape, FluidDomainBbox)) and b.subdomain!=rank):
			b.bounded=False
			
def reboundRemoteBodies(ids):
	'''
	update states of bodies handled by other workers, argument 'states' is a list of [id,state] (or [id,state,shape] conditionnaly)
	'''
	if isinstance(ids,list):
		for id in ids:
			b = O.bodies[id]
			if b and not isinstance(b.shape,GridNode): b.bounded=True
	else: #when passing numpy array we need to convert 'np.int32' to 'int'
		for id in ids:
			b = O.bodies[id.item()] 
			if b and not isinstance(b.shape,GridNode): b.bounded=True 

def updateDomainBounds(subdomains): #subdomains is the list of subdomains by body ids
	'''
	Update bounds of current subdomain, broadcast, and receive updated bounds from other subdomains
	Precondition: collider.boundDispatcher.__call__() 
	'''
	wprint( "Updating bounds: "+str(subdomains))
	if(rank==0):
		send_buff=np.zeros(6)*np.nan
	else:
		subD=O.bodies[subdomains[rank-1]].shape #shorthand to shape of current subdomain
		send_buff=np.append(subD.boundsMin,subD.boundsMax)
	recv_buff = np.empty(6*numThreads)
	timing_comm.Allgather("updateDomainBounds",send_buff,recv_buff)
	
	for r in range(1,numThreads):
		O.bodies[subdomains[r-1]].shape.boundsMin = recv_buff[6*r:6*r+3]
		O.bodies[subdomains[r-1]].shape.boundsMax = recv_buff[3+6*r:6+6*r]
		#if(VERBOSE_OUTPUT):#condition here to avoid concatenation overhead
			#mprint("Updated ", O.bodies[subdomains[r-1]].subdomain, " with min=", O.bodies[subdomains[r-1]].shape.boundsMin," and max=", O.bodies[subdomains[r-1]].shape.boundsMax)

def maskedPFacet(pf, boolArray):
	'''
	List bodies within a facet selectively, the ones marked 'True' in boolArray (i.e. already selected from another facet) are discarded
	'''
	l=[]
	for id in [pf.node1.id, pf.node2.id, pf.node3.id, pf.conn1.id, pf.conn2.id, pf.conn3.id]:
		if not boolArray[id]:
			l.append(id)
			boolArray[id]=True

def maskedPFacet(b, boolArray):
	'''
	List bodies within a facet selectively, the ones marked 'True' in boolArray (i.e. already selected from another facet) are discarded
	'''
	l=[]
	pf=b.shape
	for id in [b.id,pf.node1.id, pf.node2.id, pf.node3.id, pf.conn1.id, pf.conn2.id, pf.conn3.id]:
		if not boolArray[id]:
			l.append(id)
			boolArray[id]=True
	return l

def maskedConnection(b, boolArray):
	'''
	List bodies within a facet selectively, the ones marked 'True' in boolArray (i.e. already selected from another facet) are discarded
	'''
	l=[]
	pf=b.shape
	for id in [b.id,pf.node1.id, pf.node2.id]:
		if not boolArray[id]:
			l.append(id)
			boolArray[id]=True
	return l

def genLocalIntersections(subdomains):
	'''
	Defines sets of bodies within current domain overlapping with other domains.
	The structure of the data for domain 'k' is:
	[[id1, id2, ...],  <----------- intersections[0] = ids of bodies in domain k interacting with master domain (subdomain k itself excluded)
	 [id3, id4, ...],  <----------- intersections[1] = ids of bodies in domain k interacting with domain rank=1 (subdomain k itself excluded)
	 ...
	 [domain1, domain2, domain3, ...], <---------- intersections[k] = ranks (not ids!) of external domains interacting with domain k
	 ...
	 ]
	'''
	intersections=[[] for n in range(numThreads)]
	for sdId in subdomains:
		#grid nodes or grid connections could be appended twice or more, as they can participate in multiple pfacets and connexions
		#this bool list is used to append only once
		if (ENABLE_PFACETS):
			if rank==0: print("this one makes mpi inneficient when ENABLE_PFACETS, contact yade devs if you are interested in developping MPI for PFacets")
			appended = np.repeat([False],len(O.bodies))
		subdIdx=O.bodies[sdId].subdomain
		intrs=O.interactions.withBodyAll(sdId)
		#special case when we get interactions with current domain, only used to define interactions with master, otherwise some intersections would appear twice
		if subdIdx==rank:
			for i in intrs:
				otherId=i.id1 if i.id2==sdId else i.id2
				b=O.bodies[otherId]
				if not b:continue #in case the body was deleted
				if b.subdomain==0:
					if isinstance(b.shape,PFacet):
						intersections[0]+= maskedPFacet(b, appended); continue
					if isinstance(b.shape,GridConnection):
						intersections[0]+=maskedConnection(b, appended); continue
					#else (standalone body, normal case)
					intersections[0].append(otherId)
			if len(intersections[0])>0: intersections[subdIdx].append(0)
			continue
		# normal case
		for i in intrs:
			otherId=i.id1 if i.id2==sdId else i.id2
			b=O.bodies[otherId]
			if not b:continue #in case the body was deleted
			if b.subdomain!=rank: continue
			if b.isSubdomain : intersections[rank].append(subdIdx) #intersecting subdomain (will need to receive updated positions from there)
			else:
				if isinstance(b.shape,PFacet):
						intersections[subdIdx]+= maskedPFacet(b, appended); continue
				if isinstance(b.shape,GridConnection):
						intersections[subdIdx]+=maskedConnection(b, appended); continue
				#else (standalone body, normal case)
				intersections[subdIdx].append(otherId)
				
		#for master domain set list of interacting subdomains (could be handled above but for the sake of clarity complex if-else-if are avoided for now)
		if rank==0 and len(intersections[subdIdx])>0:
			intersections[0].append(subdIdx)
	#wprint( "found "+str(len(intrs))+" intersections"+str(intersections))
	return intersections

def updateRemoteStates(states, setBounded=False):
	'''
	update states of bodies handled by other workers, argument 'states' is a list of [id,state] (or [id,state,shape] conditionnaly)
	'''
	ids=[]
	for bst in states:
		#print bst[0],O.bodies[bst[0]]
		ids.append(bst[0])
		b=O.bodies[bst[0]]
		b.state=bst[1]
		#if SEND_SHAPES: b.shape=bst[2]
		if setBounded and not isinstance(b.shape,GridNode): b.bounded=True 
	return ids

def genUpdatedStates(b_ids):
	'''
	return list of [id,state] (or [id,state,shape] conditionnaly) to be sent to other workers
	'''
	return [[id,O.bodies[id].state] for id in b_ids] if not SEND_SHAPES else [[id,O.bodies[id].state,O.bodies[id].shape] for id in b_ids]



#############   COMMUNICATIONS   ################"

def sendRecvStates():
	#____1. get ready to receive positions from other subdomains
	
	pstates = []
	buf = [] #heuristic guess, assuming number of intersecting is ~linear in the number of rows, needs
		
	if rank!=0: #the master process never receive updated states (except when gathering)
		for otherDomain in O.subD.intersections[rank]:
			#wprint( str(": getting states from ")+str(otherDomain))
			if not USE_CPP_MPI:
				buf.append(bytearray(1<<22)) #FIXME: smarter size? this is for a few thousands states max (empirical); bytearray(1<<24) = 128 MB 
				pstates.append( comm.irecv(buf[-1],otherDomain, tag=_ID_STATE_SHAPE_))  #warning leaving buffer size undefined crash for large subdomains (MPI_ERR_TRUNCATE: message truncated)
			else:
				O.subD.mpiIrecvStates(otherDomain) #use yade's messages (coded in cpp)
				#wprint("Receivers set for ",otherDomain)
		#mprint("prepared receive: "+str(time.time()-start)); start=time.time()
	
	#____2. broadcast new positions (should be non-blocking if n>2, else lock) - this includes subdomain bodies intersecting the current one	
	reqs=[]
	for k in O.subD.intersections[rank]:
		if k==rank or k==0: continue #don't broadcast to itself... OTOH this list intersections[rank] will be used to receive
		#if len(b_ids)>0:#skip empty intersections, it means even the bounding boxes of the corresponding subdomains do not overlap
		#mprint("sending "+str(len(O.subD.intersections[k]))+" states to "+str(k))
		if not OPTIMIZE_COM:
			comm.send(genUpdatedStates(O.subD.intersections[k]), dest=k, tag=_ID_STATE_SHAPE_) #should be non-blocking if n>2, else lock?
		else:
			if not USE_CPP_MPI:
				reqs.append(comm.isend(O.subD.getStateValues(k), dest=k, tag=_ID_STATE_SHAPE_)) #should be non-blocking if n>2, else lock?
			else:
				#wprint("Send state to ",k)
				O.subD.mpiSendStates(k)
				#wprint("Sent state to ",k)
		for r in reqs: r.wait() #empty if USE_CPP_MPI
		
	#____3. receive positions and update bodies
	
	if rank==0: return #positions sent from master, done. Will receive forces instead of states
	if not USE_CPP_MPI:
		nn=0	
		for ss in pstates:
			states=ss.wait()
			if not OPTIMIZE_COM:
				updateRemoteStates(states)
			else:
				O.subD.setStateValuesFromIds(O.subD.mirrorIntersections[O.subD.intersections[rank][nn]],states)
				nn+=1
	else:
		for otherDomain in O.subD.intersections[rank]:
			#mprint("listen cpp from "+str(otherDomain)+" in "+str(O.subD.intersections[rank]))
			#subD.mpiRecvStates(otherDomain)
			O.subD.mpiWaitReceived(otherDomain)
			O.subD.setStateValuesFromBuffer(otherDomain)


def isendRecvForces():
	'''
	Communicate forces from subdomain to master
	Warning: the sending sides (everyone but master) must wait() the returned list of requests
	'''	
	O.freqs=[] #keep that one defined even if empty, it is accessed in other functions
	if ACCUMULATE_FORCES:
		if rank!=0:
			forces0=[[id,O.forces.f(id),O.forces.t(id)] for id in  O.subD.mirrorIntersections[0]]
			#wprint ("worker "+str(rank)+": sending "+str(len(forces0))+" "+str("forces to 0 "))
			#O.freqs.append(comm.isend(forces0, dest=0, tag=_FORCES_))
			comm.send(forces0, dest=0, tag=_FORCES_)
		else: #master
			receiveForces(O.subD.intersections[0])

def waitForces():
	'''
	wait until all forces are sent to master. 
	O.freqs is empty for master, and for all threads if not ACCUMULATE_FORCES
	'''
	for r in O.freqs: r.wait()


##### INITIALIZE MPI #########

# Flag used after import of this module, turned True after scene is distributed
O.splitted=False
O.splittedOnce=False #after the first split we have additional bodies (Subdomains) and engines in the merged scene, use this flag to know

def mergeScene():
	if O.splitted:
		if MERGE_W_INTERACTIONS or ERASE_REMOTE_MASTER or DISTRIBUTED_INSERT:
			O.subD.mergeOp()
			sendRecvStatesRunner.dead = isendRecvForcesRunner.dead = waitForcesRunner.dead = collisionChecker.dead = True
			O.splitted=False
			collider.doSort = True
			global NUM_MERGES; NUM_MERGES +=1; 
		else:
			if rank>0:
				# Workers
				send_buff=np.asarray(O.subD.getStateBoundsValuesFromIds([b.id for b in O.bodies if b.subdomain==rank]))
				size=np.array(len(send_buff),dtype=int)
			else:
				#Master
				send_buff=np.array([0])
				size=np.array(0,dtype=int)
			sizes=np.empty(numThreads,dtype=int)
			# Master get sizes from all workers
			timing_comm.Gather("mergeScene_sizes",size,sizes,root=0)
			
			
			if(rank==0):
				# MASTER
				# Alloc sizes for workers 
				dat=np.ones(sizes.sum(),dtype=np.float64)
				# Displacement indexes where data should be stored/received in targeted array
				# dspl should be visible by everyone
				dspl=np.empty(numThreads, dtype=int)
				dspl[0] = 0
				for i in range(1, len(sizes)):
					dspl[i] = dspl[i-1] + sizes[i-1];
			else:
				dspl=None
				dat=None

			# data sent = [data, size of data] (for each worker)
			# data recv = [allocated target_array, array of different sizes, displacement, data type]
			timing_comm.Gatherv("mergeScene_data",[send_buff, size], [dat, sizes, dspl, MPI.DOUBLE], root=0)
			if(rank==0): #master
				for worker_id in range(1, numThreads):
					# generate corresponding ids (order is the same for both master and worker)
					ids = [b.id for b in O.bodies if b.subdomain==worker_id]
					shift = dspl[worker_id];
					if (worker_id != numThreads-1):
						shift_plus_one = dspl[worker_id+1];
					else:		
						shift_plus_one = len(dat);
					O.subD.setStateBoundsValuesFromIds(ids,dat[shift: shift_plus_one]);
					reboundRemoteBodies(ids)
			# turn mpi engines off
			sendRecvStatesRunner.dead = isendRecvForcesRunner.dead = waitForcesRunner.dead = collisionChecker.dead = True
			O.splitted=False
			collider.doSort = True

def splitScene(): 

	'''
	Split a monolithic scene into distributed scenes on threads
	precondition: the bodies have subdomain no. set in user script
	'''
	if not COPY_MIRROR_BODIES_WHEN_COLLIDE: mprint("COPY_MIRROR_BODIES_WHEN_COLLIDE=False is not supported")
	if not ERASE_REMOTE: mprint("ERASE_REMOTE=False is not supported")
	if not O.splittedOnce:
		if DOMAIN_DECOMPOSITION: #if not already partitionned by the user we partition here
			if rank == 0:
				decomposition = dd.decompBodiesSerial(comm) 
				decomposition.partitionDomain(fibreList) 
		if rank == 0 or DISTRIBUTED_INSERT: 
			O.subD=Subdomain() #for storage only, this one will not be used beyond that 
			subD= O.subD #alias
			#insert "meta"-bodies
			subD.subdomains=[] #list subdomains by body ids
			if mit_mode or commSplit : O.subD.comm=comm #make sure the c++ uses the merged intracommunicator
			
			for k in range(1,numThreads):
				domainBody=Body(shape=Subdomain(ids=[b.id for b in O.bodies if b.subdomain==k]),subdomain=k) #note: not clear yet how shape.subDomainIndex and body.subdomain should interact, currently equal values
				domainBody.isSubdomain=True
				subD.subdomains.append(O.bodies.append(domainBody))
				
			O.subdomainIds = subD.subdomains 
			O.thisSubdomainId = 0 

			#tell the collider how to handle this new thing
			collider.boundDispatcher.functors=collider.boundDispatcher.functors+[Bo1_Subdomain_Aabb()]
			if FLUID_COUPLING: 
				collider.boundDispatcher.functors = collider.boundDispatcher.functors+[Bo1_FluidDomainBbox_Aabb()]
			collider.targetInterv=0
			collider.keepListsShort=True # probably not needed, O.bodies.insertAtId should turn it on automaticaly 
			O.bodies.useRedirection=True # idem
			O.bodies.allowRedirection=False
			
			#BEGIN Garbage (should go to some init(), usually done in collider.__call__() but in the mpi case we want to collider.boundDispatcher.__call__() before collider.__call__()
			collider.boundDispatcher.sweepDist=collider.verletDist;
			collider.boundDispatcher.minSweepDistFactor=collider.minSweepDistFactor;
			collider.boundDispatcher.targetInterv=collider.targetInterv;
			collider.boundDispatcher.updatingDispFactor=collider.updatingDispFactor;
			#END Garbage
		if not DISTRIBUTED_INSERT: #we send scene from master to all workers
			sceneAsString= O.sceneToString() if rank==0 else None
			sceneAsString=comm.bcast(sceneAsString,root=0)
			if rank > 0: 
				O.stringToScene(sceneAsString) #receive a scene pre-processed by master (i.e. with appropriate body.subdomain's)  
				# as long as subD.subdomains isn't serialized we need to rebuild it here since it's lost
				domainBody=None
				subdomains=[] #list of subdomains by body id
				for b in O.bodies:
					if b.isSubdomain:
						subdomains.append(b.id)
						if b.subdomain==rank: domainBody=b
				if domainBody==None: wprint("SUBDOMAIN NOT FOUND FOR RANK=",rank)
				O.subD = domainBody.shape
				O.subD.subdomains = subdomains
				
		O._sceneObj.subdomain = rank
		if mit_mode or commSplit : O.subD.comm=comm #make sure the c++ uses the merged intracommunicator
		
		O.subD.init() 
		
		
		if FLUID_COUPLING:
			fluidCoupling = utils.typedEngine('FoamCoupling') 
			fluidCoupling.comm = comm 
			fluidCoupling.getFluidDomainBbox() #triggers the communication between yade procs and Yales2/openfoam procs, get's fluid domain bounding boxes from all yales2 procs. 
			fluidCoupling.setIdList(fluidBodies) 
			#if (fluidBodies) :  # incase fluidBodies are not set to the engine directly in the user script. 
			
		updateMirrorIntersections()
		
		idx = O.engines.index(utils.typedEngine("NewtonIntegrator"))
		O.engines=O.engines[:idx+1]+[PyRunner(iterPeriod=1,initRun=True,command="sys.modules['yade.mpy'].sendRecvStates()",label="sendRecvStatesRunner")]+O.engines[idx+1:]
		
		# append force communicator before Newton
		O.engines=O.engines[:idx]+[PyRunner(iterPeriod=1,initRun=True,command="sys.modules['yade.mpy'].isendRecvForces()",label="isendRecvForcesRunner")]+O.engines[idx:]
		
		# append engine waiting until forces are effectively sent to master
		O.engines=O.engines+[PyRunner(iterPeriod=1,initRun=True,command="pass",label="waitForcesRunner")]
		O.engines=O.engines+[PyRunner(iterPeriod=1,initRun=True,command="if sys.modules['yade.mpy'].checkNeedCollide(): O.pause()",label="collisionChecker")]
		
		O.splitted=True
		O.splittedOnce = True
	
	else: 
		if (DOMAIN_DECOMPOSITION and RESET_SUBDOMAINS_WHEN_COLLIDE):
			if rank == 0:
				decomposition = dd.decompBodiesSerial(comm) 
				decomposition.partitionDomain() 
			O.subD.splitBodiesToWorkers(RESET_SUBDOMAINS_WHEN_COLLIDE)
			updateMirrorIntersections()
		if rank == 0 :
			O.interactions.clear()
			unboundRemoteBodies()
			if (ERASE_REMOTE and ERASE_REMOTE_MASTER): eraseRemote()
		sendRecvStatesRunner.dead = isendRecvForcesRunner.dead = waitForcesRunner.dead = collisionChecker.dead = False
		O.splitted = True
		
		
bodiesToImport=[]

def updateMirrorIntersections():
	global bodiesToImport
	subD=O.subD
	start = time.time()
	if (not O.splitted):
		unboundRemoteBodies()
		eraseRemote()
	collider.boundDispatcher.__call__()
	updateDomainBounds(subD.subdomains) #triggers communications
	collider.__call__() #see [1]
	unboundRemoteBodies() #in splitted stage we exploit bounds to detect bodies which are no longer part of intersections (they will be left with no bounds after what follows)

	subD.intersections=genLocalIntersections(subD.subdomains)
	#update mirror intersections so we know message sizes in advance
	subD.mirrorIntersections=[[] for n in range(numThreads)]
	if rank==0:#master domain
		for worker in range(1,numThreads):#FIXME: we actually don't need so much data since at this stage the states are unchanged and the list is used to re-bound intersecting bodies, this is only done in the initialization phase, though
			#wprint("sending mirror intersections to "+str(worker)+" ("+str(len(subD.intersections[worker]))+" bodies), "+str(subD.intersections[worker]))
			m = O.intrsctToBytes(subD,worker,False) if SEND_BYTEARRAYS else  subD.intersections[worker];
			timing_comm.send("sendIntersections",m, dest=worker, tag=_MIRROR_INTERSECTIONS_)
	else:
		# from master
		b_ids=comm.recv(source=0, tag=_MIRROR_INTERSECTIONS_)
		wprint("Received mirrors from master: ",b_ids)
		#FIXME: we are assuming that Body::id_t is 4 bytes here, not that portable...
		numInts0= int(len(b_ids)/4) if SEND_BYTEARRAYS else len(b_ids)  #ints = 4 bytes
		
		if numInts0>0:
			if SEND_BYTEARRAYS:
				O.bufferFromIntrsct(subD,0,numInts0,True)[:]=b_ids
				b_ids=np.frombuffer(b_ids,dtype=np.int32)
			else:
				subD.mirrorIntersections= [b_ids]+subD.mirrorIntersections[1:]

			reboundRemoteBodies(b_ids)
			# since interaction with 0-bodies couldn't be detected before, mirror intersections from master will
			# tell if we need to wait messages from master (and this is declared via intersections) 
			if not 0 in subD.intersections[rank]:
				temp=subD.intersections[rank]
				temp+=[0]
				subD.intersections=subD.intersections[:rank]+[temp]+subD.intersections[rank+1:]
			else:
				if not O.splittedOnce: mprint("0 already in intersections (should not happen)")
		reqs=[]
		
		#from workers
		for worker in subD.intersections[rank]:
			if worker==0: continue #already received above
			#wprint("subD.intersections["+str(rank)+"]: "+str(subD.intersections[rank]))
			#buf = bytearray(1<<22) #CRITICAL
			reqs.append([worker,comm.irecv(None, worker, tag=_MIRROR_INTERSECTIONS_)])

		for worker in subD.intersections[rank]:
			if worker==0: continue #we do not send positions to master, only forces
			#wprint("sending "+str(len(subD.intersections[worker]))+" states to "+str(worker))
			wprint("Send mirrors to: ", worker)
			m = O.intrsctToBytes(subD,worker,False) if SEND_BYTEARRAYS else  subD.intersections[worker];
			timing_comm.send("splitScene_intersections", m, dest=worker, tag=_MIRROR_INTERSECTIONS_)
		for req in reqs:
			intrs=req[1].wait()
			wprint("Received mirrors from: ", req[0], " : ",np.frombuffer(intrs,dtype=np.int32))
			if SEND_BYTEARRAYS:
				O.bufferFromIntrsct(subD,req[0],int(len(intrs)/4),True)[:]=intrs
				intrs=np.frombuffer(intrs,dtype=np.int32)
			else:
				subD.mirrorIntersections = subD.mirrorIntersections[0:req[0]]+[intrs]+subD.mirrorIntersections[req[0]+1:]
			#reboundRemoteBodies(intrs)
			
	if rank!=0:
		for l in subD.mirrorIntersections:
			if len(l)>0:
				reboundRemoteBodies(l)
	if(ERASE_REMOTE): eraseRemote() # erase the bodies which still have no bounds
		
	"""
	" NOTE: FK, what to do here:
	" 1- all threads loop on reqs, i.e the intersecting subdomains of the current subdomain.
	" 2- during this loop, check whether the current subdomain needs bodies from the intersecting ones
	" 3- in all cases, isend the ids needed, (if no ids needed send empty array)
	" 3.1- build a ranks list "requestedSomethingFrom" to loop later on to receive data
	" 4- loop again on reqs to get the ids needed by other subdomains (with blocking recv as we used isend)
	" 5- if the data recved is empty (nothing requested), do nothing. Else isend the bodies (c++)
	" 6- loop on "requestedSomethingFrom" ranks and recv the bodies (blocking, c++, using MPI_Probe to know the message size)
	" 7- comm.barrier(), just in case
	"""
	if(COPY_MIRROR_BODIES_WHEN_COLLIDE or MERGE_W_INTERACTIONS):
		requestedSomethingFrom=[]
		bodiesToImport=[[] for worker in range(numThreads)]
		sent=[]
		if rank>0: #master doesn't need bodies
			#reqs=subD.intersections[rank]
			for worker in subD.intersections[rank]:
				#worker=req[0]
				#if(worker==0):continue
				for mirrorBodyId in subD.mirrorIntersections[worker]:
					if O.bodies[mirrorBodyId]==None:
						bodiesToImport[worker]+=[mirrorBodyId]
				if(len(bodiesToImport[worker])>0):
					requestedSomethingFrom.append(worker)
				wprint("I request ids ",bodiesToImport[worker], " from ",worker)
				sent.append(comm.isend(bodiesToImport[worker], worker, tag=_MIRROR_INTERSECTIONS_))
		wprint("will wait requests from ",subD.intersections[rank])
		for worker in subD.intersections[rank]:
			if worker!=0:
				wprint("waiting requests from ",worker)
				requestedIds=comm.recv(source=worker, tag=_MIRROR_INTERSECTIONS_)
				if(len(requestedIds)>0):
					wprint("I will now send ",requestedIds," to ",worker)
					subD.sendBodies(worker,requestedIds)
		
		for worker in requestedSomethingFrom:
			subD.receiveBodies(worker)
		for s in sent: s.wait()
		subD.completeSendBodies();
	if not collider.keepListsShort: collider.doSort = True
	collider.__call__()
	collider.execTime+=int((time.time()-start)*1e9)
	collider.execCount+=1
	try:
		collisionChecker.execTime-=int((time.time()-start)*1e9)
	except:
		pass

	#maxVelocitySq is normally reset in NewtonIntegrator in the same iteration as bound dispatching, since Newton will not run before next iter in our case we force that value to avoid another collision detection at next step
	utils.typedEngine("NewtonIntegrator").maxVelocitySq=0.5


def eraseRemote(): 
	if rank>0 or ERASE_REMOTE_MASTER: # suppress external bodies from scene
		#numBodies = len(O.bodies)
		#for id in range(numBodies):
		#mprint("will erase ",[b.id for b in O.bodies if (not b.bounded and b.subdomain!=rank)])
		for b in O.bodies:
			if not b.bounded and b.subdomain != rank:
				connected = False #a gridNode could be needed as part of interacting facet/connection even if not overlaping a specific subdomain. Assume connections are always bounded for now, we thus only check nodes.
				if isinstance(b.shape,GridNode):
					for f in b.shape.getPFacets():
						if f.bounded: connected = True
					for c in b.shape.getConnections():
						if c.bounded: connected = True
				if not connected:
					O.bodies.erase(b.id)


##### RUN MPI #########
def mpirun(nSteps,np=numThreads,withMerge=False):
	'''
	Parallel version of O.run() using MPI domain decomposition.
	
	Parameters
        ----------
        nSteps : int
            The numer of steps to compute
        np : int
            number of iterations
        withMerge : bool
            wether subdomains should be merged into master at the end of the run (default False). If True the scene in the master process is exactly in the same state as after O.run(nSteps,True). The merge can be time consumming, it is recommended to activate only if post-processing or other similar tasks require it.
	'''
	stack=inspect.stack()
	global userScriptInCheckList, LOAD_SIM
	if len(stack[3][1])>12 and stack[3][1][-12:]=="checkList.py":
		userScriptInCheckList=stack[1][1]
	caller_name = stack[2][3]
	
	if (np>numThreads):  
		if not mit_mode: autoInitialize(np) #this will set numThreads
		else: mprint("number of cores can't be increased after first call to mpirun")
	if(mit_mode and rank==0 and not caller_name=='execfile'): #if the caller is the user's script, everyone already calls mpirun and the workers are not waiting for a command.
		for w in range(1,numThreads):
			comm.send("yade.mpy.mpirun(nSteps="+str(nSteps)+",withMerge="+str(withMerge)+")",dest=w,tag=_MASTER_COMMAND_)
			wprint("Command sent to ",w)
	initStep = O.iter
	if not O.splitted:
		wprint("splitting")
		splitScene()
		wprint("splitted")
	if YADE_TIMING:
		O.timingEnabled=True
	if not (MERGE_SPLIT):
		O.run(nSteps,True) #a pyrunner will pause, or trigger collider and continue, depending on WITH_BODY_COPY flag
		if withMerge: mergeScene() #will be useful to see evolution in QGLViewer, for instance
	else: #merge/split or body_copy for each collider update
		collisionChecker.dead=True
		while (O.iter-initStep)<nSteps:
			O.step()
			if checkNeedCollide():
				mergeScene()
				splitScene()
		mergeScene()
	if YADE_TIMING and rank<=MAX_RANK_OUTPUT:
		timing_comm.print_all()
		from yade import timing
		time.sleep((numThreads-rank)*0.002) #avoid mixing the final output, timing.stats() is independent of the sleep
		mprint( "#####  Worker "+str(rank)+"  ######")
		timing.stats() #specific numbers for -n4 and gabion.py


def sendTerminateMessage(): 
	'''Used for coupled mpi simulations to notify non yade-mpi processes of termination, (they call mpi.finalize after this broadcast) , e.g: yade-openfoam, yade-yales2 etc'''  
	data = numpy.arange(3, dtype = 'i') 
	worldRank = worldComm.Get_rank()
	worldComm.Bcast([data,MPI.INT], root=0)


def killMPI(): 
	MPI.Finalize()
