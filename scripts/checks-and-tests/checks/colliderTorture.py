# bruno.chareyre@grenoble-inp.fr
# This script checks collider correctness and performance in bouncing spheres with heavy erase/insert along with iterations

import time

newton.gravity = (0, -10, 0)
newton.damping = 0

O.engines=[ #this reproduces the yade default, not present during checkList.py execution
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb(),Bo1_Box_Aabb()],verletDist=0.3,label="collider"),
        InteractionLoop(
            [Ig2_Facet_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()], # <====== disable sphere-sphere contacts, they will be virtual interactions
            [Ip2_FrictMat_FrictMat_FrictPhys()],
            [Law2_ScGeom_FrictPhys_CundallStrack(label="law")],
            label="interactionLoop"
        ),
        NewtonIntegrator(label="newton")
    ]

O.dynDt = False
O.dt = 3e-3

O.bodies.append(box(center=[0, 0, 0], extents=[500, 0, 500], color=[0, 0, 1], fixed=True))
N = 40
# bunch of non-interacting spheres crossing each other in the air
O.bodies.append(sphere((N / 3, 6, N / 3), 6))
O.bodies.append(sphere((6 + N / 3, 12, 6 + N / 3), 6))
for i in range(N):
	for j in range(N):
		O.bodies.append(sphere((i, 2 + cos(0.6 * (i + j)), j), 1))
		O.bodies.append(sphere((i + 0.1, 6 + 1.5 * cos(2 * (i + j) + 1), j + 0.3), 1))

O.saveTmp()

buffer = []


def erase(chunkSize):
	chunk = [b for b in O.bodies][-chunkSize:]
	buffer.append(chunk)
	for b in chunk:
		id = b.id
		O.bodies.erase(b.id)
		if not isinstance(b.shape, Sphere):
			print("problem with", id, "you erased too much")


def insert(chunkSize):
	global buffer
	if len(buffer) < 1:
		print("problem")
	for b in buffer[-1]:
		O.bodies.append(Body(shape=b.shape, state=b.state, material=b.material))  #clone, we can't rei-insert
	buffer = buffer[:-1]


#by default 1/3 of the total removed incrementaly then inserted
def cycle(maxBufferSize=int(N / 3), chunkSize=int(3 * N), doCheck=True, iterations=1):
	global buffer
	times = []
	t1 = time.time()
	O.run(iterations, True)
	times.append(time.time() - t1)
	if doCheck:
		print("== init run ==")
		check()
	t1 = time.time()
	while len(buffer) < maxBufferSize:
		erase(chunkSize)
		O.run(iterations, True)
	times.append(time.time() - t1)
	if doCheck:
		print("== delete ==")
		check()
	t1 = time.time()
	while len(buffer) > 0:
		insert(chunkSize)
		#buffer=buffer[:-1]
		O.run(iterations, True)
	times.append(time.time() - t1)
	if doCheck:
		print("== insert ==")
		check()
	return sum(times)


ints1 = []
ints2 = []


def check():
	global ints1, ints2
	collider.__call__()
	sig1 = signature()
	ints1 = O.interactions.all()
	collider.doSort = True
	collider.__call__()
	sig2 = signature()
	ints2 = O.interactions.all()

	if sig1 == sig2:
		print("signatures match", sig1)
	else:
		if O.bodies.enableRedirection:
			raise YadeCheckError("signatures mismatch " + str(sig1) + " vs " + str(sig2))
		else:
			print("WARN: this is a bug with enableRedirection=False: signatures mismatch " + str(sig1) + " vs " + str(sig2))


def signature():
	sig = 0
	count = 0
	for i in O.interactions.all():
		# identify overlapping boxes disregarding verletDist, this is the minimal list that any AABB collider should find
		# it means result of the check will not depend on collider optimizations giving more or less virtual interactions
		d = O.bodies[i.id1].state.pos - O.bodies[i.id2].state.pos
		sumRad = yade.math.Real(0)
		if i.id1 > 0:
			sumRad += O.bodies[i.id1].shape.radius
		if i.id2 > 0:
			sumRad += O.bodies[i.id2].shape.radius
		away = abs(d[0]) >= sumRad or abs(d[1]) >= sumRad or abs(d[2]) >= sumRad

		if (not away or i.isReal):  #isReal is for sphere-box
			sig += i.id1 + i.id2
			count += 1
	return count, sig


def testFraction(fraction, withInsert=True):
	collider.__call__()  #init
	times = []
	n = max(1, int(2 * N * N * fraction))
	erase(n)
	if withInsert:
		insert(n)
	t1 = time.time()
	collider.__call__()
	times.append(time.time() - t1)
	if not withInsert:
		insert(n)
	# now do the timing with initSort
	collider.__call__()  #init
	erase(n)
	if withInsert:
		insert(n)
	t1 = time.time()
	collider.doSort = True
	collider.__call__()
	times.append(time.time() - t1)
	if not withInsert:
		insert(n)
	return times


for enableRedirection in [False, True]:
	print("================")
	print("redirection:", enableRedirection)
	print("================")

	O.reload()
	O.bodies.enableRedirection = enableRedirection
	cycle(10, 100, doCheck=True, iterations=50)  # remove spheres by batches of 100 x 10 times, then insert again, and check corectness
	O.reload()
	O.bodies.enableRedirection = enableRedirection
	O.step()

	Nc = 3
	Niter = 1
	substeps = 2
	start = time.time()
	for k in range(7):
		t = 0
		for k in range(Nc):
			t += cycle(substeps, int(0.4 * N * N * 2 / substeps), doCheck=False, iterations=Niter)
		print("time for ", Nc, " erase/insert cycles in", 2 * Niter * substeps * Nc, "iterations:", t, "sec")
	print("====== TOTAL CYCLING TIME", time.time() - start, "for ", len([b for b in O.bodies]), "real bodies in len(O.bodies)=", len(O.bodies), " =======")
	print("")
	print("Now insert/erase fractions of the total")
	for fraction in [0, 0.001, 0.01, 0.04, 0.1, 0.25, 0.4, 0.7]:
		t = testFraction(fraction, True)
		print("fraction=", fraction, "collider time:", t[0], "vs. initSort", t[1])
		#if (t[0]>3*t[1]): raise YadeCheckError("collider run in more than twice the initSort time, ratio:"+str(times[0]/times[1]))
