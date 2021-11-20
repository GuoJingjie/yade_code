import os
import sys
import pickle
import csv
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

#-------------------#
# Constants
#-------------------#
bigSep = "\n====================== "
sep = "---------------------- "

#-------------------#
# Variables
#-------------------#
d_ad = 0


#-------------------#
# Defining utils
#-------------------#
def color_gradient(grad_nb, p):
	p.colors = []
	for i in range(grad_nb):
		c = (i / float(grad_nb))
		p.colors.append((p.r[0] * c + p.r[1] * (1 - c), p.v[0] * c + p.v[1] * (1 - c), p.b[0] * c + p.b[1] * (1 - c)))


def average(qT, t):
	""" Average qT over time.

	Parameters:
	- qT : Contains the quantity. : List of floats.
	"""
	n_time = len(t) - 1
	# Finding the first value to take into account.
	i_deb = 0
	while i_deb < len(t) and t[i_deb] < pPP.mean_begin_time:
		i_deb += 1
	if i_deb > n_time - 1:
		print('WARNING average: End of simulation before start of averaging.')
		print('WARNING average: Taking only last profile.')
		i_deb = n_time - 1
	# Initialisation
	q = qT[i_deb]
	# Averaging
	i = i_deb + 1
	while i < n_time + 1 and t[i] < pPP.mean_end_time:
		q += qT[i]
		i += 1
	q /= (i - i_deb)
	return q


def average_profile(qT, t, n=False):
	""" Average qT profiles over time.

	Parameters:
	- qT : Contains all the profiles in the time. : List of lists of floats.
	- n : Enable normalise (for histograms) : bool
	"""
	n_time = len(t) - 1
	# Finding the first value to take into account.
	i_deb = 0
	while i_deb < len(t) and t[i_deb] < pPP.mean_begin_time:
		i_deb += 1
	if i_deb > n_time - 1:
		print('WARNING average_profile: End of simulation before start of averaging.')
		print('WARNING average_profile: Taking only last profile.')
		i_deb = n_time - 1
	# Initialisation
	q = qT[i_deb][:]
	# Averaging
	i = i_deb + 1
	while i < n_time + 1 and t[i] < pPP.mean_end_time:
		for j in range(0, len(q)):
			q[j] += qT[i][j]
		i += 1
	if n:
		summ = sum(q)
		q = [v / summ for v in q]
	else:
		q = [v / (i - i_deb) for v in q]
	return q


def ponderate_average_profile(pT, qT, t):
	""" Average qT profiles over time.

	Parameters:
	- qT : Contains all the profiles in the time. : List of lists of floats.
	- n : Enable normalise (for histograms) : bool
	"""
	n_time = len(t) - 1
	# Finding the first value to take into account.
	i_deb = 0
	while i_deb < len(t) and t[i_deb] < pPP.mean_begin_time:
		i_deb += 1
	if i_deb > n_time - 1:
		print('WARNING average_profile: End of simulation before start of averaging.')
		print('WARNING average_profile: Taking only last profile.')
		i_deb = n_time - 1
	# Initialisation
	p = pT[i_deb][:]
	q = []
	print(len(p), " ", len(qT[i_deb]))
	for j in range(len(p)):
		q.append(p[j] * qT[i_deb][j])
	# Averaging
	i = i_deb + 1
	while i < n_time + 1 and t[i] < pPP.mean_end_time:
		for j in range(0, len(q)):
			q[j] += pT[i][j] * qT[i][j]
			p[j] += pT[i][j]
		i += 1

	for j in range(len(q)):
		if p[j] > 0.0:
			q[j] /= p[j]

	return q


def adim(q, star):
	""" non-dimensionalize the quantitie q by star.

	Parameters:
	- q : list
	- star : float
	"""
	return [a / star for a in q]


def integration(phi, y, dx):
	""" Integrate y along x ponderate by phi.

	"""
	q = 0
	for j in range(len(y)):
		q += phi[j] * y[j] * dx

	return q


def computeTransportLayerThickness(qsz, dz):
	""" Computes the thickness of the transport layer based on the law proposed by Duran et al 2012. 

	"""

	# Computing qsat
	qsat = 0
	for j in range(len(qsz)):
		qsat += max(0.0, qsz[j]) * dz

	if qsat > 1e-7:
		# Computing zbar
		zbar = 0
		for j in range(len(qsz)):
			zbar += max(0.0, qsz[j]) * (j * dz) * dz
		zbar /= qsat

		# Computing lambda
		l = 0
		for j in range(len(qsz)):
			l += pow((j * dz) - zbar, 2) * max(0.0, qsz[j]) * dz
		l = sqrt(l / qsat)
	else:
		l = 0.0

	return l


#-------------------#
# Import measures functions
#-------------------#
execfile("params_post_proc.py")
color_gradient(len(sys.argv) - 1, pP1D)
execfile("common/simulationPyRunners.py")
execfile("common/measures.py")


def post_process(dr):
	# Update name_value
	name_value = str(eval(pPP.param))

	# Post Processing
	for p in pP1D.post_process:
		for key in p:
			print(p[key])
			data[key] = eval(p[key])

	### Storing 2D data
	if pP2D.plot_enable:
		bval = eval(pPP.param)
		if not (bval in bdata):
			bdata[bval] = {}
			for key in pP2D.measures:
				bdata[bval][key] = []
		for key in pP2D.measures:
			print(pP2D.measures[key])
			bdata[bval][key].append(eval(pP2D.measures[key]))

	### Ploting 1D data
	if pP1D.plot_enable:
		m = pP1D.markers.pop()
		c = pP1D.colors.pop()
		### Ploting figures
		print(sep + "Ploting data.")
		# Plots
		for key in pP1D.plots:
			errx = None
			erry = None
			if len(pP1D.plots[key]) > 2:
				if pP1D.plots[key][2][0] != "":
					errx = data[pP1D.plots[key][2][0]]
				if pP1D.plots[key][2][1] != "":
					erry = data[pP1D.plots[key][2][1]]
			for i in range(len(pP1D.plots[key][0])):
				x = pP1D.plots[key][0][i]
				y = pP1D.plots[key][1][i]
				me = int(max(1.0, pP1D.me * len(data[x])))
				if pP1D.mec == None:
					mec = c
				else:
					mec = pP1D.mec
				axs[key].errorbar(
				        data[x],
				        data[y],
				        xerr=errx,
				        yerr=erry,
				        color=c,
				        marker=m,
				        markevery=me,
				        markerfacecolor=c,
				        markeredgewidth=pP1D.mew,
				        markeredgecolor=mec,
				        markersize=pP1D.ms,
				        label=r"$" + pPP.name_param + "=" + name_value + "$"
				)
		# PlotsT
		for key in pP1D.plotsT:
			space = pP1D.plotsT[key][2]
			for i in range(len(data["time"])):
				for j in range(len(pP1D.plotsT[key][0])):
					x = pP1D.plotsT[key][0][j]
					y = pP1D.plotsT[key][1][j]
					me = int(max(1.0, pP1D.me * len(data[x])))
					if pP1D.mec == None:
						mec = c
					else:
						mec = pP1D.mec
					if i < 1:
						axsT[key].plot(
						        [v + i * space for v in data[x][i]],
						        data[y],
						        color=c,
						        marker=m,
						        markevery=me,
						        markerfacecolor=c,
						        markeredgewidth=pP1D.mew / len(data["time"]),
						        markeredgecolor=mec,
						        markersize=pP1D.ms / len(data["time"]),
						        label=r"$" + pPP.name_param + "=" + name_value + "$"
						)
					else:
						axsT[key].plot(
						        [v + i * space for v in data[x][i]],
						        data[y],
						        color=c,
						        marker=m,
						        markevery=me,
						        markerfacecolor=c,
						        markeredgewidth=pP1D.mew / len(data["time"]),
						        markeredgecolor=mec,
						        markersize=pP1D.ms / len(data["time"])
						)
		# Orientations
		for key in pP1D.orientations:
			for j in range(len(pP1D.orientations[key][0])):
				xyz = pP1D.orientations[key][0][j]
				X = data[xyz[0]]
				Y = data[xyz[1]]
				Z = data[xyz[2]]
				C_name = pP1D.orientations[key][1][j]
				C = data[C_name]
				#tmp = axsO[key].scatter3D(X, Y, Z, c=C, s=[Vector3(X[j], Y[j], Z[j]).norm()*80.0 for j in range(len(C))], cmap="Greys")
				#if C_name in pPP.plots_names:
				#	tmp = figsO[key].colorbar(tmp)
				#	tmp.ax.set_ylabel(pPP.plots_names[C_name])
				for i in range(len(X)):
					norm = min(Vector3(X[i], Y[i], Z[i]).norm(), 1.0)
					axsO[key].quiver3D(
					        0.0,
					        0.0,
					        0.0,
					        X[i],
					        Y[i],
					        Z[i],
					        colors=[(0, 0, 0, norm), (0, 0, 0, norm), (0, 0, 0, norm)],
					        linewidth=2.0,
					        length=norm
					)


def plot_external_data():
	for ext_key in pP1D.plotsExtPath:
		path = pP1D.plotsExtPath[ext_key]
		# Selecting marker and color
		m = pP1D.markers.pop()
		c = pP1D.colors.pop()
		# Getting data
		data = {}
		# Getting the extension
		ext = path.split(".")[-1]
		if ext == "py":
			execfile(path)
		elif ext == "csv":
			with open(path) as f:
				dialect = csv.Sniffer().sniff(f.read(1024), delimiters=';')
				f.seek(0)
				reader = csv.DictReader(f, dialect=dialect)
				data = {key: [] for key in reader.fieldnames}
				for row in reader:
					for key in reader.fieldnames:
						data[key].append(float(row[key].replace(",", ".")))
		elif ext == "txt":
			with open(path) as f:
				lines = [line.rstrip().split() for line in f.readlines()]
				keys = lines[0]
				data = {key: [] for key in keys}
				for line in lines[1:]:
					for i in range(len(keys)):
						data[keys[i]].append(float(line[i]))
		# Plots Ext
		for key in pP1D.plotsExt:
			for i in range(len(pP1D.plotsExt[key][0])):
				x = pP1D.plotsExt[key][0][i]
				y = pP1D.plotsExt[key][1][i]
				me = int(max(1.0, pP1D.me * len(data[x])))
				if pP1D.mec == None:
					mec = c
				else:
					mec = pP1D.mec
				axs[key].plot(
				        data[x],
				        data[y],
				        color=c,
				        marker=m,
				        markevery=me,
				        markerfacecolor=c,
				        markeredgewidth=pP1D.mew,
				        markeredgecolor=mec,
				        markersize=pP1D.ms,
				        label=ext_key,
				        linestyle='None'
				)


#-------------------#
# Creating 1D Figures
#-------------------#
if pP1D.plot_enable:
	# Plots
	figs = {}
	axs = {}
	for key in pP1D.plots:
		figs[key] = plt.figure()
		axs[key] = plt.gca()
		plt.xlabel(pPP.plots_names[pP1D.plots[key][0][0]])
		plt.ylabel(pPP.plots_names[pP1D.plots[key][1][0]])
	for key in pP1D.alims:
		if pP1D.alims[key][0]:
			axs[key].set_xlim(pP1D.alims[key][0][0], pP1D.alims[key][0][1])
		if pP1D.alims[key][1]:
			axs[key].set_ylim(pP1D.alims[key][1][0], pP1D.alims[key][1][1])
	# PlotsT
	figsT = {}
	axsT = {}
	for key in pP1D.plotsT:
		figsT[key] = plt.figure()
		axsT[key] = plt.gca()
		plt.xlabel(pPP.plots_names[pP1D.plotsT[key][0][0]])
		plt.ylabel(pPP.plots_names[pP1D.plotsT[key][1][0]])
	# Orientations
	figsO = {}
	axsO = {}
	for key in pP1D.orientations:
		figsO[key] = plt.figure()
		axsO[key] = plt.gca(projection='3d')
		axsO[key].set_aspect('equal', 'box')
		axsO[key].set_xlabel("$x$")
		axsO[key].set_ylabel("$y$")
		axsO[key].set_zlabel("$z$")
		axsO[key].set_xticklabels([])
		axsO[key].set_yticklabels([])
		axsO[key].set_zticklabels([])
		axsO[key].w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
		axsO[key].w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
		axsO[key].w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
	for key in pP1D.alimsO:
		if pP1D.alimsO[key][0]:
			axsO[key].set_xlim(pP1D.alimsO[key][0][0], pP1D.alimsO[key][0][1])
		if pP1D.alimsO[key][1]:
			axsO[key].set_ylim(pP1D.alimsO[key][1][0], pP1D.alimsO[key][1][1])
		if pP1D.alimsO[key][2]:
			axsO[key].set_zlim(pP1D.alimsO[key][2][0], pP1D.alimsO[key][2][1])

# Declaring batch data storage
bdata = {}
badata = {}

#-------------------#
# Post Processing 1D
#-------------------#
for dr in sys.argv[1:]:
	print(bigSep + dr)
	os.chdir(dr)
	execfile("params.py")
	os.chdir("..")
	d_ad = eval(pPP.d_ad)
	with open(dr + '/data.pkl', 'rb') as f:
		data = pickle.load(f)
	print(sep + "Post processing 1D.")
	post_process(dr)

#-------------------#
# Post Processing External
#-------------------#
pP1D.r = pPP.ext_r
pP1D.b = pPP.ext_b
pP1D.v = pPP.ext_v
color_gradient(len(pP1D.plotsExtPath), pP1D)
plot_external_data()

#-------------------#
# Post Processing 2D
#-------------------#
if pP2D.plot_enable:
	print(sep + "Post processing 2D, Plotting and Saving.")
	### Post processing
	for bval in bdata:
		for p in pP2D.post_process:
			for key in p:
				print(p[key])
				bdata[bval][key] = eval(p[key])
	### Additional data
	for p in pP2D.add:
		for key in p:
			print(p[key])
			badata[key] = eval(p[key])
	### Sorting data
	params = []
	params_val = []
	for p, v in bdata.items():
		params.append(p)
		params_val.append(v)
	params, params_val = zip(*sorted(zip(params, params_val)))
	for key in pP2D.plots:
		plt.figure()
		plt.xlabel(pPP.plots_names[pP2D.plots[key][0][0]])
		plt.ylabel(pPP.plots_names[pP2D.plots[key][1][0]])
		if pP2D.alims[key][0]:
			plt.xlim(pP2D.alims[key][0][0], pP2D.alims[key][0][1])
		if pP2D.alims[key][1]:
			plt.ylim(pP2D.alims[key][1][0], pP2D.alims[key][1][1])
		# Plotting
		batch_markers = pP2D.markers[:]
		pP2D.r = pPP.n_r
		pP2D.v = pPP.n_v
		pP2D.b = pPP.n_b
		color_gradient(len(params), pP2D)
		# Normal plots
		for i in range(len(params)):
			p = params[i]
			v = params_val[i]
			m = batch_markers.pop()
			c = pP2D.colors.pop()
			for j in range(len(pP2D.plots[key][0])):
				x = pP2D.plots[key][0][j]
				y = pP2D.plots[key][1][j]
				me = int(max(1.0, pP2D.me * len(v[x])))
				lab = ""
				if pPP.name_param != "":
					lab += r'$' + pPP.name_param + "=" + str(p) + "$"
				if len(pP2D.plots[key]) > 2:
					if lab != "":
						lab += " "
					lab += pP2D.plots[key][2][j]

				if pP2D.mec == None:
					mec = c
				else:
					mec = pP2D.mec
				style = "-"
				if j > 0:
					style = ":"
				plt.plot(
				        v[x],
				        v[y],
				        color=(c[0] / (j + 1), c[1] / (j + 1), c[2] / (j + 1)),
				        marker=m,
				        markevery=me,
				        markeredgewidth=pP2D.mew,
				        markeredgecolor=mec,
				        markerfacecolor=c,
				        markersize=pP2D.ms / (j + 1),
				        label=lab,
				        linestyle=style
				)
				plt.legend(fancybox=True, framealpha=0.5, loc='center left', bbox_to_anchor=(1.0, 0.5))

		pP2D.r = pPP.add_r
		pP2D.v = pPP.add_v
		pP2D.b = pPP.add_b
		color_gradient(len(pP2D.plot_adds.keys()), pP2D)
		if (key in pP2D.plot_adds):
			# Additional data plots
			m = batch_markers.pop()
			c = pP2D.colors.pop()
			for j in range(len(pP2D.plot_adds[key][0])):
				x = pP2D.plot_adds[key][0][j]
				y = pP2D.plot_adds[key][1][j]
				me = int(max(1.0, pP2D.me * len(v[x])))
				lab = ""
				if len(pP2D.plot_adds[key]) > 2:
					if lab != "":
						lab += " "
					lab += pP2D.plot_adds[key][2][j]

				if pP2D.mec == None:
					mec = c
				else:
					mec = pP2D.mec
				plt.plot(
				        badata[x],
				        badata[y],
				        color=(c[0] / (j + 1), c[1] / (j + 1), c[2] / (j + 1)),
				        marker=m,
				        markevery=me,
				        markeredgewidth=pP2D.mew,
				        markeredgecolor=mec,
				        markerfacecolor=c,
				        markersize=pP2D.ms / (j + 1),
				        label=lab
				)
				plt.legend(fancybox=True, framealpha=0.5, loc='center left', bbox_to_anchor=(1.0, 0.5))

		pP2D.r = pPP.ext_r
		pP2D.v = pPP.ext_v
		pP2D.b = pPP.ext_b
		color_gradient(len(pP2D.plotsExtPath.keys()), pP2D)
		if (key in pP2D.plotsExt):
			for key_path in pP2D.plotsExtPath:
				path = pP2D.plotsExtPath[key_path]
				# Selecting marker and color
				m = batch_markers.pop()
				c = pP2D.colors.pop()
				# Getting data
				data = {}
				# Getting the extension
				ext = path.split(".")[-1]
				if ext == "py":
					execfile(path)
				elif ext == "csv":
					with open(path) as f:
						dialect = csv.Sniffer().sniff(f.read(1024), delimiters=';')
						f.seek(0)
						reader = csv.DictReader(f, dialect=dialect)
						data = {k: [] for k in reader.fieldnames}
						for row in reader:
							for k in reader.fieldnames:
								data[k].append(float(row[k].replace(",", ".")))
				elif ext == "txt":
					with open(path) as f:
						lines = [line.rstrip().split() for line in f.readlines()]
						keys = lines[0]
						data = {key: [] for key in keys}
						for line in lines[1:]:
							for i in range(len(keys)):
								data[keys[i]].append(float(line[i]))
				for i in range(len(pP2D.plotsExt[key][0])):
					x = pP2D.plotsExt[key][0][i]
					y = pP2D.plotsExt[key][1][i]
					me = int(max(1.0, pP2D.me * len(data[x])))
					if pP2D.mec == None:
						mec = c
					else:
						mec = pP2D.mec
					plt.plot(
					        data[x],
					        data[y],
					        color=c,
					        marker=m,
					        markevery=me,
					        markerfacecolor=c,
					        markeredgewidth=pP2D.mew,
					        markeredgecolor=mec,
					        markersize=pP2D.ms,
					        label=key_path,
					        linestyle='None'
					)
					plt.legend(fancybox=True, framealpha=0.5, loc='center left', bbox_to_anchor=(1.0, 0.5))
		# Adding grid
		plt.grid()
		# Saving figure
		if pPP.save_figs:
			plt.savefig(pPP.save_fig_dir + pPP.name_case + "_" + pPP.name_param + "_" + key + ".pdf", bbox_inches="tight")

	for key in pP2D.loglogs:
		plt.figure()
		plt.xlabel(pPP.plots_names[pP2D.loglogs[key][0][0]])
		plt.ylabel(pPP.plots_names[pP2D.loglogs[key][1][0]])
		if pP2D.alims[key][0]:
			plt.xlim(pP2D.alims[key][0][0], pP2D.alims[key][0][1])
		if pP2D.alims[key][1]:
			plt.ylim(pP2D.alims[key][1][0], pP2D.alims[key][1][1])
		# Plotting
		batch_markers = pP2D.markers[:]
		pP2D.r = pPP.n_r
		pP2D.v = pPP.n_v
		pP2D.b = pPP.n_b
		color_gradient(len(params), pP2D)
		# Normal plots
		for i in range(len(params)):
			p = params[i]
			v = params_val[i]
			m = batch_markers.pop()
			c = pP2D.colors.pop()
			for j in range(len(pP2D.loglogs[key][0])):
				x = pP2D.loglogs[key][0][j]
				y = pP2D.loglogs[key][1][j]
				me = int(max(1.0, pP2D.me * len(v[x])))
				lab = ""
				if pPP.name_param != "":
					lab += r'$' + pPP.name_param + "=" + str(p) + "$"
				if len(pP2D.loglogs[key]) > 2:
					if lab != "":
						lab += " "
					lab += pP2D.loglogs[key][2][j]

				if pP2D.mec == None:
					mec = c
				else:
					mec = pP2D.mec
				plt.loglog(
				        v[x],
				        v[y],
				        color=c,
				        marker=m,
				        markevery=me,
				        markeredgewidth=pP2D.mew,
				        markeredgecolor=mec,
				        markerfacecolor=c,
				        markersize=pP2D.ms / (j + 1),
				        label=lab
				)
				plt.legend(fancybox=True, framealpha=0.5, loc='center left', bbox_to_anchor=(1.0, 0.5))

		pP2D.r = pPP.add_r
		pP2D.v = pPP.add_v
		pP2D.b = pPP.add_b
		color_gradient(len(pP2D.plot_adds.keys()), pP2D)
		if (key in pP2D.plot_adds):
			# Additional data plots
			m = batch_markers.pop()
			c = pP2D.colors.pop()
			for j in range(len(pP2D.plot_adds[key][0])):
				x = pP2D.plot_adds[key][0][j]
				y = pP2D.plot_adds[key][1][j]
				me = int(max(1.0, pP2D.me * len(v[x])))
				lab = ""
				if len(pP2D.plot_adds[key]) > 2:
					if lab != "":
						lab += " "
					lab += pP2D.plot_adds[key][2][j]

				if pP2D.mec == None:
					mec = c
				else:
					mec = pP2D.mec
				plt.loglog(
				        badata[x],
				        badata[y],
				        color=c,
				        marker=m,
				        markevery=me,
				        markeredgewidth=pP2D.mew,
				        markeredgecolor=mec,
				        markerfacecolor=c,
				        markersize=pP2D.ms / (j + 1),
				        label=lab
				)
				plt.legend(fancybox=True, framealpha=0.5, loc='center left', bbox_to_anchor=(1.0, 0.5))

		pP2D.r = pPP.ext_r
		pP2D.v = pPP.ext_v
		pP2D.b = pPP.ext_b
		color_gradient(len(pP2D.plotsExtPath.keys()), pP2D)
		if (key in pP2D.loglogsExt):
			for key_path in pP2D.plotsExtPath:
				path = pP2D.plotsExtPath[key_path]
				# Selecting marker and color
				m = batch_markers.pop()
				c = pP2D.colors.pop()
				# Getting data
				data = {}
				# Getting the extension
				ext = path.split(".")[-1]
				if ext == "py":
					execfile(path)
				elif ext == "csv":
					with open(path) as f:
						dialect = csv.Sniffer().sniff(f.read(1024), delimiters=';')
						f.seek(0)
						reader = csv.DictReader(f, dialect=dialect)
						data = {key: [] for key in reader.fieldnames}
						for row in reader:
							for k in reader.fieldnames:
								data[k].append(float(row[k].replace(",", ".")))
				elif ext == "txt":
					with open(path) as f:
						lines = [line.rstrip().split() for line in f.readlines()]
						keys = lines[0]
						data = {key: [] for key in keys}
						for line in lines[1:]:
							for i in range(len(keys)):
								data[keys[i]].append(float(line[i]))
				for i in range(len(pP2D.loglogsExt[key][0])):
					x = pP2D.loglogsExt[key][0][i]
					y = pP2D.loglogsExt[key][1][i]
					me = int(max(1.0, pP2D.me * len(data[x])))
					if pP2D.mec == None:
						mec = c
					else:
						mec = pP2D.mec
					plt.loglog(
					        data[x],
					        data[y],
					        color=c,
					        marker=m,
					        markevery=me,
					        markerfacecolor=c,
					        markeredgewidth=pP2D.mew,
					        markeredgecolor=mec,
					        markersize=pP2D.ms,
					        label=key_path,
					        linestyle='None'
					)
					plt.legend(fancybox=True, framealpha=0.5, loc='center left', bbox_to_anchor=(1.0, 0.5))
		# Adding grid
		plt.grid()
		# Saving figure
		if pPP.save_figs:
			plt.savefig(pPP.save_fig_dir + pPP.name_case + "-loglog" + "_" + pPP.name_param + "_" + key + ".pdf", bbox_inches="tight")

if pP1D.plot_enable:
	## Adding legends
	for key in axs:
		axs[key].legend(fancybox=True, framealpha=0.5, loc='center left', bbox_to_anchor=(1.0, 0.5))
	for key in axsT:
		axsT[key].legend(fancybox=True, framealpha=0.5, loc='center left', bbox_to_anchor=(1.0, 0.5))

	## Adding grids
	for key in axs:
		axs[key].grid()
	for key in axsT:
		axsT[key].grid()

	## Adding Patchs
	for key in pP1D.patchs:
		p = pP1D.patchs[key]
		rect = plt.Rectangle(eval(p['pos']), eval(p['w']), eval(p['h']), facecolor='w', edgecolor='k', hatch='/', alpha=0.3)
		axs[key].add_patch(rect)

	### Converting xlabel with radian writing
	#axs["rotx"].set_xticklabels([r"$" + format(r/np.pi, ".2g")+ r"\pi$" for r in axs["rotx"].get_xticks()])
	#axs["roty"].set_xticklabels([r"$" + format(r/np.pi, ".2g")+ r"\pi$" for r in axs["roty"].get_xticks()])
	#axs["rotz"].set_xticklabels([r"$" + format(r/np.pi, ".2g")+ r"\pi$" for r in axs["rotz"].get_xticks()])

	### Saving figures
	if pPP.save_figs:
		print(bigSep + "Saving figures.")
		for key in figs:
			figs[key].savefig(pPP.save_fig_dir + pPP.name_case + "_" + pPP.name_param + "_" + key + ".pdf", bbox_inches="tight")
		for key in figsT:
			figsT[key].savefig(pPP.save_fig_dir + pPP.name_case + "_" + pPP.name_param + "_" + key + "T.pdf", bbox_inches="tight")
		for key in figsO:
			for i in range(4):
				axsO[key].view_init(i * 90.0 / 3.0, -90.0)
				figsO[key].savefig(pPP.save_fig_dir + pPP.name_case + "_" + pPP.name_param + "_" + key + str(i) + ".pdf", bbox_inches="tight")
			axsO[key].view_init(45, -45.0)
			figsO[key].savefig(pPP.save_fig_dir + pPP.name_case + "_" + pPP.name_param + "_" + key + "4" + ".pdf", bbox_inches="tight")

### Showing figures
if pPP.show_figs:
	plt.show()
