# 2019 © Remi Monthiller <remi.monthiller@gmail.com>
#-------------------#
# Setting latex
#-------------------#
plt.rc('text', usetex=True)
#-------------------#
# Setting Font
#-------------------#
font = {
        'family': 'normal',
        #        'weight' : 'bold',
        'size': 24
}
plt.rc('font', **font)

# Dictionaries of data
phiAB = {
        (1.0, 1.0): 0.64,
        (1.5, 1.0): 0.0,
        (2.0, 1.0): 0.6,
        (2.5, 1.0): 0.0,
        (3.0, 1.0): 0.56,
}
mut018AB = {
        (1.0, 1.0): tan(0.475),
        (1.5, 1.0): tan(0.545),
        (2.0, 1.0): tan(0.615),
        (2.5, 1.0): tan(0.615),
        (3.0, 1.0): tan(0.615),
}
mut1AB = {}
musAB = { # Forte incertitude spheres
        (1.0, 1.0): tan(0.45),
        (1.5, 1.0): tan(0.60),
        (2.0, 1.0): tan(0.71),
        (2.5, 1.0): tan(0.71),
        (3.0, 1.0): tan(0.71),
        (2.0, 2.0): tan(0.815),
        (3.0, 2.0): tan(0.755),
        (3.0, 3.0): tan(0.785),
        (4.0, 2.0): tan(0.795),
        (5.0, 1.0): tan(0.715),
        (5.0, 2.0): tan(0.775)
}
mudAB = {
        (1.0, 1.0): tan(0.375),
        (1.05, 1.0): tan(0.385),
        (1.1, 1.0): tan(0.385),
        (1.2, 1.0): tan(0.415),
        (1.3, 1.0): tan(0.425),
        (1.4, 1.0): tan(0.435),
        (1.5, 1.0): tan(0.465),
        (2.0, 1.0): tan(0.575),
        (2.5, 1.0): tan(0.575),
        (3.0, 1.0): tan(0.575),
        (2.0, 2.0): tan(0.645),
        (3.0, 2.0): tan(0.635),
        (3.0, 3.0): tan(0.635),
        (4.0, 2.0): tan(0.655),
        (5.0, 1.0): tan(0.59),
        (5.0, 2.0): tan(0.675)
}


# Basic plot parameters
class pPP:
	#-------------------#
	# Plot Names
	#-------------------#
	show_figs = False
	#-------------------#
	# Saving figures
	#-------------------#
	save_fig_dir = "./post_proc/"
	save_figs = True
	if save_figs and not os.path.isdir(save_fig_dir):
		os.mkdir(save_fig_dir)
	#-------------------#
	# Adimensionalisation
	#-------------------#
	d_ad = "pP.dvs"
	d_ad_name = "d_{vs}"
	#-------------------#
	# Mean operation
	#-------------------#
	mean_begin_time = 200.0
	mean_end_time = 1000.0
	#-------------------#
	# Plot visuals
	#-------------------#
	# Color gradient
	r = [0.5, 0.5]
	v = [0.0, 0.5]
	b = [0.0, 0.25]
	# Color gradient for normal plots
	n_r = [0.5, 0.0]
	n_v = [0.5, 0.5]
	n_b = [0.5, 0.0]
	# Color gradient for external data
	ext_r = [0.0, 0.5]
	ext_v = [0.5, 0.0]
	ext_b = [0.5, 0.0]
	# Color gradient for additional data
	add_r = [0.0, 0.5]
	add_v = [0.5, 0.0]
	add_b = [0.5, 0.0]
	# Plot markers
	markers = [
	        "$\mathbf{H}$", "$\mathbf{G}$", "$\mathbf{F}$", "$\mathbf{E}$", "$\mathbf{D}$", "$\mathbf{C}$", "$\mathbf{B}$", "$\mathbf{A}$", "d", "*",
	        r"$\mathbf{\gamma}$", r"$\mathbf{\lambda}$", r"$\mathbf{\delta}$"
	]
	# Marker density
	me = 0.05
	# Marker edge width
	mew = 1.5
	# Marker edge color
	mec = None  #(1.0, 1.0, 1.0)
	# Marker size
	ms = 10.0
	#-------------------#
	# Plot Names
	#-------------------#
	plots_names = {
	        "z": r"$z^* =  \frac{z}{" + d_ad_name + "}$",
	        "time": r"$t$ (s)",
	        "dirocc": r"Probability",
	        "mtheta": r"$\theta (rad)$",
	        "mphi": r"$\phi (rad)$",
	}
	# Plot params
	name_case = "wet"
	param = "pP.A"
	name_param = "A"


# 1D plot parameters
class pP1D:
	plot_enable = True
	#-------------------#
	# Variable Utils
	#-------------------#
	# Parameters for spatial averaging
	z_av_min = "pP.dvs*5"
	z_av_max = "pP.dvs*9"
	#-------------------#
	# Post Processing
	#-------------------#
	rot_i = "-1"
	post_process = [
	        {
	                "dirocc": "data['dirs'][" + rot_i + "][3]",
	                "dirx": "data['dirs'][" + rot_i + "][0]",
	                "diry": "data['dirs'][" + rot_i + "][1]",
	                "dirz": "data['dirs'][" + rot_i + "][2]",
	                "mdirc": "data['mdirs'][" + rot_i + "][3]",
	                "mdirx": "data['mdirs'][" + rot_i + "][0]",
	                "mdiry": "data['mdirs'][" + rot_i + "][1]",
	                "mdirz": "data['mdirs'][" + rot_i + "][2]",
	                "z": "[z/d_ad for z in data['ori'][" + rot_i + "][0]]",
	                "mtheta": "data['ori'][" + rot_i + "][1]",
	                "vtheta": "[t/2.0 for t in data['ori'][" + rot_i + "][2]]",
	                "mphi": "data['ori'][" + rot_i + "][3]",
	                "vphi": "[p/2.0 for p in data['ori'][" + rot_i + "][4]]",
	        }
	]
	#-------------------#
	# Plot Visuals
	#-------------------#
	r = pPP.r
	v = pPP.v
	b = pPP.b
	colors = []
	markers = pPP.markers[:]
	me = pPP.me
	mew = pPP.mew
	mec = pPP.mec
	ms = pPP.ms
	#-------------------#
	# Plots
	#-------------------#
	# Axis limits
	alims = {
	        "mtheta": [[0.0, pi], []],
	        "mphi": [[-pi / 2.0, pi / 2.0], []],
	}
	# Normal plots
	plots = {
	        "mtheta": [["mtheta"], ["z"], ["vtheta", ""]],
	        "mphi": [["mphi"], ["z"], ["vphi", ""]],
	}
	# Attempt to plot profiles as function of the time
	plotsT = {}
	# Plot external data path
	plotsExtPath = {
	        #"Experimental data: 6,\nFrey et al. 2014":"./exp-data/Frey2014_EXP6.py"
	}
	# Plot external data
	plotsExt = {
	        #"vx":[["vx"], ["z"]],
	        #"qsx":[["qsx"], ["z"]],
	        #"phi":[["phi"], ["z"]],
	}
	# Plot external data
	alimsO = {
	        "dirs": [[-0.5, 1.5], [-1.0, 1.0], [-1.0, 1.0]],
	        "mdirs": [[-0.5, 1.5], [-1.0, 1.0], [-1.0, 1.0]],
	}
	orientations = {
	        "dirs": [[["dirx", "diry", "dirz"]], ["dirocc"]],
	        "mdirs": [[["mdirx", "mdiry", "mdirz"]], ["mdirc"]],
	}
	# Plot a patch
	patchs = {}


class pP2D:
	plot_enable = False
	# Measuring 2D data in the 1D data
	measures = {}
	post_process = [
	        {},
	        {},
	]
	add = [
	        {},
	]
	#-------------------#
	# Plot visuals
	#-------------------#
	r = pPP.r
	v = pPP.v
	b = pPP.b
	colors = []
	markers = pPP.markers[:]
	me = pPP.me
	mew = pPP.mew
	ms = pPP.ms
	#-------------------#
	# Plots
	#-------------------#
	alims = {}
	plots = {}
	loglogs = {}
	plot_adds = {}
	plotsExtPath = {}
	plotsExt = {}
	loglogsExt = {}
