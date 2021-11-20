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
	        "mean_phi": r"$\bar{\phi}$",
	        "mean_vx": r"$\bar{U^p_x}^* = \frac{\bar{U^p_x}}{\sqrt{(\rho_s/\rho_f - 1) g " + d_ad_name + "}}$",
	        "vx": r"${U^p_x}^* = \frac{U^p_x}{\sqrt{(\rho_s/\rho_f - 1) g" + d_ad_name + "}}$",
	        "mean_vfx": r"$\bar{U^f_x}^* = \frac{\bar{U^f_x}}{\sqrt{(\rho_s/\rho_f - 1) g " + d_ad_name + "}}$",
	        "vfx": r"${U^f_x}^* = \frac{U^p_x}{\sqrt{(\rho_s/\rho_f - 1) g " + d_ad_name + "}}$",
	        "mean_qsx": r"$\bar{q_s}^*(z)$",
	        "qs": r"${Q_s}^*$",
	        "qf": r"${Q_f}^*$",
	        "shields": r"$\theta$",
	        "sh": r"$\theta$",
	        "shmu": r"$\frac{\theta}{\mu}$",
	        "shsqmu": r"$\frac{\theta}{\sqrt{\mu}}$",
	        "z": r"$z^* =  \frac{z}{" + d_ad_name + "}$",
	        "time": r"$t$ (s)",
	        "mean_z_phi": r"$\phi_{max}$",
	        "var_z_phi": r"$\sigma_\phi$",
	        "lt": r"$\lambda$",
	        "mean_lm": r"$lm^* = \frac{lm}{d_{vs}}$",
	        "mean_re": r"$R_{xy}^f$",
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
	post_process = [
	        {
	                # Adimensionalisation
	                "z": "[i * pF.dz/d_ad for i in range(pN.n_z)]",
	                #"vx":"[adim(l, sqrt(-pM.g[2] * d_ad)) for l in data['vx']]", # For dry cases
	                "vx": "[adim(l, sqrt((pP.rho/pF.rho - 1.0) * -pM.g[2] * d_ad)) for l in data['vx']]",  # For wet cases
	                "vfx": "[adim(l, sqrt((pP.rho/pF.rho - 1.0) * -pM.g[2] * d_ad)) for l in data['vfx']]",  # Only for wet cases
	                "lm": "[adim(l, d_ad) for l in data['lm']]",  # Only for wet cases
	                "shields": "[s for s in data['shields']]",  # Only for wet cases
	                "shmu": "[s/mudAB[(pP.A,pP.B)] for s in data['shields']]",  # Only for wet cases
	                "shsqmu": "[s/sqrt(mudAB[(pP.A,pP.B)]) for s in data['shields']]",  # Only for wet cases
	        },
	        {
	                # Averaging.
	                "mean_phi": "average_profile(data['phi'], data['time'])",
	                "mean_vx": "ponderate_average_profile(data['phi'], data['vx'], data['time'])",
	                "mean_vfx": "average_profile(data['vfx'], data['time'])",  # Only for wet cases
	                "mean_lm": "average_profile(data['lm'], data['time'])",  # Only for wet cases
	                "mean_re": "average_profile(data['ReS'], data['time'])",  # Only for wet cases
	        },
	        {
	                # Flows
	                "mean_qsx":
	                        "[data['mean_phi'][i] * data['mean_vx'][i] for i in range(len(data['mean_phi']))]",
	                "qs":
	                        "[integration(data['phi'][i], data['vx'][i], pF.dz/d_ad) for i in range(len(data['time']))]",
	                "qf":
	                        "[integration([1.0 - p for p in data['phi'][i]], data['vfx'][i], pF.dz/d_ad) for i in range(len(data['time']))]",  # Only for wet cases
	                "mean_z_phi":
	                        "[np.mean(data['phi'][i][int(" + z_av_min + "/pF.dz):int(" + z_av_max + "/pF.dz)]) for i in range(len(data['time']))]",
	                "var_z_phi":
	                        "[sqrt(np.var(data['phi'][i][int(" + z_av_min + "/pF.dz):int(" + z_av_max + "/pF.dz)])) for i in range(len(data['time']))]",
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
	        "lm": [[], []],  # Only for wet cases
	        "re": [[], []],  # Only for wet cases
	        "qsx": [[], []],
	        "phi": [[], []],
	        "vx": [[], []],
	        "vfx": [[], []],  # Only for wet cases
	        "qs": [[], []],
	        "sh": [[], []],  # Only for wet cases
	        "mean_z_phi": [[], []],
	}
	# Normal plots
	plots = {
	        "lm": [["mean_lm"], ["z"]],  # Only for wet cases
	        "re": [["mean_re"], ["z"]],  # Only for wet cases
	        "qsx": [["mean_qsx"], ["z"]],
	        "phi": [["mean_phi"], ["z"]],
	        "vx": [["mean_vx"], ["z"]],
	        "vfx": [["mean_vfx"], ["z"]],  # Only for wet cases
	        "qs": [["time"], ["qs"]],
	        "sh": [["time"], ["shields"]],  # Only for wet cases
	        "mean_z_phi": [["time"], ["mean_z_phi"]],
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
	alimsO = {}
	orientations = {}
	# Plot a patch
	patchs = {
	        "phi": {
	                "pos": "(-500, " + z_av_min + "/d_ad)",
	                "w": "1000",
	                "h": "(" + z_av_max + "-" + z_av_min + ")/d_ad"
	        },
	}


class pP2D:
	plot_enable = False
	# Measuring 2D data in the 1D data
	measures = {
	        "qs": "average(data['qs'], data['time'])",
	        "lt": "computeTransportLayerThickness(data['mean_qsx'], pF.dz)",
	        "qf": "average(data['qf'], data['time'])",  # Only for wet cases
	        "sh": "average(data['shields'], data['time'])",  # Only for wet cases
	        # Phi
	        "phi_max": "average(data['mean_z_phi'], data['time'])",
	        "A": "average(data['A'], data['time'])",
	}
	post_process = [
	        {
	                "A_fit": "np.linspace(1.0, 3.0, 100)",
	        },
	        {
	                "fit":
	                        "np.polynomial.polynomial.polyval(bdata[bval]['A_fit'], np.polynomial.polynomial.polyfit(bdata[bval]['A'], bdata[bval]['phi_max'], 1))"
	        },
	]
	add = [
	        {
	                "theta3/2": "[s**(3.0/2.0) for s in badata['sh']]"  # Only for wet cases
	        },
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
	alims = {
	        "qs(sh)": [[], []],  # Only for wet cases
	        "lt(sh)": [[], []],  # Only for wet cases
	}
	plots = {
	        "lt(sh)": [["sh"], ["lt"]],  # Only for wet cases
	        "qs(sh)": [["sh"], ["qs"]],  # Only for wet cases
	        "qs(shmu)": [["shmu"], ["qs"]],  # Only for wet cases
	        "qs(shsqmu)": [["shsqmu"], ["qs"]],  # Only for wet cases
	        "phi_max(A)": [["A", "A_fit"], ["phi_max", "fit"], ["simulations", "fit"]],
	}
	loglogs = {
	        "qs(sh)": [["sh"], ["qs"]],  # Only for wet cases
	}
	plot_adds = {
	        "qs(sh)": [["sh"], ["theta3/2"], [r"$\theta^{\frac{3}{2}}$"]],  # Only for wet cases
	}
	plotsExtPath = {}
	plotsExt = {}
	loglogsExt = {}
