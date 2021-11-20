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
	        },
	        {
	                # Averaging.
	                "mean_phi": "average_profile(data['phi'], data['time'])",
	                "mean_vx": "ponderate_average_profile(data['phi'], data['vx'], data['time'])",
	        },
	        {
	                # Flows
	                "mean_qsx": "[data['mean_phi'][i] * data['mean_vx'][i] for i in range(len(data['mean_phi']))]",
	                "qs": "[integration(data['phi'][i], data['vx'][i], pF.dz/d_ad) for i in range(len(data['time']))]",
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
	        "qsx": [[], [4, 9]],
	        "phi": [[], [4, 9]],
	        "vx": [[], [4, 9]],
	}
	# Normal plots
	plots = {
	        "qsx": [["mean_qsx"], ["z"]],
	        "phi": [["mean_phi"], ["z"]],
	        "vx": [["mean_vx"], ["z"]],
	        "qs": [["time"], ["qs"]],
	        "sh": [["time"], ["shields"]],
	}
	# Attempt to plot profiles as function of the time
	plotsT = {}
	# Plot external data path
	plotsExtPath = {"Experimental data: 20,\nFrey et al. 2014": "./data-exp/Frey2014_EXP20.py"}
	# Plot external data
	plotsExt = {
	        "vx": [["vx"], ["z"]],
	        "qsx": [["qsx"], ["z"]],
	        "phi": [["phi"], ["z"]],
	}
	# Plot external data
	alimsO = {}
	orientations = {}
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
