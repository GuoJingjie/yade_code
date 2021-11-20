# 2019 © Remi Monthiller <remi.monthiller@gmail.com>
#!/usr/bin/python

#-------------------#
# Import
#-------------------#
# Plot modules
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
# Readdata modules
import os
import sys
import pickle
import csv

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
        'size': 18
}
plt.rc('font', **font)

#-------------------#
# Plot parameters
#-------------------#
# Basics
enable_show_figs = True
enable_save_figs = False
directory_save_figs = "./post_proc/"
# Average time steps number
av_dt_nb = 50

#-------------------#
# Data
#-------------------#

# Read Simulation Data
data = {}
case_directory = "./case-wet/"
with open(case_directory + 'data.pkl', 'rb') as f:
	data = pickle.load(f)
	print("Avalable data : " + str(data.keys()))

# Getting vx
if len(data['vx']) < av_dt_nb:
	vx = np.array(data['vx'])
	# Take only last time step :
	vx_average = np.array(data['vx'][-1])
	print("WARNING : No averaging but took only last time step.")
else:
	vx = np.array(data['vx'])
	# Average over the last time steps
	vx_average = np.average(np.array(data['vx']), axis=0)

# Getting phi
if len(data['phi']) < av_dt_nb:
	vx = np.array(data['phi'])
	# Take only last time step :
	vx_average = np.array(data['phi'][-1])
	print("WARNING : No averaging but took only last time step.")
else:
	vx = np.array(data['phi'])
	# Average over the last time steps
	vx_average = np.average(np.array(data['phi']), axis=0)

# Getting z
z = np.array(data['z'][-1])

# Other Data
m_A = np.array([1.0, 1.4, 1.8, 2.2, 2.6, 3.0])
m_alpha = np.array([0.375, 0.445, 0.545, 0.595, 0.595, 0.595])
m_mu = np.tan(m_alpha)
alpha_A = poly.Polynomial(poly.polyfit(m_A, m_alpha, 3))
mu_A = poly.Polynomial(poly.polyfit(m_A, m_mu, 3))

#-------------------#
# Plot
#-------------------#

fig = plt.figure()
ax = plt.gca()
plt.xlabel(r"$\overline{U_x^p}$")
plt.ylabel(r"$z$")
x = vx_average
y = z
m = "."  # Marker
c = (0.5, 0.0, 0.0)  # Color (red, green, blue)
me = int(max(1.0, 0.05 * len(x)))  # Markevery
ax.errorbar(x, y, yerr=0, color=c, marker=m, markevery=me, markerfacecolor=c, markeredgewidth=0, markersize=10)

#-------------------#
# Save
#-------------------#

if enable_save_figs:
	fig.savefig(pPP.save_fig_dir + "simple_post_proc" + ".pdf", bbox_inches="tight")

#-------------------#
# Show
#-------------------#

if enable_show_figs:
	plt.show()
