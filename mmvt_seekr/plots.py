"""
make_model.py
Simulation Enabled Estimation of Kinetic Rates (SEEKR) is a tool that facilitates the preparation, running and analysis of multiscale MD/BD/Milestoning  simulations for the calculation of protein-ligand binding kinetics.

extracts transition information from the simulation output files and creates the milestoning model

Parameters
		----------
		milestone_filename : string Required 
				name of the XML file containing all information regarding simulation directories, parameters, etc.

		Returns
		-------
		model : class 
				contains all required information for milestoning analysis

		max_steps : int 
				total number of s
"""
import matplotlib.pyplot as plt
import pickle 
import numpy as np
from cycler import cycler


def plot_n_conv(conv_values, conv_intervals):
	fig, ax = plt.subplots()
	new_colors = [plt.get_cmap('tab20')(1.* i/(_get_colormap(conv_values))) for i in range(_get_colormap(conv_values))]
	plt.rc('axes', prop_cycle=(cycler('color', new_colors)))
	#ax = fig.add_subplot(1,1,1,)

	for i in range(conv_values.shape[0]):
		for j in range(conv_values.shape[1]):
				if np.sum(conv_values[i][j][:]) != 0:
					label_string = 'Src: '+str(i) +',' + 'Dest: '+str(j)
					ax.plot(np.multiply(conv_intervals,2e-6), conv_values[i][j][:], label = label_string ,
						linestyle='-', marker="o", markersize = 1)
	plt.xlabel('time (ns)')
	plt.ylabel('N/T')
	#plt.legend(loc = 'right')
	#box = fig.get_position()
	#plt.set_position([box.x0,box.y0, box.width * 0.8, box.height])
	plt.legend(loc ='center left', bbox_to_anchor=(1, 0.5), ncol = 2)  
	plt.grid(b=True,axis = 'y', which = 'both')
	return fig, ax

def plot_r_conv(conv_values, conv_intervals):
	fig, ax = plt.subplots()
	new_colors = [plt.get_cmap('tab20')(1.* i/(_get_colormap(conv_values))) for i in range(_get_colormap(conv_values))]
	plt.rc('axes', prop_cycle=(cycler('color', new_colors)))
	#ax = fig.add_subplot(1,1,1,)

	for i in range(conv_values.shape[0]):
		for j in range(conv_values.shape[1]):
				if np.sum(conv_values[i][j][:]) != 0:
					label_string = 'anchor ' +str(i) + ',' + 'Milestone '+str(j)
					ax.plot(np.multiply(conv_intervals,2e-6), conv_values[i][j][:], label = label_string ,
						linestyle='-', marker="o", markersize = 1)
	plt.xlabel('time (ns)')
	plt.ylabel('R/T')
	#plt.legend(loc = 'right')
	#box = fig.get_position()
	#plt.set_position([box.x0,box.y0, box.width * 0.8, box.height])
	plt.legend(loc ='center left', bbox_to_anchor=(1, 0.5), ncol = 2)  
	plt.grid(b=True,axis = 'y', which = 'both')
	return fig, ax

def _get_colormap(conv_values):
	cmap_length = 0
	for i in range(conv_values.shape[0]):
		for j in range(conv_values.shape[1]):
				if np.sum(conv_values[i][j][:]) != 0:
					cmap_length +=1
	return cmap_length

def plot_p_equil(conv_values, conv_intervals):
	fig, ax = plt.subplots()
	new_colors = [plt.get_cmap('tab20')(1.* i/conv_values.shape[0]) for i in range(conv_values.shape[0])]
	plt.rc('axes', prop_cycle=(cycler('color', new_colors)))
#f, (ax, ax2, ax3) = plt.subplots(3, 1, sharex=True)
	for i in range(conv_values.shape[0]):
		label_string = 'anchor ' +str(i)
		ax.plot(np.multiply(conv_intervals,2e-6), p_equil_conv[i][:], 
						 label = label_string ,linestyle='-', marker="o", markersize = 1)
	plt.xlabel('time (ns)')
	plt.ylabel('p equil')
	plt.legend(loc ='center left', bbox_to_anchor=(1, 0.5), ncol = 2)  
	plt.grid(b=True,axis = 'y', which = 'both')
	return fig, ax

def plot_k_conv(conv_values, conv_intervals):
	fig, ax = plt.subplots()
	ax.plot(np.multiply(conv_intervals,2e-6), conv_values, linestyle='-', marker="o", markersize = 1)
	plt.ylabel('k off (s^-1)')
	plt.xlabel('time (ns)')
	#plt.legend(loc ='center left', bbox_to_anchor=(1, 0.5))
	return fig, ax

def MCMC_conv(running_avg, running_std):
	fig = plt.figure()
	ax = fig.add_subplot(2,1,1,)
	ax.plot(running_avg)
	ax.set_ylabel('Average off rate (1/s)')
	ax.set_xlabel('MCMC Samples')
	ax2 = fig.add_subplot(2,1,2,)
	ax2.plot(running_std)
	ax2.set_ylabel('off rate st. dev. (1/s)')
	ax2.set_xlabel('MCMC Samples')

	plt.savefig('MCMC_conv.png', format ='png', dpi=300 )
	pickle.dump(fig, open('MCMC.fig.pickle', 'wb'))
	plt.show()
	return fig, ax1, ax2 