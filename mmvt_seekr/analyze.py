#!/usr/bin/env python

#===============================================
# MODULE DOCSTRING
#=================================================

"""
MMVT SEEKR module with main analysis functions

calculates rates, convergence, minimum simulation estimates from a milestone model 
"""

import random, os
from subprocess import check_output
from pprint import pprint
from math import exp, log
import numpy as np
from scipy import linalg as la
from scipy.stats import gamma as gamma
from itertools import islice
import xml.etree.ElementTree as ET

R_GAS_CONSTANT = 0.0019872 # in kcal/mol*K
k_boltz = 1.3806488e-23




def analyze_kinetics(model, bound_indices, max_steps =[None], verbose=False,):
	'''main function to perform all kinetics analyses.
	Given a Model() object with its statistics filled out, it will return an estimate of the kinetic
	value, given all or a subset of its statistics.

	Parameters
	-----------
	model : object
		Model() object containing all milestone information and transition statistics
	bound_indices : list
		list of intiger indices of all milestones to be treated as bound states
	max_steps : list
		list of maximum step numbers that are used ti limit amount of data used in kinetics calculations
		can be set by the user and is also used for convergence analysis
	verbose: bool
		enable verbose printing- provides additional details such as values of transition rate matrix


	Returns
	---------
	p_equil: numpy array
		array of equilibrium probability of each voronoi cell 

	N : numpy array
		array of number of transitions between all milestones

	R : numpy array
		array of transition rates for all milestones

	T : numpy array
		array of mean first passage times from bound state to all other states
	
	T_tot : numpy array
		array of total simulation time spent in each cell

	Q : numpy array
		Transition rate matrix

	N_conv : list
		list of transition count arrays (N) for each convergence interval

	R_conv : list
		list of transition rate arrays (R) for each convergence interval

	k_cell : numpy array
		array of transition rates between Voronoi cells

	'''
	counts = {}; times = {}; total_counts = {}; total_cell_counts = {}; total_times = {}; avg_times = {}; trans = {}; total_cell_times = {}; T_a = {}
	end_indeces = [];
	N = {}
	if verbose: print('max_steps', max_steps) 
	for site in model.sites:
		for anchor in site.anchors:
			if anchor.md == True and anchor.directory:
				if verbose: print('Anchor', anchor.fullname)

				#unpack the list of max steps for each anchor
				if max_steps[0] == None:
					anchor_max_steps = None
				else:
					anchor_max_steps = max_steps[int(anchor.index)]


				this_counts, this_total_counts, this_total_times, this_avg_times = anchor._get_md_transition_statistics(model.md_time_factor, 
						anchor_max_steps)
				this_cell_counts, this_cell_time = anchor._get_md_vt_collisions(model.md_time_factor, anchor_max_steps)
				
				total_counts = _add_dictionaries(total_counts, this_total_counts)
				if verbose: print('counts',  this_counts)
				total_cell_counts = _add_dictionaries(total_cell_counts, this_cell_counts)
				total_cell_times[int(anchor.index)] = this_cell_time
				if verbose: print('times', this_total_times)
				if verbose: print('cell times', total_cell_times)
				#total_times = add_dictionaries(total_times, this_total_times)
				#total_cell_times = add_dictionaries(total_cell_times, this_cell_times)
				for src_key in list(this_counts.keys()):
					if src_key in list(counts.keys()):
						counts[src_key] = _add_dictionaries(counts[src_key], this_counts[src_key])
					else:
						counts[src_key] = this_counts[src_key]
				#print "len(transitions)", len(milestone.transitions)
			for milestone in anchor.milestones:
				if milestone.end == "true": # then its an end milestone, and there will be no transitions out of it
					end_indeces.append(int(milestone.id))
			for src_key in list(this_total_times.keys()):
				if src_key in list(times.keys()):
					times[src_key] = _add_dictionaries(times[src_key], this_total_times[src_key])
				else:
					times[src_key] = this_total_times[src_key]



## Calculate Voronoi cell equilibrium probability ##
	k_cell = np.zeros((len(total_cell_times),len(total_cell_times)))  
	k_mod = np.zeros((len(total_cell_times),len(total_cell_times))) 
	



	for cell in list(total_cell_counts.keys()):
		for new_cell in list(total_cell_counts[cell].keys()):
			if new_cell == -1: continue #skip transitions to bound state milestone
			elif new_cell not in end_indeces: #hard code for testing
				#print cell, ", ",  new_cell
				#print total_cell_counts[cell][new_cell]
				#print total_cell_times[cell]
				k_cell[cell][new_cell] = (float(total_cell_counts[cell][new_cell])/float(total_cell_times[cell]))
				#print k_cell[cell][new_cell]
		 
	#print "k_cell"; pprint(k_cell)
 
## Create the steady state flux matrix##

	for i in range(len(k_cell[1])):
		for j in range(len(k_cell[1])):
			if i == j:
				k_mod[i][j] = -(np.sum(k_cell[i]))
			else:
				k_mod[i][j] = k_cell[j][i]

## Substitute redundant equation with normalization condition

	k_mod[-1][:] = 1
	#print "k_mod"; pprint(k_mod)

	p_0= np.zeros((len(total_cell_times)), dtype="float")  
	p_0[-1] = 1.0

	#pprint(p_0)

	#print "k_cell:", np.shape(k_cell)
	#print "p_0:", np.shape(p_0)
	#k_cell_trans = np.transpose(k_cell)




## Calculate the equilibrium probabilities for each voronoi cell
	p_equil = np.linalg.solve(k_mod, p_0)
	if verbose: pprint(p_equil)
	#print np.shape(p_equil)
	

	p_equil_ref = p_equil[-1]
	#print p_equil[-1]
	#print p_equil_ref
	#print range(len(p_equil))
	delta_G = np.zeros(len(p_equil))
	for i in range(len(p_equil)):
		#print i
		#print p_equil[i] 
		delta_G[i] = -model.temperature * R_GAS_CONSTANT * log(p_equil[i] / p_equil_ref)
		#print delta_G[i]


	#print("Delta G: "); pprint(delta_G)
## Using the V cell equilibrium probabilities, calculate the rate matrix, Q
	if verbose: print("counts: ", counts)
	if verbose: print("times: ", times)

	T_a = np.zeros(len(p_equil))
	for cell in total_cell_times:
		T_a[cell] = p_equil[cell]/total_cell_times[cell]
	T_tot = 1/np.sum(T_a)


	N = np.zeros((len(p_equil)+1,len(p_equil)+1))
	N_conv = np.zeros((len(p_equil)+1,len(p_equil)+1))
	for anchor in list(counts.keys()):
		for src in list(counts[anchor].keys()):
			for dest in list(counts[anchor][src].keys()):
				#print p_equil[int(anchor)]
				#print counts[anchor][src][dest]
				#print total_cell_times[int(anchor)]
				N[src][dest] = p_equil[int(anchor)] * float(counts[anchor][src][dest])/ total_cell_times[int(anchor)]
				if max_steps[0] != None: 
					if  total_cell_times[int(anchor)] >= max_steps[int(anchor)] * model.md_time_factor: 
						N_conv[src][dest] = float(counts[anchor][src][dest])/ total_cell_times[int(anchor)]
					else:
						N_conv[src][dest] = np.nan      

	if verbose: print("N:", N)

	R = np.zeros(len(p_equil)+1)
	R_conv = np.zeros((len(p_equil)+1,len(p_equil)+1))
	for anchor in list(times.keys()):
		for src in list(times[anchor].keys()): 
			R[src] += (p_equil[int(anchor)] * times[anchor][src]/ total_cell_times[int(anchor)])
			if max_steps[0] != None:
				if total_cell_times[int(anchor)] >= max_steps[int(anchor)] * model.md_time_factor:
					R_conv[int(anchor)][src] = times[anchor][src]/ total_cell_times[int(anchor)] 
				else:
					R_conv[int(anchor)][src] = np.nan

	if verbose: print("R:", R)


	Q = np.zeros((len(p_equil)+1,len(p_equil)+1))
	for i in range(len(N[0])):
		for j in range(len(N[0])):
			Q[i][j] = N[i][j]/R[i]
			

	for i in range(len(Q[0])):
		Q[i][i] = -np.sum(Q[i])



	if verbose: print("")
	if verbose: print(Q)

# Calculate MFPT 
	
	T= _calc_MFPT_vec(Q)


	total_sim_time = 0
	for i in list(total_cell_times.keys()):
		if verbose: print(i, total_cell_times[i]*1e9, "ns")
		total_sim_time += total_cell_times[i]

	if verbose: print("Total simulation time: " ,  total_sim_time*1e9, "ns")


	return p_equil, N, R, T, T_tot, Q, N_conv, R_conv, k_cell, 


def _get_index_dict(trans_dict):
	'''
	Currently not used
	given a transition (or count) dictionary, will return a transition (or count) matrix
	'''
	index_dict = {}
	trans_dict_keys = trans_dict.keys()
	trans_dict_keys = [key.replace('inf', 'inf_0') for key in trans_dict_keys]
	n = len(trans_dict_keys)
	#print trans_dict_keys
	i = 0
	trans_dict_keys = sorted(trans_dict_keys, key=lambda keystring:keystring.split('_')[0]) # first sort by site
	trans_dict_keys = sorted(trans_dict_keys, key=lambda keystring:int(keystring.split('_')[1])) # then sort by milestone
	#print trans_dict_keys
	#count_dict_keys = sorted(count_dict.keys(), key=lambda keystring:keystring )#int(keystring.split('_')[1])) # then sort by milestone
	#count_dict_keys.sort()
	for key in trans_dict_keys:
		if key == 'inf_0': key = 'inf'
		index_dict[i] = key
		i += 1

def calc_kon_from_bd(model, bound_indices, Q):
	'''
	Calculate the on rate for the model by incorporating BD simulation data with 
	MMVT md simulation data. k_on is calculated using the Northrup, Allsion, McCammon method

	Parameters
	----------
	model: object
		the milestoning model object containing all transition information
	bound_indices : list
		list of integer indices for all milestones to be treated as bound states
	Q : numpy array
		the transition rate matrix from the MMVT simulations

	Returns
	--------
	k_on : float
		the calculated association rate constant

	'''
	# b-surface milestone
	inf_index = -1
	bd_time = 0.0
	n= Q.shape[0]

	#print("Q", Q)
	#get preliminary transition matrix, K
	K_prelim = _rate_mat_to_prob_mat(Q) #converts rate matrix Q back to a probability matrix and incubation time vector
	K_prelim_mod = K_prelim
	K_prelim_mod[-1,:] = 0.0
	#print("K prelim", K_prelim)
	K = np.zeros((n+1,n+1))
	K[:-1,:-1] = K_prelim_mod

	#print("K", K)

	#extract BD milestone statistics for infinity state
	bd_counts, bd_total_counts, bd_total_times, bd_avg_times = model.bd_milestone._get_bd_transition_statistics(results_filename ="results.xml", bd_time=bd_time)

	src_key = model.bd_milestone.index
	#print("src_key", src_key)
	#print(bd_counts)

	
	bd_counts[src_key][inf_index]= bd_counts[src_key].pop('inf')
	#print(bd_counts)

	#add BD escape probability to K
	for key in bd_counts[src_key].keys():
		K[src_key][int(key)]=float(bd_counts[src_key][key])/float(bd_total_counts[src_key])
	#add infinity state to K
	K[inf_index][:] = 0.0
	K[inf_index][inf_index] = 1.0
	#modify K for bound/sink states
	for key in bound_indices:
		K[int(key),:] = 0.0
		K[int(key)][int(key)] =1 

	#print("K  mod", K.shape, K)
	K_trans = np.transpose(K)
	K = K_trans
	#print("k trans", K)

	#extract BD statistics from B surface calculation
	b_surface_counts, b_surface_total_counts, b_surface_total_times, b_surface_avg_times = model.b_surface._get_bd_transition_statistics(results_filename="results.xml", bd_time=bd_time)
	#print(b_surface_counts)
	#src_key = b_surface_counts[0]
	src_key = model.b_surface.index
	b_surface_trans = {src_key:{}}
	b_surface_counts[src_key][inf_index]= b_surface_counts[src_key].pop('inf')
	#print("b surf counts", b_surface_counts)

	for dest_key in b_surface_counts[src_key].keys():
		b_surface_trans[src_key][dest_key] = float(b_surface_counts[src_key][dest_key]) / float(b_surface_total_counts[src_key])

	#b_surface_trans = _process_trans_for_bd(b_surface_trans, inf_index)
	q0 = _trans_dict_to_q0_vector(b_surface_trans,K)

	#print("q0", q0)
	beta = _get_beta_from_K_q0(K, q0, bound_indices)
	#print("beta", beta)
	k_b = _run_compute_rate_constant(results_filename=os.path.join("b_surface", "results.xml"), browndye_bin_dir="")
	#print( "k(b):", k_b)
	k_on = k_b * beta

	return k_on


def _process_trans_for_bd(trans_dict, inf_index):
  '''
  Not currently used
  given a transition (or count) dictionary, will return a transition (or count) matrix
  '''
  new_trans_dict = {}
  src_key = trans_dict.keys()[0] #get the source milestone
  new_src_key = src_key.split('_')[2]
  for key, val in trans_dict[src_key].items():
    print(key)
    key.replace('inf', inf_index)
    #new_key = key.split('_')[1]
    #print(new_key)
    #new_trans_dict[new_src_key] = trans_dict[key] 


  print(trans_dict)


  return trans_dict


def _trans_dict_to_q0_vector(trans_dict, K,):
	'''given a transition matrix, will return a vector of starting fluxes based on b-surface stats.'''
	n = K.shape[0]
	#print("b trans_dict", trans_dict)
	src_key = list(trans_dict)[0] # the key that refers to the b-surface trans dict
	q0_vector = np.zeros((n,1))
	for key in trans_dict[src_key].keys():
			q0_vector[int(key),0] = trans_dict[src_key][key]
	return q0_vector

def _calc_MFPT_vec(Q):
	Q_hat = Q[:-1,:-1]

	#if verbose: print Q_hat

	I = np.zeros(len(Q_hat[0]))
		#I[:] = np.sqrt(len(Q_hat[0]))
	I[:] = 1
	
	T= la.solve(Q_hat,-I)

		#MFPT = T[0]
		#if verbose: print "MFPT =", T, "fs"

		#k_off = 1e15/MFPT

		#print "T", T  

	return T

def _rate_mat_to_prob_mat(Q):
	n = Q.shape[0]
	K= np.matrix(np.zeros((n,n)))
	sum_vector = np.zeros(n)
	
	for i in range(n):
		sum_vector[i]= -1* Q[i,i]

	#print("sum vec", sum_vector)

	for i in range(n):
		for j in range(n):
			if i==j: continue
			K[i,j] = Q[i,j] / sum_vector[i] 	

	return K

def _run_compute_rate_constant(results_filename, browndye_bin_dir=""):
	'runs the Browndye program compute_rate_constant to find the value k(b)'
	process_trajectories = os.path.join(browndye_bin_dir, "compute_rate_constant")
	cmd = "%s < %s" % (process_trajectories, results_filename)
	output_string = check_output(cmd, shell=True) # run the Browndye process command
	root = ET.fromstring(output_string) # read the XML string
	rate = root[0] # the first <rate> tag, because it doesn't really matter which one we choose
	rate_constant = rate.find('rate-constant')
	k_on_tag = rate_constant.find('mean')
	reaction_probability = rate.find('reaction-probability')
	beta_tag = reaction_probability.find('mean')
	assert (k_on_tag != None) and (beta_tag != None), "Alert: Unable to find rate constant <mean> tags after running compute_rate_constant on file %s" % results_filename
	k = float(k_on_tag.text)
	beta = float(beta_tag.text)
	k_b = k / beta # the flux to the b-surface
	return k_b

def _get_beta_from_K_q0(K, q0, bound_indices):
	'given a transition matrix K and starting vector q_0, returns a beta value for k-ons'
	K_inf = np.matrix(K) ** 99999999
	#print("K inf", K_inf)
	#n,m = K_inf.shape
	q_inf = np.dot(K_inf, q0)
	#print("q inf", q_inf)
	#print "q_inf:", q_inf
	beta = 0.0
	for bound_index in bound_indices:
		beta += q_inf[bound_index, 0]
	return beta

def monte_carlo_milestoning_error(model, bound_indices, Q0, N_pre, R_pre, p_equil, T_tot, num = 1000, skip = 100, stride =1,  verbose= False):
	'''Calculates an error estimate by sampling a distribution of rate matrices assumming 
	a poisson (gamma) distribution with parameters Nij and Ri using Markov chain Monte Carlo
		
	Enforces detailed Balance-- using a modified version of Algorithm 4 form Noe 2008 for rate matrices.--  
	Distribution is:  p(Q|N) = p(Q)p(N|Q)/p(N) = p(Q) PI(q_ij**N_ij * exp(-q_ij * Ri))


	Parameters
	----------
	model : object
		milestoning model object containing all transition and milestone information
	bound_indices : list
		list of indices of all milestones to be treated as bound states
	Q0 : numpy array
		Initial rate matrix (from MLE estimate) on which perturbations will be made to obtain a distribution
	N_pre : numpy array
		inital count matrix used to generate rate matrix perturbations from a gamma distribution
	R_pre : numpy arrray
		initial transition time matrix used to generate rate matrix perturbations from a gamma distribution
	p_equil : numpy array
		equilibrium probabilities oof each Voronoi cell (generated from MLE estimate)
	T_tot : numpy array
		total simulation time in each Voronoi cell
	num : int
		number of rate matrix (Q) samples to be generated
	skip : int
		number of inital rate matrix samples to skip for "burn in"
	stride : int
		frequency at which rate matrix samples are recorded- larger frequency reduces correlation between samples
	verbose : bool
		allow additional verbosity/printing

	Returns
	-------
	k_off_list : list
		list of calculated k off values calculated from rate matrix samples
	running_avg : list
		average k off calculated at each sampling interval (used for plotting convergence)
	running_std : list
		standard deviation of k off calculated at each interval (used for plotting convergence)
	k_on_list : list
		list of calculated k on values calculated from rate matrix samples
	k_on_avg : list
		average k on calculated at each sampling interval (used for plotting convergence)
	k_on_std : list
		standard deviation of k on calculated at each interval (used for plotting convergence)

	'''
	m = N_pre.shape[0] #get size of count matrix
	Q = Q0
	Q_mats = []
	N = []
	R = []
	k_off_list = []
	running_avg = []
	running_std = []
	k_on_list = []
	k_on_avg = []
	k_on_std =[]
	
	
	N = N_pre
	R = R_pre

	if verbose: print("Q", Q.shape)
	if verbose: print(Q)
	P = np.zeros((m,m))
	Q_test = np.zeros((m,m))
	tau = np.zeros(m)


	if verbose: print("N", N)
	if verbose: print("R", R)
	
	Qnew = Q
	if verbose: print("collecting ", num, " MCMC samples from ", num*(stride)+ skip, " total moves")  
	for counter in range(num*(stride)+ skip):
		if verbose: print("MCMC stepnum: ", counter)
		Qnew = np.zeros((m,m)) #np.matrix(np.copy(T))
		for i in range(m): # rows
			for j in range(m): # columns
				Qnew[i,j] = Q[i,j]


		for i in range(m): # rows
			for j in range(m): # columns
				if i == j: continue
				if Qnew[i,j] == 0.0: continue
				if Qnew[j,j] == 0.0: continue


				Q_gamma = 0
				delta = Qnew[i,j]
				while ((delta) >= (Qnew[i,j])):# or (Qnew[i,j] - Q_gamma) >= abs(Qnew[i,i])):
					Q_gamma = gamma.rvs(a=N[i,j], scale = 1/R[i],)
					
					delta =  Qnew[i,j] - Q_gamma
	

				log_p_Q_old = N[i,j] * log(Qnew[i,j])  - Qnew[i,j] * R[i] #+ -Qnew[i,i] * R[i]

				log_p_Q_new = N[i,j] * log(Qnew[i,j] - delta) - (Qnew[i,j] - delta) * R[i] #+ -(Qnew[i,i] + delta) * R[i]

				if verbose: print("log P(Q_new)", log_p_Q_new)
				if verbose: print("log P(Q_old)", log_p_Q_old)

				r2 = random.random()  
				p_acc =  log_p_Q_new - log_p_Q_old
				#if verbose: print("p_acc", p_acc, "r", log(r2))
					
				if log(r2) <= p_acc: #log(r) can be directly compared to log-likeliehood acceptance, p_acc
					#if verbose: print("performing non-reversible element shift...")
						

					Qnew[i,i] = (Qnew[i,i]) + delta
					Qnew[i,j] = Qnew[i,j] - delta


					if verbose: print(Qnew)

		if counter > skip and counter % stride == 0:
				T_err = _calc_MFPT_vec(Qnew)
				k_off_list.append(1/T_err[0])
				running_avg.append(np.average(k_off_list))
				running_std.append(np.std(k_off_list))
				k_on_list.append(calc_kon_from_bd(model,bound_indices, Qnew))
				k_on_avg.append(np.average(k_on_list))
				k_on_std.append(np.std(k_on_list))	
 
		Q = Qnew
	if verbose: print("final MCMC matrix", Q)
	return k_off_list, running_avg, running_std, k_on_list, k_on_avg, k_on_std  

def check_milestone_convergence(model, bound_indices, conv_stride, skip, max_steps, verbose=False,):
	'''
	Calculates the key MMVT quantities N, R, and Q as a function of simulation time
	to estimate which milestones have been sufficiently sampled. 

	Quantities are pulled from the data at step intervals determined by the conv_stride value 
	with the option to skip steps from the beginning of the data with the skip variable

	Parameters
	-----------
	model : object
		milestoning model object containing all milestone and transition information
	bound indices : list
		list of integer indices of milestones to be treated as bound states
	conv_stride : int
		frequency at which a sample is pulled for N and R from the simulation data
	skip : int
		number of steps to be skipped from the beginning of each simulation
	max_steps : int
		total number of simulation steps to be used in the calculation
	verbose: bool
		allow additional printing

	Returns
	-------
	N_conv: list
		list of transition count matrix N for each convergence interval
	R_conv : list
		list of transition time matrix R for each convergence interval 
	k_cell_conv : list
		list of cell transition rates for each convergence interval
	p_equil_conv : list
		list of cell equilibrium probabilities for each convergence interval
	k_off_conv : list
		list of calculated off rate at each convergence interval
	k_on_conv : list
		list og calculated on rate at each convergence interval
	conv_intervals : list
		list of step numbers used for each convergence sample

	'''
	n_anchors = 0
	for site in model.sites:
		n_anchors += site.num_anchors   
	conv_intervals = np.arange(conv_stride, max_steps, conv_stride)
	conv_intervals = conv_intervals + skip
	max_step_list = np.zeros(n_anchors)
	N_conv = np.zeros((15,15,len(conv_intervals)))
	R_conv = np.zeros((15,15,len(conv_intervals)))
	k_off_conv = np.zeros(len(conv_intervals))
	k_on_conv = np.zeros(len(conv_intervals))
	k_cell_conv = np.zeros((15,15,len(conv_intervals)))
	p_equil_conv = np.zeros((15,len(conv_intervals)))

	for interval_index in range(len(conv_intervals)):
		max_step_list[:] = conv_intervals[interval_index]
		p_equil, N, R, T, T_tot, Q, n_conv, r_conv, k_cell = analyze_kinetics(model, bound_indices, max_steps=max_step_list, verbose=verbose,)

		MFPT = T[0]
		k_off = 1/MFPT
		k_on = calc_kon_from_bd(model, bound_indices, Q)

		for index, x in np.ndenumerate(n_conv):
				N_conv[index[0]][index[1]][interval_index]=x
		for index2,y in np.ndenumerate(r_conv):
			R_conv[index2[0]][index2[1]][interval_index]= y 
		for index3,z in np.ndenumerate(k_cell):
			k_cell_conv[index3[0]][index3[1]][interval_index]= z
		for index4,j in np.ndenumerate(p_equil):
			p_equil_conv[index4[0]][interval_index]= j   
		k_off_conv[interval_index]=k_off
		k_on_conv[interval_index] = k_on

	return N_conv, R_conv, k_cell_conv, p_equil_conv, k_off_conv, k_on_conv, conv_intervals, 

def _add_dictionaries(dict1, dict2):
	'''
	adds the values numerically within each dictionary
	NOTE: dict1 is updated and returned BY REFERENCE
	'''
	new_dict = dict1
	for key in list(dict2.keys()):
		if key in list(dict1.keys()):
			dict1[key] += dict2[key]
		else:
			dict1[key] = dict2[key]

	return dict1
def _make_windows(seq, n):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result
        
def _calc_window_rmsd(conv_values):
	#print(conv_values)
	average = np.mean(conv_values)
	#print(average)
	test =conv_values - average
	#print(test)
	new_test = np.delete(test, 0)
	#print(new_test)
	RMSD = np.sqrt(np.sum(new_test)**2/(len(conv_values) -1))
	#print(RMSD)
	return RMSD, average


def _find_conv_min(conv_values, conv_intervals, window_size, cutoff, conv_windows):

	conv_times = np.zeros((conv_values.shape[0], conv_values.shape[1]))
	#conv_times[:] = np.nan
	for i in range(conv_values.shape[0]):
		for j in range(conv_values.shape[1]):
			if np.sum(conv_values[i][j][:]) != 0:
				#print("I, J", i,j)
				conv_times[i][j]= np.nan
				rmsd_list = []
				windows = _make_windows(conv_values[i][j][:], window_size)
				index = 0
				conv_counter = 0
				for w in windows:
					RMSD, window_average = _calc_window_rmsd(w)
					#N_rmsd_list.append(RMSD)
					#print(index)
					#print(index)
					#print(conv_counter, conv_windows)
					if RMSD <= (cutoff * window_average):
						#consider anchor converged
						#print("RMSD less than cutoff", index, RMSD, window_average)
						#print("Convergence meeasures:", conv_counter,"of",  conv_windows)
						conv_counter += 1
						#print("counter", conv_counter, "index", index)
						if conv_counter == conv_windows:
							#print("Entry %i , %i converged at index: %i" %(i,j,index))
							max_int = index + window_size
							min_time = conv_intervals[max_int]
							#print(min_time)
							conv_times[i][j]= min_time
							break
					else: conv_counter = 0
					index += 1
					
				if np.isnan(conv_times[i][j]): print("Entry %i, %i did not meet minimum convergence criteria of %i windows" %(i,j,conv_windows))
	return conv_times


def calc_RMSD_conv(model, N_conv, R_conv, conv_intervals, window_size, cutoff, conv_windows):
	''' Estimate the convergence of sampling for each milestone using a sliding window RMSD cutoff
	Milestones are considered converged when the values of BOTH N and R remain below the cutoff threshold for 
	a specified number of windows.

	The cutoff is a percentage of the magnutude of the corresponding value (e.g. less than 5% of the magnutude of N)

	Parameters
	-----------
	model : object
		milestoning model object containing all milestone information and transition statistics
	N_conv: list
		list of transition count matrix N for each convergence interval
	R_conv : list
		list of transition time matrix R for each convergence interval 
	conv_intervals : list
		list of stepnumbers for which convergence samples were extracted and calculated
	window size : int
		number of convergence samples used in a window RMSD calculation
	cutoff : float
		RMSD must be less than the cutoff times the magnitude of the quantity (N or R) to be considered converged
	conv_windows : int
		number of consecutive windows the RMSD must remain below the cutoff to be considered converged

	Return
	-------
	min_anchor_times :list
		list of times (in stepnumbers) for each voronoi cell where the convergence criteria are satisfied
	'''

	min_anchor_times = np.zeros((len(N_conv[0])))
	#print(min_anchor_times)
	print("Calculating N/T convergence")
	n_conv_times = _find_conv_min(N_conv, conv_intervals, window_size, cutoff, conv_windows)
	print("Calculating R/T convergence")
	r_conv_times = _find_conv_min(R_conv, conv_intervals, window_size, cutoff, conv_windows)

	for site in model.sites:
		for anchor in site.anchors:
			n_steps = 0
			if np.isnan(r_conv_times[int(anchor.index)][:]).any():
				continue
			r_steps = max(r_conv_times[int(anchor.index)][:])
			for milestone in anchor.milestones:
				if np.isnan(n_conv_times[int(milestone.id)][:]).any():
					continue
				else:
					if max(n_conv_times[int(milestone.id)][:]) > n_steps:               
						n_steps = max(n_conv_times[int(milestone.id)][:])
			min_anchor_times[int(anchor.index)] = max(n_steps, r_steps) 
			print("anchor", anchor.index, min_anchor_times[int(anchor.index)])


	return min_anchor_times
