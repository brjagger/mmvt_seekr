"""
analyze.py
Simulation Enabled Estimation of Kinetic Rates (SEEKR) is a tool that facilitates the preparation, running and analysis of multiscale MD/BD/Milestoning  simulations for the calculation of protein-ligand binding kinetics.

Performs kinetic analysis, including calculation of rate matrix, MFPT (on and off rates), and Milestone free energy

Parameters
		----------
		calc_type: string default= "off" 
				type of calculation, options "on" or "off"
		model: object Required
			the SEEKR milestoning model
		bound_dict: dictionary Required
			dictionary of bound states in the milestoning model



		Returns
		-------
		model : class 
				contains all required information for milestoning analysis
"""
import random
from pprint import pprint
from math import exp, log
import numpy as np
from scipy import linalg as la
from scipy.stats import gamma as gamma
from itertools import islice


R_GAS_CONSTANT = 0.0019872 # in kcal/mol*K
k_boltz = 1.3806488e-23




def analyze_kinetics(calc_type, model, bound_indices, max_steps =[None], verbose=False,):
	'''main function to perform all kinetics analyses.
	Given a Model() object with its statistics filled out, it will return an estimate of the kinetic
	value, given all or a subset of its statistics.
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


				this_counts, this_total_counts, this_total_times, this_avg_times = anchor.get_md_transition_statistics(model.md_time_factor, 
						anchor_max_steps)
				this_cell_counts, this_cell_time = anchor.get_md_vt_collisions(model.md_time_factor, anchor_max_steps)
				
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


#      if anchor.bd == True and anchor.directory:
#        this_counts, this_total_counts, this_total_times, this_avg_times = milestone.get_bd_transition_statistics(bd_time=bd_time)
#        print 'TIME', this_avg_times
#        total_counts = add_dictionaries(total_counts, this_total_counts)
#        total_times = add_dictionaries(total_times, this_total_times)
#        for src_key in this_counts.keys():
#          if src_key in counts.keys():
#            counts[src_key] = add_dictionaries(counts[src_key], this_counts[src_key])
#          else:
#            counts[src_key] = this_counts[src_key]

 # for src_key in total_times.keys(): # construct the average incubation times
 #   R_cell[src_key] = total_times[src_key] / total_cell_times[src_key]

#  for src_key in counts.keys():
#    temp = {}
#    for dest_key in counts[src_key].keys(): # make the transition probability dictionary
			#temp[dest_key] = float(counts[src_key][dest_key]) / float(total_counts[src_key])
			#N_cell[src_key][dest_key]



## Calculate Voronoi cell equilibrium probability ##
	#for cell_src_key in total_cell_counts.keys():
	#total_cell_times = {0: 84668000, 1: 86494000, 2:86504000, 3:92662000, 4:43832000, 5:42430000, 6:43522000} #hard coded -- includes times from combined face sampling
	#print "total_cell_counts"; pprint(total_cell_counts)  
	#print "total times"; pprint(total_times)
	#print "cell times"; pprint(total_cell_times)
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
	
	T= calc_MFPT_vec(Q)


	total_sim_time = 0
	for i in list(total_cell_times.keys()):
		if verbose: print(i, total_cell_times[i]*1e9, "ns")
		total_sim_time += total_cell_times[i]

	print("Total simulation time: " ,  total_sim_time*1e9, "ns")


	return p_equil, N, R, T, T_tot, Q, N_conv, R_conv, k_cell, 


def get_index_dict(trans_dict):
	'''given a transition (or count) dictionary, will return a transition (or count) matrix'''
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
	# b-surface milestone
	n= Q.shape[0]

	#get preliminary transition matrix, K
	K, avg_t = rate_mat_to_prob_mat(Q) #converts rate matrix Q back to a probability matrix and incubation time vector
	
	#add infinity state to K
	K[n+1][n+1] = 1.0
	#modify K for bound/sink states
	for key in bound_indices:
		K[int(key),:] = 0.0
		K[int(key)][int(key)] =1 


	#extract BD statistics from B surface calculation
	b_surface_counts, b_surface_total_counts, b_surface_total_times, b_surface_avg_times = model.b_surface.get_bd_transition_statistics(
		results_filename="results.xml", bd_time=bd_time)
	src_key = b_surface_counts.keys()[0]
	b_surface_trans = {src_key:{}}


	for dest_key in b_surface_counts[src_key].keys():
		b_surface_trans[src_key][dest_key] = float(b_surface_counts[src_key][dest_key]) / float(b_surface_total_counts[src_key])

	q0 = trans_dict_to_q0_vector(b_surface_trans)

	beta = get_beta_from_K_q0(K, q0, bound_indices)

	k_b = run_compute_rate_constant(results_filename=os.path.join("b_surface", "results.xml"), browndye_bin_dir="")
	print( "k(b):", k_b)
	k_on = k_b * beta

	return k_on

def trans_dict_to_q0_vector(trans_dict, K):
	'''given a transition matrix, will return a vector of starting fluxes based on b-surface stats.'''
	n = K.shape[0]
	src_key = trans_dict.keys()[0] # the key that refers to the b-surface trans dict
	q0_vector = np.zeros((n,1))
	for i in range(n): # move down the vector element by element
		dest_key = i
		if dest_key in trans_dict[src_key].keys():
			q0_vector[i,0] = trans_dict[src_key][dest_key]
	return q0_vector

def calc_MFPT_vec(Q):
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


def rate_mat_to_prob_mat(Q, calc_type='on'):
	''' converts a rate matrix Q into probability matrix (kernel) K and an incubation time vector'''
	n = Q.shape[0]
	P = np.matrix(np.zeros((n,n)))
	K = np.matrix(np.zeros((n,n)))
	sum_vector = np.zeros(n)
	avg_t = np.zeros(n)
	for i in range(n): # first make the sum of all the rates
		for j in range(n):
			if j == i: continue
			sum_vector[j] += Q[i,j]


	for i in range(n):
		for j in range(n):
			if j == i: continue
			if sum_vector[j] == 0:
				K[i,j] = 0.0
				#avg_t[j] = 0.0
			else:
				K[i,j] = Q[i,j] / sum_vector[j]

		if sum_vector[i] != 0.0:
			avg_t[i] = 1.0 / sum_vector[i]
		else: # then the sum vector goes to zero, make the state go to itself
			if calc_type == "on":
				K[i,i] = 1.0
			elif calc_type == "off":
				K[i,i] = 0.0
			#avg_t[i] = something???

	return K, avg_t

def run_compute_rate_constant(results_filename, browndye_bin_dir=""):
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

def get_beta_from_K_q0(K, q0, bound_indices):
	'given a transition matrix K and starting vector q_0, returns a beta value for k-ons'
	K_inf = np.matrix(K) ** 99999999
	#n,m = K_inf.shape
	q_inf = np.dot(K_inf, q0)
	#print "q_inf:", q_inf
	beta = 0.0
	for bound_index in bound_indices:
		beta += q_inf[bound_index, 0]
	return beta

def monte_carlo_milestoning_error(Q0, N_pre, R_pre, p_equil, T_tot, num = 1000, skip = 100, stride =1,  verbose= False):
	'''Samples distribution of rate matrices assumming a poisson (gamma) distribution with parameters Nij and Ri using Markov chain Monte Carlo
		Enforces detailed Balance-- using a modified version of Algorithm 4 form Noe 2008 for rate matrices.--  
	Distribution is:  p(Q|N) = p(Q)p(N|Q)/p(N) = p(Q) PI(q_ij**N_ij * exp(-q_ij * Ri))
	'''
	m = N_pre.shape[0] #get size of count matrix
	Q = Q0
	Q_mats = []
	N = []
	R = []
	k_off_list = []
	running_avg = []
	running_std = []
	
	
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
				T_err = calc_MFPT_vec(Qnew)
				k_off_list.append(1/T_err[0])
				running_avg.append(np.average(k_off_list))
				running_std.append(np.std(k_off_list))
 
		Q = Qnew
	if verbose: print("final MCMC matrix", Q)
	return k_off_list, running_avg, running_std 

def check_milestone_convergence(model, bound_indices, conv_stride, skip, max_steps, calc_type, verbose=False,):
	'''

	'''
	n_anchors = 0
	for site in model.sites:
		n_anchors += site.num_anchors   
	conv_intervals = np.arange(conv_stride, max_steps, conv_stride)
	conv_intervals = conv_intervals + skip
	max_step_list = np.zeros(n_anchors)
	N_conv = np.zeros((15,15,len(conv_intervals)))
	R_conv = np.zeros((15,15,len(conv_intervals)))
	k_conv = np.zeros(len(conv_intervals))
	k_cell_conv = np.zeros((15,15,len(conv_intervals)))
	p_equil_conv = np.zeros((15,len(conv_intervals)))

	for interval_index in range(len(conv_intervals)):
		max_step_list[:] = conv_intervals[interval_index]
		p_equil, N, R, T, T_tot, Q, n_conv, r_conv, k_cell = analyze_kinetics(calc_type, model, bound_indices, max_steps=max_step_list, verbose=verbose,)

		MFPT = T[0]
		k_off = 1/MFPT

		for index, x in np.ndenumerate(n_conv):
				N_conv[index[0]][index[1]][interval_index]=x
		for index2,y in np.ndenumerate(r_conv):
			R_conv[index2[0]][index2[1]][interval_index]= y 
		for index3,z in np.ndenumerate(k_cell):
			k_cell_conv[index3[0]][index3[1]][interval_index]= z
		for index4,j in np.ndenumerate(p_equil):
			p_equil_conv[index4[0]][interval_index]= j   
		k_conv[interval_index]=k_off

	return N_conv, R_conv, k_cell_conv, p_equil_conv, k_conv, conv_intervals, 

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
	min_anchor_times = np.zeros((len(N_conv[0])))
	#print(min_anchor_times)
	print("Calculating N/T convergence")
	n_conv_times = _find_conv_min(N_conv, conv_intervals, window_size, cutoff, conv_windows)
	print("Calculating R/T convergence")
	r_conv_times = _find_conv_min(R_conv, conv_intervals, window_size, cutoff, conv_windows)
	
	#print(n_conv_times)
	#print(r_conv_times)



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