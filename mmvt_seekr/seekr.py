"""
seekr.py
control script for generating all input files for SEEKR calculation

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


##Parse SEEKR input file###

import argparse
import os, re, shutil, sys
from itertools import islice

inputs = { #contains all DEFAULT parameters for seekr input file
	###Genaral Parameters###
	'empty_rootdir':False,
	'watermodel':'',
	'master_temperature':298,
	'equil_steps':1000000,
	'equil_write_freq':10000,
	'prod_steps':10000000,
	'prod_write_freq':100000,

}

def _parse_seekr_input(inp_filename):
	new_inputs = {}
	listopen = False # this gives us the ability to read parameters across several lines
	ourlist = []
	for line in open(inp_filename,'r'):
		splitline = line.strip().split() + ['#']

		if not listopen: splitlist = []
		for i in range(len(splitline)):
			if splitline[i][0] == '[':
				splitline[i] = splitline[i][1:].strip()
				listopen = True
				if splitline[i] and splitline[i][-1] == ']':
					listopen = False
					break
			if splitline[i] and splitline[i][-1] == ']':
				listopen = False
				splitlist.append(splitline[i][:-1])
				break

			if splitline[i].startswith('#'): # the remainder of this line is a comment
				break
			splitlist.append(splitline[i].strip())

		if not listopen:
			if len(splitlist) < 2: 
				continue
			word1 = splitlist[0].lower()
			word2 = " ".join(splitlist[1:])
			if word1[0] == '#': # this line is a comment
				continue
			new_inputs[word1] = word2 # populate the input dictionary with this value
	return new_inputs

def _get_inputs(args,):
# combine the default input with the user-defined input
#inputscript = sys.argv[1] # the input file for SEEKR
	if args['input_filename']:
		input_filename = args['input_filename']
		new_inp = _parse_seekr_input(input_filename)
		#print(new_inp)
		inputs.update(new_inp)
		#print(inputs)
	else:
		unittest.main()
		exit()
	return inputs

def _get_sys_params(inp):
	sys_params = {#list of parameters needed to generate filetree
		'project_name':inp['project_name'],
		'rootdir':inp['rootdir'],
		'system_pdb_filename':inp['system_pdb_filename'], #Protein-Ligand complex
		'system_parm_filename':inp['system_parm_filename'], #for AMBER FF
		'system_psf_filename':inp['system_psf_filename'],   # for CHARMM FF
		'lig_pqr_filename':inp['lig_pqr_filename'], #for BD simulations
		'rec_dry_pqr_filename':inp['rec_dry_pqr_filename'], #for BD simulations
		'empty_rootdir':_boolean(inp['empty_rootdir']),
		}
	return sys_params

def _get_md_settings(inp):
	md_settings= { # settings for the md module
	##General shared MD settings
		'ff':inp['ff'], # the forcefield to use
		'watermodel':inp['watermodel'],
		'master_temperature':inp['master_temperature'],
	##Restrained equilibration settings
		'equil_settings':{
			'namd_settings':{
				'numsteps':inp['equil_steps'],
				'write_freq':inp['equil_write_freq']}
				},
	##MMVT Production specific settings
		'prod_settings':{
			'namd_settings':{
				'numsteps':inp['prod_steps'],
				'write_freq':inp['prod_write_freq']}
				},
		}

	return md_settings

def _boolean(arg):
	if str(arg).lower() in ("false","", "0", "0.0", "[]", "()", "{}"):
		return False
	else:
		return True

def _parse_milestone_inputs(inp):
	milestones = []
	for key in sorted(inp.keys()):
		if re.match("milestone[0-9]+", key):
			milestone_dict = {'key': key}
			for param in inp[key].split(','):
				split_param = param.strip().split()
				if len(split_param) < 2: continue
				milestone_key = split_param[0]
				milestone_value = ' '.join(split_param[1:]) 
				milestone_dict[milestone_key] = milestone_value
			milestones.append(milestone_dict)
	return milestones

def _generate_milestone_lists(milestones, anchor_list):
	for milestone in milestones:
		milestone_pairs = []
		milestone_vals_string = milestone['milestone_values']
		raw_milestone_vals = milestone_vals_string.split()
		#print("milestone values",  raw_milestone_vals)
		pairs = _make_milestone_pairs(raw_milestone_vals)
		index = 0
		for p in pairs:
			milestone_pairs.append(p)
			if index not in anchor_list:
				#print(index, p)
				anchor_list.append(p)
			else: 
				anchor_list[index].append(p)
			index +=1
		#print(milestone_pairs)
		milestone['milestone_pair_list'] = milestone_pairs
		
	return milestones, anchor_list

def _make_milestone_pairs(seq, n=2):
	"Returns a sliding window (of width n) over data from the iterable"
	"   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
	it = iter(seq)
	result = tuple(islice(it, n))
	if len(result) == n:
		yield result
	for elem in it:
		result = result[1:] + (elem,)
		yield result

def _distance_milestones(milestone):

	return
	

def _generate_filetree(inp, sys_params):
	print("Creating file tree at location: %s" % inp['rootdir'])
	if sys_params['empty_rootdir'] and os.path.isdir(inp['rootdir']): # then we're in test mode, delete anything that is in the rootdir
		print("'empty_rootdir' set to True. Deleting all subdirectories/files in rootdir:", inp['rootdir'])
		shutil.rmtree(inp['rootdir']) # remove all the files within the grid subdirectory
	if not os.path.isdir(inp['rootdir']): # if the rootdir doesn't exist yet
		os.mkdir(inp['rootdir']) # then create it
	return


def prepare_seekr():
	anchor_list = []
	parser = argparse.ArgumentParser(description="top level SEEKR program used to prepare and generate all input files necessary for a SEEKR calculation")
	parser.add_argument('input_filename', metavar='INPUT_FILENAME', type = str, help="name of SEEKR input file")
	args = parser.parse_args() #parse command line arguments
	args = vars(args) #converts to dictionary

	_get_inputs(args,)
	sys_params = _get_sys_params(inputs) #parse general system parameters (structures, forcefield files, directory names from input file)
	md_settings = _get_md_settings(inputs) #parse parameters for MD simulations from input file
	milestones = _parse_milestone_inputs(inputs) #parse milestone CV parameters and milestone values from input file
	milestones, anchor_list = _generate_milestone_lists(milestones, anchor_list) #generate upper/lower bounds of a SINGLE CV for each anchor/Voronoi cell
	#print(anchor_list)
	_generate_filetree(inputs, sys_params) #creates/clears top level ditectory
	#filetree_settings_all=dict(filetree_settings.items() + sys_params.items())
	#filetree.main()




	#for milestone in milestones:
	#	print('Milestone pairs', milestone['milestone_pair_list'])




if __name__ == "__main__": prepare_seekr()


