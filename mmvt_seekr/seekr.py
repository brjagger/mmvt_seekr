"""
seekr.py
main script for generating all input files for SEEKR calculation

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
import xml.etree.cElementTree as ET
from xml.dom.minidom import Document
from xml.dom import minidom
import pickle
import filetree, md, bd
import pdb2 as pdb

inputs = { #contains all DEFAULT parameters for seekr input file
	###Genaral Parameters###
	'empty_rootdir':False,
	'watermodel':'',
	'master_temperature':298,
	'equil_steps':1000000,
	'equil_write_freq':10000,
	'prod_steps':10000000,
	'prod_write_freq':100000,
	'cell_shape': 'box',
	'extendedsystem_filename':'',
	'system_bin_coordinates':'',
	'eval_stride' : 10,	
	'md_time_factor' : 2e-15,
	'bd_time_factor' : 1,
	'empty_pqrxml_path':"./empty.pqrxml",
	'browndye_bin_dir':'',
	'ion1rad':None,
	'ion2rad':None,
	'ion1conc':None,
	'ion2conc':None,
	'ions':[],
	'lpbe_npbe':'lpbe',
	'apbs_executable':'apbs',
	'inputgen_executable':"/software/repo/moleculardynamics/apbs/2018.2.1/share/apbs/tools/manip/inputgen.py",
	'inputgen_fadd':'100',
	'inputgen_gmemceil':'64000',
	'inputgen_resolution':'0.5',
	'inputgen_cfac':'4.0',
	'fhpd_numtraj':'1000',
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
		'system_rst_filename':inp['system_rst_filename'],
		'system_bin_coordinates':inp['system_bin_coordinates'],
		'extendedsystem':inp['extendedsystem_filename'],
		'lig_pqr_filename':inp['lig_pqr_filename'], #for BD simulations
		'rec_dry_pdb_filename':inp['rec_dry_pdb_filename'],
		'rec_dry_pqr_filename':inp['rec_dry_pqr_filename'], #for BD simulations
		#'rec_pdb_filename':inp['rec_pdb_filename'],
		'empty_rootdir':_boolean(inp['empty_rootdir']),
		'md_time_factor': inp['md_time_factor'],
		'bd_time_factor' : inp['bd_time_factor'],
		}
	if inp['ff'] == 'amber':
		sys_params.update({'system_params_filename':inp['system_parm_filename']}) #for AMBER FF	
	if inp['ff'] == 'charmm':
		sys_params.update({'system_params_filename':inp['system_psf_filename']})  # for CHARMM FF

	return sys_params

def _get_md_settings(inp, md_file_paths):
	md_settings= { # settings for the md module
	##General shared MD settings
		'ff':inp['ff'], # the forcefield to use
		'watermodel':inp['watermodel'],
		'master_temperature':inp['master_temperature'],
		'cell_shape':inp['cell_shape'],
	##Restrained equilibration settings
		'equil':_boolean(inp['equil']),
		'equil_settings':{				
				'ensemble': inp['equil_ensemble'],
				'namd_settings':{
					'write_freq':inp['equil_write_freq'],
					'numsteps':inp['equil_steps'],
					'extendedsystem':inp['extendedsystem_filename']
					}	
				},
	##MMVT Production specific settings
		'prod_settings':{
				'ensemble' : inp['prod_ensemble'],
				'namd_settings': {
					'numsteps':inp['prod_steps'],
					'write_freq':inp['prod_write_freq'],
					'eval_stride' : inp['eval_stride'],
					}
				},
		'md_file_paths':md_file_paths, # file paths to the MD directories in the anchor file		
		}

	return md_settings

def _get_bd_settings(inp, sys_params, struct):

	#parser = pdb2.Big_PDBParser()
	#bd_receptor_dry_pqr=parser.get_structure('bd_receptor_dry_pqr', sys_params['bd_rec_pqr_filename'], pqr=True)


	bd_settings={ # settings for the bd model
	'rec_struct':struct['receptor_dry_pqr'],
	'lig_struct': struct['ligand'],
	'temperature':inp['master_temperature'],
	'bd_centerx':inp['bd_centerx'],
	'bd_centery':inp['bd_centery'],
	'bd_centerz':inp['bd_centerz'],
	'threads':int(inp['bd_threads']),
	'fhpd_numtraj':inp['fhpd_numtraj'],
	# reaction sites information in milestone_settings below
	'browndye_bin_dir':inp['browndye_bin_dir'],
	#'bd_file_paths':bd_file_paths,
	'b_surface_path':os.path.join(sys_params['rootdir'], 'b_surface'),
	'bd_milestone_path':os.path.join(sys_params['rootdir'], 'bd_milestone'),
	#'bd_radius':inp['bd_radius'],
	#'prods_per_anchor':inp['bd_prods_per_anchor'],
	#'starting_surfaces':[], # x,y,z,radius
	#'ending_surfaces':[], # x,y,z,radius
    'n-trajectories':inp['n-trajectories'],
	'b_surface_ending_surfaces':[], # when the ligand is started from the b-surface, where it can possibly end
	'siteids':[],
	#'starting_conditions':'configs', # may be 'spheres' or 'configs'. 'spheres': the ligands are started at random points along the starting_surfaces and then classified into states. 'configs' take configurations from 'lig_configs'
	#'increments':[],
	#'milestones':milestone_list,
	#'milestone_pos_rot_list':milestone_pos_rot_list,
	'empty_pqrxml_path':inp['empty_pqrxml_path'],
	'apbs_settings':{
	  'apbs_executable':inp['apbs_executable'],
	  'ions':[],
	  #'ion1rad':inp['ion1rad'], # negative ion
	  #'ion2rad':inp['ion2rad'], # positive ion
	  #'ion1conc':inp['ion1conc'],
	  #'ion2conc':inp['ion2conc'],
	  'temp':inp['master_temperature'],
	},
	'inputgen_settings':{
	  'inputgen_executable':inp['inputgen_executable'],
	  'fadd':inp['inputgen_fadd'],
	  'gmemceil':inp['inputgen_gmemceil'],
	  'resolution':inp['inputgen_resolution'],
	  'cfac':inp['inputgen_cfac'],
	},
	}

	print("STRUCT COM", struct['rec_com'])

	for key in sorted(inp.keys()):
		if re.match("ion[0-9]+$", key): # then this is an ion
	  		ion_dict = {'key':key}
	  		for param in inp[key].split(','): # then pull every string in the list and make a dictionary
	  			split_param = param.strip().split()
	  			if len(split_param) < 2: continue
	  			ion_key = split_param[0]
	  			ion_value = ' '.join(split_param[1:])
	  			ion_dict[ion_key] = ion_value
  			bd_settings['apbs_settings']['ions'].append(ion_dict)
	  
	if inp['ion1conc'] and inp['ion1rad']: # if the old syntax is used
		assert len(bd_settings['apbs_settings']['ions'])==0, "Cannot use 'ion#' parameter at the same time that you are using 'ion1conc' and 'ion1rad' parameters."
		bd_settings['apbs_settings']['ions'].append({'concentration': inp['ion1conc'], 'charge': '-1.0', 'radius': inp['ion1rad'], 'name': 'ion1', 'key': 'ion1'})
		if inp['ion2conc'] and inp['ion2rad']:
	  		bd_settings['apbs_settings']['ions'].append({'concentration': inp['ion2conc'], 'charge': '1.0', 'radius': inp['ion2rad'], 'name': 'ion2', 'key': 'ion2'})
	
	print("bd_settings['apbs_settings']['ions']:", bd_settings['apbs_settings']['ions'])


	return bd_settings

def _boolean(arg):
	if str(arg).lower() in ("false","", "0", "0.0", "[]", "()", "{}"):
		return False
	else:
		return True

def _parse_milestone_inputs(inp):
	milestones = []
	for key in sorted(inp.keys()):
		if re.match("milestone_group[0-9]+", key):
			milestone_dict = {'key': key}
			milestone_name = key
			for param in inp[key].split(','):
				split_param = param.strip().split()
				if len(split_param) < 2: continue
				milestone_key = milestone_name + '_' + split_param[0]
				milestone_value = ' '.join(split_param[1:]) 
				milestone_dict[milestone_key] = milestone_value
			milestones.append(milestone_dict)
			#print(milestones)
	return milestones

def _generate_milestone_lists(milestones,):
	for milestone in milestones:
		milestone_pairs = []
		anchor_list = []
		key = milestone['key']
		milestone_vals_string = milestone['%s_milestone_values' % key]
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
		milestone['%s_pair_list' %milestone['key']] = milestone_pairs
		
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
	
def _generate_filetree(inp, sys_params):
	print("Creating file tree at location: %s" % inp['rootdir'])
	if sys_params['empty_rootdir'] and os.path.isdir(inp['rootdir']): # then we're in test mode, delete anything that is in the rootdir
		print("'empty_rootdir' set to True. Deleting all subdirectories/files in rootdir:", inp['rootdir'])
		shutil.rmtree(inp['rootdir']) # remove all the files within the grid subdirectory
	if not os.path.isdir(inp['rootdir']): # if the rootdir doesn't exist yet
		os.mkdir(inp['rootdir']) # then create it
	return

def _get_filetree_settings(anchor_list):
	filetree_settings = {
		'anchor_list':anchor_list,
	}

	return filetree_settings


def _write_milestone_file(anchor_list, temperature, md_time_factor, bd_time_factor, milestone_filename ):

	ourdoc = Document() # create xml document
	root = ourdoc.createElement("root") # create the xml tree
	ourdoc.appendChild(root)

	xmltemp = ourdoc.createElement("temperature") # include temperature
	root.appendChild(xmltemp)
	xmltemptext = ourdoc.createTextNode(str(temperature))
	xmltemp.appendChild(xmltemptext) 

	# MD time factor
	xmlmd_time_factor = ourdoc.createElement("md_time_factor") # include temperature
	root.appendChild(xmlmd_time_factor)
	xmlmd_time_factor_text = ourdoc.createTextNode(str(md_time_factor))
	xmlmd_time_factor.appendChild(xmlmd_time_factor_text)
	
	# BD time factor
	xmlbd_time_factor = ourdoc.createElement("bd_time_factor") # include temperature
	root.appendChild(xmlbd_time_factor)
	xmlbd_time_factor_text = ourdoc.createTextNode(str(bd_time_factor))
	xmlbd_time_factor.appendChild(xmlbd_time_factor_text)
	
	site_counter = 0

	##Currently only configured for a single site! TODO
	xmlsite = ourdoc.createElement("site")
	root.appendChild(xmlsite)
	xmlsitename = ourdoc.createElement("name")
	xmlsite.appendChild(xmlsitename)
	xmlsitetext = ourdoc.createTextNode("site_0")
	xmlsitename.appendChild(xmlsitetext)

	anchor_counter = 0
	for anchor in anchor_list:
		xmlanchor = ourdoc.createElement("anchor")
		xmlsite.appendChild(xmlanchor)
		xmlanchorname = ourdoc.createElement("name")
		xmlanchor.appendChild(xmlanchorname)
		xmlanchortext = ourdoc.createTextNode(str(anchor['name']))
		xmlanchorname.appendChild(xmlanchortext)

		xmlanchorindex = ourdoc.createElement("index")
		xmlanchor.appendChild(xmlanchorindex)
		xmltext = ourdoc.createTextNode(str(anchor_counter))
		xmlanchorindex.appendChild(xmltext)

		xmlanchordir = ourdoc.createElement("directory")
		xmlanchor.appendChild(xmlanchordir)
		xmltext = ourdoc.createTextNode(str(anchor['directory']))
		xmlanchordir.appendChild(xmltext)

		xmlmd = ourdoc.createElement("md")
		xmlanchor.appendChild(xmlmd)
		xmltext = ourdoc.createTextNode("True")
		xmlmd.appendChild(xmltext)

		xmlbd = ourdoc.createElement("bd")
		xmlanchor.appendChild(xmlbd)
		xmltext = ourdoc.createTextNode("False")
		xmlbd.appendChild(xmltext)

		for milestone_group in anchor['milestone_params']:
			xmlmilestone = ourdoc.createElement("milestone")
			xmlanchor.appendChild(xmlmilestone)

			xmlelement = ourdoc.createElement("id")
			xmlmilestone.appendChild(xmlelement)
			xmltext = ourdoc.createTextNode(str(milestone_group['lower_milestone_index']))
			xmlelement.appendChild(xmltext)

			xmlelement = ourdoc.createElement("group")
			xmlmilestone.appendChild(xmlelement)
			xmltext = ourdoc.createTextNode(str(milestone_group['milestone_group']))
			xmlelement.appendChild(xmltext)

			xmlelement = ourdoc.createElement("value")
			xmlmilestone.appendChild(xmlelement)
			xmltext = ourdoc.createTextNode(str(milestone_group['lower_bound']))
			xmlelement.appendChild(xmltext)

			xmlmilestone = ourdoc.createElement("milestone")
			xmlanchor.appendChild(xmlmilestone)

			xmlelement = ourdoc.createElement("id")
			xmlmilestone.appendChild(xmlelement)
			xmltext = ourdoc.createTextNode(str(milestone_group['upper_milestone_index']))
			xmlelement.appendChild(xmltext)

			xmlelement = ourdoc.createElement("group")
			xmlmilestone.appendChild(xmlelement)
			xmltext = ourdoc.createTextNode(str(milestone_group['milestone_group']))
			xmlelement.appendChild(xmltext)

			xmlelement = ourdoc.createElement("value")
			xmlmilestone.appendChild(xmlelement)
			xmltext = ourdoc.createTextNode(str(milestone_group['upper_bound']))
			xmlelement.appendChild(xmltext)




		anchor_counter += 1


	xml_string = ourdoc.toprettyxml(indent="  ")
	#print 'xml_string:', xml_string
	milestone_file = open(milestone_filename,'w')
	milestone_file.write(xml_string)
	milestone_file.close()

	return


def _group_milestones_to_anchor(milestones, anchor_dirlist, md_file_paths):
	anchor_index = 0
	anchor_list = []
	for anchor in anchor_dirlist:
		anchor_milestone_params = []
		milestone_params = md._gen_anchor_milestone_params(milestones, anchor_index)
		for milestone in milestone_params:
			anchor_milestone_params.append(milestone)
		anchor_list.append({'name' : anchor, 'directory' : md_file_paths[anchor_index]['prod'], 'milestone_params' : anchor_milestone_params,})
		anchor_index += 1

	return anchor_list

def _load_structures(inp, sys_params,):
	'''load protein/ligand structures needed for calculation'''
	parser = pdb.Big_PDBParser()
	print("now loading structures")

	# Read and/or create pickle files for the structures to save I/O time
	ligand_pkl_filename = os.path.join(inp['rootdir'], "ligand.pkl")
	receptor_pkl_filename = os.path.join(inp['rootdir'], "receptor.pkl")
	receptor_pkl_dry_filename = os.path.join(inp['rootdir'], "receptor_dry.pkl")
	receptor_pkl_dry_pqr_filename = os.path.join(inp['rootdir'], "receptor_dry_pqr.pkl")

	ligand=pickle_or_load(sys_params['lig_pqr_filename'], ligand_pkl_filename, struc_name="ligand", pqr=True)
	#receptor=pickle_or_load(sys_params['rec_pdb_filename'], receptor_pkl_filename, struc_name="receptor", pqr=False)
	receptor_dry=pickle_or_load(sys_params['rec_dry_pdb_filename'], receptor_pkl_dry_filename, struc_name="receptor_dry", pqr=False)
	receptor_dry_pqr=pickle_or_load(sys_params['rec_dry_pqr_filename'], receptor_pkl_dry_pqr_filename, struc_name="receptor_dry_pqr", pqr=True)

	struct={ # all parameters pertaining to structure
		'ligand':ligand,
		#'receptor':receptor,
		'receptor_dry':receptor_dry, # or create a function that will remove all waters, complicated by ions
		'receptor_dry_pqr':receptor_dry_pqr,
		'rec_com':pdb.center_of_mass(receptor_dry), # have to take into account the center of mass of the receptor itself
		'lig_com':pdb.center_of_mass(ligand), # have to take into account the center of mass of the ligand itself
	}

	return struct

def pickle_or_load(filename, picklename, struc_name="pickle",pqr=False):
	'''for large files, instead of parsing, they can be saved and loaded much more quickly as a pickle. '''
	parser = pdb.Big_PDBParser()
	if os.path.exists(picklename) and os.path.getmtime(picklename) > os.path.getmtime(filename): # if the pickle has been most recently modified
		# load the pickle
		print("reading pickle:", picklename)
		our_file=open(picklename, 'rb')
		our_obj=pickle.load(our_file)
		our_file.close()
	else:
		# then load the file itself and save the pickle
		our_obj=parser.get_structure(struc_name, filename, pqr=pqr, conventional=False) # load the file
		print("writing pickle:", picklename)
		our_file=open(picklename, 'wb')
		pickle.dump(our_obj, our_file, protocol=-1)
		our_file.close()
	return our_obj

def prepare_seekr():
	parser = argparse.ArgumentParser(description="top level SEEKR program used to prepare and generate all input files necessary for a SEEKR calculation")
	parser.add_argument('input_filename', metavar='INPUT_FILENAME', type = str, help="name of SEEKR input file")
	args = parser.parse_args() #parse command line arguments
	args = vars(args) #converts to dictionary

	_get_inputs(args,)
	sys_params = _get_sys_params(inputs) #parse general system parameters (structures, forcefield files, directory names from input file)
	md_milestones = _parse_milestone_inputs(inputs) #parse milestone CV parameters and milestone values from input file
	milestones, md_anchor_list = _generate_milestone_lists(md_milestones) #generate upper/lower bounds of a SINGLE CV for each anchor/Voronoi celli

	# TODO generate anchor lists for multiple milestone CV's
	# TODO BD milestones
	#print(milestones)
	_generate_filetree(inputs, sys_params) #creates/clears top level ditectory
	filetree_settings = _get_filetree_settings(md_anchor_list)
	md_filetree_settings_all = {**filetree_settings, **sys_params}
	anchor_dirlist, md_file_paths = filetree.md_filetree(md_filetree_settings_all)
	md_settings = _get_md_settings(inputs, md_file_paths) #parse parameters for MD simulations from input file
	md_settings_all = {**md_settings, **sys_params, **filetree_settings,}
	md.main(md_settings_all, milestones)
	print(milestones[-1])
	milestone_filename= os.path.join(sys_params['rootdir'], 'milestones.xml') 
	anchor_list = _group_milestones_to_anchor(milestones, anchor_dirlist, md_file_paths,)
	print('Anchor List',anchor_list)
	_write_milestone_file(anchor_list, md_settings['master_temperature'], 
		sys_params['md_time_factor'], sys_params['bd_time_factor'],milestone_filename)

	structures = _load_structures(inputs, sys_params)
	bd_settings = _get_bd_settings(inputs, sys_params, structures)
	bd_milestone = anchor_list[-1]
	print(bd_milestone)
	#bd_milestone_pair = bd_milestone['%s_pair_list' %bd_milestone['key']][-1]
	bd_lower_bound = bd_milestone['milestone_params'][0]['lower_bound']
	bd_lower_bound_index = bd_milestone['milestone_params'][0]['lower_milestone_index']
	bd_index = bd_milestone['milestone_params'][0]['upper_milestone_index']
	b_surf_distance = bd_milestone['milestone_params'][0]['upper_bound']
	#bd_milestone_index = 
	bd_settings.update({'b_surf_distance' : b_surf_distance,
		'bd_lower_bound' : bd_lower_bound,
		'bd_lower_bound_index' : bd_lower_bound_index,
		'bd_index': bd_index,
		})

	print("BD b surface distance", b_surf_distance)
	bd.main(bd_settings)
	#_write_milestone_file(anchor_list, md_settings['master_temperature'], 
	#	sys_params['md_time_factor'], sys_params['bd_time_factor'],milestone_filename)





if __name__ == "__main__": prepare_seekr()


