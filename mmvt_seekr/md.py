import os
from shutil import copyfile
import namd_inputs

def _make_relative_path(oldpath, fromdir):
	'''given a relative 'oldpath' to a file, will construct a relative path from a 'fromdir' directory to the file'''
	oldabspath = os.path.abspath(oldpath)
	fromdirpath = os.path.abspath(fromdir)
	sep = os.path.sep
	oldpathlist = oldabspath.split(sep)
	fromdirlist = fromdirpath.split(sep)
	same = True
	for i in range(len(fromdirlist)): # for every directory on up to the fromdir
		if same and fromdirlist[0]==oldpathlist[0]: # if we are the same up to this point, and these two members are the same
			fromdirlist.pop(0)
			oldpathlist.pop(0)
		elif fromdirlist[0]!=oldpathlist[0]: # as soon as we have a divergence in the path tree
			same=False
			backsteps = len(fromdirlist) # we have to take a few steps back from the 'fromdir'
			backlist = ['..']*backsteps # create a list of '..' backwardses
			relpathlist = backlist + oldpathlist # append the '..'s to the remaining oldpathlist
			break

	return sep.join(relpathlist) # and join them with the '/' to make our relative path


def _prep_building(settings):
	i = settings['index']
	rootdir = settings['rootdir']
	system_pdb_filename = settings['system_pdb_filename']
	ff = settings['ff']
	system_params_filename = settings['system_params_filename']
	system_rst_filename = settings['system_rst_filename']
	building = settings['md_file_paths'][i]['building']	


	copyfile(system_pdb_filename, os.path.join(building,'holo.pdb'))
	if ff == 'amber':
		prmtop = os.path.join(building,'holo.parm7')
		inpcrd = os.path.join(building,'holo.rst7')
		copyfile(system_params_filename, prmtop)
		copyfile(system_rst_filename, inpcrd)
	if ff == 'charmm':
		copyfile(system_params_filename, os.path.join(building,'holo.psf'))

	return prmtop, inpcrd

# def _prep_equil(md_settings):
# 	inp, namd_params = namd_inputs._make_input(stage='equil',get_cell=False, settings=md_settings)
# 	input_name = os.path.join(md_settings['md_file_paths'][i]['equil'], 'equil_1.namd')
# 	inp.save(input_name)

		

# 	return

def prep(settings, stage, inpname, outname='', ):
	print("creating %s files" % stage)
	
	i = settings['index']
	path = settings['md_file_paths'][i][stage]
	ff = settings['ff'] # forcefield
	anchor_list = settings['anchor_list']
	temperature = settings['master_temperature']

	if not outname: outname = stage # provide the default outname
	get_cell = False

	if stage == 'equil':
		stage_settings = settings['equil_settings']
		if not stage_settings['namd_settings']['extendedsystem']: # if they didn't define a starting XSC file
			get_cell = settings['cell_shape'] # create the periodic boundary conditions
		stage_settings = settings['equil_settings']

	if stage == 'prod':
		stage_settings = settings['prod_settings']

	ensemble = stage_settings['ensemble']
	namd_settings = stage_settings['namd_settings'] # extra settings to send to the NAMD input file
	
	#settings['whoami'] = str(pos_milestone.index)
	#settings['siteid'] = str(pos_milestone.siteid)
	#settings['grid_edge_rad'] = '0.0' # temporary value until this can be more effectively predicted; will speed up TCL script evaluation
	if ff == 'amber':
		namd_settings['parmfile']=_make_relative_path(settings['prmtop'],path) # Find the relative location of the prmtop/inpcrd to be used in the mins
		namd_settings['ambercoor']=_make_relative_path(settings['inpcrd'],path)
	elif ff == 'charmm':
		namd_settings['amber'] = 'no'
		namd_settings['coordinates'] = '../holo_wet.pdb'
		namd_settings['structure'] = '../building/holo_wet.psf'
		namd_settings['paratypecharmm'] = 'on'
		counter = 1
		for parameter in settings['charmm_settings']['parameters']:
			if counter == 1: # we need to make sure the namd input file can get more than one parameter file
				i = ""
			else:
				i = str(counter)
			namd_settings['parameters%s' % i] = parameter
			counter += 1
			#namd_settings['parameters2'] = '/extra/moana/rswift1/Permeability/Colvar/A/1F/par_all36_lipid.prm'

	#fhpd_file = 'fhpd.txt'
	namd_settings['watermodel'] = settings['watermodel']
	if stage == 'equil':
		holo= settings['system_pdb_filename']
		namd_settings['inpfilename'] = inpname
		namd_settings['outfilename'] = outname + "_1"
		print(settings['equil_settings']['namd_settings']['write_freq'])
		write_freq= settings['equil_settings']['namd_settings']['write_freq']
		ensemble = settings['equil_settings']['ensemble']

		inp, namd_params=namd_inputs._make_input(holo, stage, ff,write_freq, ensemble=ensemble,  get_cell=get_cell, 
			settings=namd_settings) #...
		input_name = os.path.join(path, '%s_1.namd' % (stage))
		inp.save(input_name)
		param_filename=(os.path.join(path, 'namd_parameters.pkl'))
		param_file = open(param_filename, 'wb')
		pickle.dump(namd_params, param_file)
		param_file.close()
		return os.path.join(stage, namd_settings['outfilename']) # return ens_equil


# 	if stage == 'fwd_rev':
# 		prelim_string_template = '''set TEMP $temperature
# set NUM_REVERSALS $num_reversals
# global REV_FILENAME_BASE; set REV_FILENAME_BASE "REV_COMPLETED.txt"
# global FWD_FILENAME_BASE; set FWD_FILENAME_BASE "FWD_COMPLETED.txt"
# global RESTART_FREQ; set RESTART_FREQ $restart_freq
# global RUN_FREQ; set RUN_FREQ $run_freq
# set ENS_EQUIL_FIRST_FRAME $begin
# set ENS_EQUIL_STRIDE $stride
# set LAUNCHES_PER_CONFIG $launches_per_config
# set FRAME_CHUNK_SIZE $frame_chunk_size
# set UMBRELLA_GLOB_DIR "../ens_equil/" ;# the umbrella sampling directory
# set UMBRELLA_GLOB_NAME "ens_equil_0_?.dcd" ;# the umbrella sampling trajs
# set nr [numReplicas]
# set replica_id [myReplica] ;# get the ID of this replica'''
# 		hedron = positions_orient.get_hedron(settings['quat_hedron'])
# 		settings['care_about_self'] = 'False'
# 		settings['phase'] = 'reverse'
# 		settings['abort_on_crossing'] = 'True'
# 		#settings['fhpd_file'] = '../reverse/fhpd.txt'
# 		settings['max_num_steps'] = stage_settings['max_num_steps']
# 		settings['milestone_string'] = make_milestone_list(settings['raw_milestone_list'],hedron)
# 		namd_settings['tclforces_vars'] = tclforces(settings).get_output()

# 		num_frames = settings['ensemble_equil_settings']['namd_settings']['numsteps']
# 		dcd_write_freq = int(settings['ensemble_equil_settings']['namd_settings']['dcdfreq'])
# 		begin = stage_settings['extract_first'] # new parameter
# 		end = int(num_frames)/dcd_write_freq + 1
# 		stride = stage_settings['extract_stride'] # new parameter
# 		launches_per_config = stage_settings['launches_per_config']
# 		frame_chunk_size = stage_settings['frame_chunk_size']
# 		restart_freq = stage_settings['restart_freq'] # New parameter
# 		run_freq = stage_settings['run_freq'] # new parameter
# 		num_reversals = int((end - begin)/stride)

# 		namd_settings['inpfilename'] = inpname
# 		namd_settings['outfilename'] = outname+".$replica_id"
# 		namd_settings['coordinates'] = "../holo_wet.pdb" # this gets overwritten...
# 		namd_settings['ambercoor'] = '' # this has to be disabled when coordinates are presented
# 		namd_settings['velocities'] = ''
# 		if bool(stage_settings['extract_xst']):
# 			namd_settings['extendedsystem'] = "../ens_equil/ens_equil_0_1.restart.xsc"
# 		namd_settings['bincoordinates'] = ""
# 		namd_settings['binvelocities'] = ""
# 		namd_settings['temperature'] = temperature
# 		namd_settings['replicaUniformPatchGrids'] = 'on'
# 		namd_settings['numsteps'] = ''
# 		namd_settings['restartfreq'] = '$RESTART_FREQ'
# 		#namd_settings['id'] = i
# 		prelim_settings = {'temperature':temperature, 'num_reversals':num_reversals, 'restart_freq':restart_freq, 'run_freq':run_freq, 'begin':begin, 'stride':stride, 'launches_per_config':launches_per_config, 'frame_chunk_size':frame_chunk_size}
# 		prelim_string = Adv_template(prelim_string_template,prelim_settings).get_output() # fill in the missing values
# 		namd_settings['prelim_string'] = prelim_string
# 		post_string_settings = {}
# 		namd_settings['post_string'] = File_template(fwd_rev_template_location, post_string_settings).get_output()

# 		inp, namd_params=namd_inputs.make_input(holo, ff, stage, temperature, write_freq, fixed=fixed, ensemble=ensemble, get_cell=get_cell, constraints=const, settings=namd_settings) #...
# 		input_name = os.path.join(path, 'fwd_rev1.namd')
# 		inp.save(input_name) 
# 		param_filename=(os.path.join(path, 'namd_parameters.pkl'))   
# 		param_file = open(param_filename, 'wb')
# 		pickle.dump(namd_params, param_file)
# 		param_file.close()
# 		#make_fhpd_script(path, 'forward_template.namd', settings['fhpd_file'], number_of_neighbors=number_of_neighbors)
# 		return os.path.join(stage, namd_settings['outfilename'])


def main(md_settings):
	'''called by seekr, executes other necessary commands'''
	anchor_list = md_settings['anchor_list']
	ff= md_settings['ff']
	for i in range(len(anchor_list)):
		md_settings['index'] = i

		if ff == 'amber':
			md_settings['prmtop'], md_settings['inpcrd']= _prep_building(md_settings)
		else:
			print("only amber ff currently implemented")
		equil_out = prep(md_settings, stage='equil', inpname='', )
		#inp, namd_params = namd_inputs._make_input(stage='equil',get_cell=False, settings=md_settings)
		#input_name = os.path.join(md_settings['md_file_paths'][i]['equil'], 'equil_1.namd')
		#inp.save(input_name)	

		#f settings['min']: prep(settings, holo, stage='min',  inpname='') # if we are supposed to do minimizations
		#f settings['temp_equil']: last_out=prep(settings, holo, stage='temp_equil', inpname=os.path.join('..',last_out), temperatures=settings['temp_equil_settings']['temp_range']) # if we are supposed to do temperature equilibrations
		#f settings['ens_equil']: last_out=prep(settings, holo, stage='ens_equil', inpname=os.path.join('..',last_out)) # if we are supposed to do temperature equilibrations
		#f settings['ensemble_equil_settings']['colvars']: # if running colvars
		#colvars.main(settings)

		
		#prep(settings,holo, stage='reverse', inpname=os.path.join('..',last_out)) # this can't be run until the ensemble is complete
		#prep(settings,holo, stage='forward', inpname=os.path.join('..',last_out)) # this can't be run until the reverse stage is complete
		#rep(settings, holo, stage='fwd_rev', inpname=os.path.join('..',last_out)) # we can put the two stages together
	return