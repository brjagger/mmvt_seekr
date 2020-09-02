import os
from shutil import copyfile
import namd_inputs
import pickle
from adv_template import Adv_template, File_template

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

def _gen_anchor_milestone_params(milestones, anchor_index, end=False):
	anchor_milestone_params= []

	for milestone in milestones:
		milestone_param_dict = {}
		milestone_pair = milestone['%s_pair_list' %milestone['key']][anchor_index]
		milestone_param_dict['milestone_group'] = milestone['key']
		milestone_param_dict['milestone_pair'] = milestone_pair
		milestone_param_dict['lower_bound'] = float(milestone_pair[0])
		milestone_param_dict['upper_bound'] = float(milestone_pair[1])
		milestone_param_dict['lower_anchor'] = (anchor_index -1)
		milestone_param_dict['upper_anchor'] = (anchor_index +1)
		milestone_param_dict['lower_milestone_index'] = (anchor_index) ## may need to be changed for multidimensional milestones
		milestone_param_dict['upper_milestone_index'] =  (anchor_index +1)
		milestone_param_dict['curr_anchor'] = anchor_index
		if end:
			milestone_param_dict['lower_end'] = 'false'
			milestone_param_dict['upper_end'] = 'true'
		else:
			milestone_param_dict['lower_end'] = 'false'
			milestone_param_dict['upper_end'] = 'false'
		anchor_milestone_params.append(milestone_param_dict)

	return anchor_milestone_params

def _prep_building(settings):
	i = settings['index']
	rootdir = settings['rootdir']
	system_pdb_filename = settings['system_pdb_filename']
	ff = settings['ff']
	system_params_filename = settings['system_params_filename']
	system_rst_filename = settings['system_rst_filename']
	system_bincoor_filename = settings['system_bin_coordinates']
	system_extended_system_filename = settings['extendedsystem']
	building = settings['md_file_paths'][i]['building']	


	copyfile(system_pdb_filename, os.path.join(building,'holo.pdb'))
	if ff == 'amber':
		prmtop = os.path.join(building,'holo.parm7')
		inpcrd = os.path.join(building,'holo.rst7')
		bincoor= os.path.join(building,'holo.coor')
		xsc = os.path.join(building,'holo.xsc')
		copyfile(system_params_filename, prmtop)
		copyfile(system_rst_filename, inpcrd)
		copyfile(system_bincoor_filename, bincoor)
		copyfile(system_extended_system_filename, xsc)
	if ff == 'charmm':
		copyfile(system_params_filename, os.path.join(building,'holo.psf'))

	return prmtop, inpcrd, bincoor, xsc

# def _prep_colvars(settings, milestones, filename='colvar.script'):
# 	print(milestones)
# 	for milestone in milestones:
# 		#colvar_type = milestone['colvar_type']
# 		milestone['colvar_filename'] = filename
# 		namd_inputs._make_colvars_input(milestone)


# 	return filename

def prep(settings, milestones, stage, inpname, outname='', ):
	#print("creating %s files" % stage)
	
	i = settings['index']
	path = settings['md_file_paths'][i][stage]
	ff = settings['ff'] # forcefield
	anchor_list = settings['anchor_list']
	temperature = settings['master_temperature']


	if not outname: outname = stage # provide the default outname
	get_cell = False

	if stage == 'equil':
		stage_settings = settings['equil_settings']
		#_prep_colvars(settings, milestones)
		if not stage_settings['namd_settings']['extendedsystem']: # if they didn't define a starting XSC file
			get_cell = settings['cell_shape'] # create the periodic boundary conditions
		stage_settings = settings['equil_settings']

	if stage == 'prod':
		stage_settings = settings['prod_settings']

	ensemble = stage_settings['ensemble']
	namd_settings = stage_settings['namd_settings'] # extra settings to send to the NAMD input file
	namd_settings['master_temperature'] = temperature

	#settings['whoami'] = str(pos_milestone.index)
	#settings['siteid'] = str(pos_milestone.siteid)
	#settings['grid_edge_rad'] = '0.0' # temporary value until this can be more effectively predicted; will speed up TCL script evaluation
	if ff == 'amber':
		namd_settings['parmfile']=_make_relative_path(settings['prmtop'],path) # Find the relative location of the prmtop/inpcrd to be used in the mins
		namd_settings['ambercoor']=_make_relative_path(settings['inpcrd'],path)
		#namd_settings['bincoordinates']= _make_relative_path(settings['bincoor'],path)
		#namd_settings['extendedsystem']= _make_relative_path(settings['xsc'],path)
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

	namd_settings['watermodel'] = settings['watermodel']
	if stage == 'equil':
		print("creating equilibration files")
		holo= settings['system_pdb_filename']
		namd_settings['inpfilename'] = inpname
		namd_settings['outfilename'] = outname + "_1"
		namd_settings['bincoordinates']= _make_relative_path(settings['bincoor'],path)
		namd_settings['extendedsystem']= _make_relative_path(settings['xsc'],path)
		#print(settings['equil_settings']['namd_settings']['write_freq'])
		write_freq= settings['equil_settings']['namd_settings']['write_freq']
		ensemble = settings['equil_settings']['ensemble']

		colvars, colvars_params = namd_inputs._make_colvars_input(milestones, stage, i)
		colvars_name = os.path.join(path,'colvars.script')
		colvars.save(colvars_name)


		#print(namd_settings)
		inp, namd_params=namd_inputs._make_input(holo, stage, ff,write_freq, ensemble=ensemble,  get_cell=get_cell, 
			settings=namd_settings) #...
		input_name = os.path.join(path, '%s_1.namd' % (stage))
		inp.save(input_name)
		param_filename=(os.path.join(path, 'namd_parameters.pkl'))
		param_file = open(param_filename, 'wb')
		pickle.dump(namd_params, param_file)
		param_file.close()

		return os.path.join(stage, namd_settings['outfilename']) # return ens_equil


	if stage == 'prod':
		print("creating production files")
		holo= settings['system_pdb_filename']
		namd_settings['inpfilename'] = inpname
		namd_settings['outfilename'] = outname + "_1"
		#print(settings['equil_settings']['namd_settings']['write_freq'])
		write_freq= settings['prod_settings']['namd_settings']['write_freq']
		ensemble = settings['prod_settings']['ensemble']

		colvars, colvars_params = namd_inputs._make_colvars_input(milestones, stage, i)
		colvars_name = os.path.join(path,'colvars.script')
		colvars.save(colvars_name)

		post_string_settings = {
			'prod_steps' : namd_settings['numsteps'],
			'eval_stride' : namd_settings['eval_stride'],
		}

		anchor_milestone_params = _gen_anchor_milestone_params(milestones, i)
		##TODO fix for multiple milestones
		post_string_settings.update(anchor_milestone_params[0])

		post_string_template = '''set crossing 		0
set whoami 		none
set incubation_time 		0
set max_steps 		$prod_steps
set EVAL_STRIDE 		$eval_stride

set CURR_ANCHOR 		$curr_anchor
set LOWER_ANCHOR 		$lower_anchor
set LOWER_MILESTONE 		$lower_milestone_index
set LOWER_BOUND 		$lower_bound
set UPPER_ANCHOR 		$upper_anchor
set UPPER_MILESTONE 		$upper_milestone_index
set UPPER_BOUND 		$upper_bound





for {set stepnum 0} {$stepnum < $max_steps} {incr stepnum $eval_stride} {
  run $EVAL_STRIDE
  set cv_val [cv colvar milestone1 value]

  if {$cv_val <= $LOWER_BOUND} {
    if {$rev_last} {
      print "trajectory stuck"
      }
          puts "SEEKR: Cell Collision: current: $CURR_ANCHOR, new: $LOWER_ANCHOR, stepnum: $stepnum"
          if {$whoami != $LOWER_MILESTONE} {
            puts "SEEKR: Milestone Transition: anchor: $CURR_ANCHOR, source: $whoami, destination: $LOWER_MILESTONE, stepnum: $stepnum, incubation time: $incubation_time"
            set incubation_time 0
            }
        revert
        rescalevels -1
        checkpoint
        set rev_last True
        set whoami $LOWER_MILESTONE

    } elseif {$cv_val >= $UPPER_BOUND} {
      if {$rev_last} {
      print "trajectory stuck"
      }
        puts "SEEKR: Cell Collision: current: $CURR_ANCHOR, new: $UPPER_ANCHOR, stepnum: $stepnum"
        if {$whoami != $UPPER_MILESTONE} {
         puts "SEEKR: Milestone Transition: anchor: $CURR_ANCHOR, source: $whoami, destination: $UPPER_MILESTONE, stepnum: $stepnum, incubation time: $incubation_time"
         set incubation_time 0
         }
        revert
        rescalevels -1
        checkpoint
        set rev_last True
        set whoami $UPPER_MILESTONE

    } else {
       checkpoint
       set rev_last False
       }
   incr incubation_time $EVAL_STRIDE
}
'''

		post_string = Adv_template(post_string_template,post_string_settings).get_output() # fill in the missing values
		namd_settings['post_string'] = post_string
		inp, namd_params=namd_inputs._make_input(holo, stage, ff,write_freq, ensemble=ensemble,  get_cell=get_cell, 
			settings=namd_settings) #...
		input_name = os.path.join(path, '%s_1.namd' % (stage))
		inp.save(input_name)
		param_filename=(os.path.join(path, 'namd_parameters.pkl'))
		param_file = open(param_filename, 'wb')
		pickle.dump(namd_params, param_file)
		param_file.close()
		return os.path.join(stage, namd_settings['outfilename'])



def main(md_settings, milestones):
	'''called by seekr, executes other necessary commands'''
	anchor_list = md_settings['anchor_list']
	ff= md_settings['ff']
	for i in range(len(anchor_list)):
		md_settings['index'] = i
		print("creating files for anchor ", i)

		if ff == 'amber':
			md_settings['prmtop'], md_settings['inpcrd'], md_settings['bincoor'], md_settings['xsc']= _prep_building(md_settings)
		else:
			print("only amber ff currently implemented")
		equil_out = prep(md_settings, milestones, stage='equil', inpname=os.path.join('..','building/holo'), )
		prod_out = prep(md_settings, milestones, stage='prod', inpname=os.path.join('..',equil_out), )
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