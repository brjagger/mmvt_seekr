# Parameters from each dictionary updated in the final NAMD parameters in the order defined below...
default_namd_input_params = {
'caption':'default/',
'date':str(datetime.date.today()),

'inpdir':'',
'inpfilename':'dummy',
'outdir':'',
'outputname':'dummy',
'firsttimestep':'0',
'timestep':'2.0',
'numsteps':'1',

'amber':'yes',
'parmfile':'dummy.prmtop',
'ambercoor':'dummy.inpcrd',
'readexclusions':'yes',
'scnb':'2.0',
'exclude':'scaled1-4',
'_1_4scaling':'0.833333', # NOTE: the '-' has been changed to a '_'. Also, variable may not begin with numerical char
'watermodel':'',

'gromacs':'off',
'grotopfile':'',
'grocoorfile':'',


'coordinates':'',
'structure':'',
'parameters':'CHARMM_parameter_file',
'paratypexplor':'on',
'paratypecharmm':'off',
'velocities':'',
'binvelocities':'$inpname.restart.vel',
'bincoordinates':'$inpname.restart.coor',
'cwd':'',
'temperature':'298',
'watermodel':'',

'outfilename':'$outname',
'binaryoutput':'yes',
'restartname':'$outname.restart',
'restartfreq':'1000',
'restartsave':'no',
'binaryrestart':'yes',
'dcdfile':'$outname.dcd',
'dcdfreq':'1000',
'dcdunitcell':'',
'veldcdfile':'',
'veldcdfreq':'',
'forcedcdfile':'',
'forcedcdfreq':'',

'outputenergies':'1000',
'mergecrossterms':'',
'outputmomenta':'',
'outputpressure':'',
'outputtiming':'10000',
'usegrouppressure':'no',

'extendedsystem':'$inpname.restart.xsc',


'cutoff':'8',
'switching':'off',
'switchdist':'',
'zeromomentum': 'on',
'ljcorrection': '',

'pme':'yes',
'pmegridsizex':'',
'pmegridsizey':'',
'pmegridsizez':'',
'pmegridspacing':'',

'rigidbonds':'all',
'rigiditerations':'100',
'rigidtolerance':'1e-8',
'usesettle':'on',

'constraints':'',
'consref':'restrain_backbone_ref.pdb',
'conskfile':'restrain_backbone_ref.pdb',
'conskcol':'O',
'constraintscaling':'1.0',

'cellbasisvector1':'',
'cellbasisvector2':'',
'cellbasisvector3':'',
'cellorigin':'',
'xstfile':'$outname.xst'  ,
'xstfreq':'10000',
'wrapwater':'on',
'wrapall':'off',
'wrapnearest':'on',

'minimization':'off',
'minimize':'5000',
'fixedatoms': 'off',
'fixedatomsfile': '',
'fixedatomscol': '',

'nonbondedfreq':'1',
'fullelectfrequency':'1',

'langevin':'on',
'langevinfile':'',
'langevincol':'O',
'langevintemp':'300',
'langevindamping':'5',
'langevinhydrogen':'no',

'useflexiblecell':'no',
'langevinpiston':'on',
'langevinpistontarget':'1.01325',
'langevinpistonperiod':'100',
'langevinpistondecay':'50',
'langevinpistontemp':'300',
'useconstantarea':'no',

'tclforces':'off',
'tclforcesscript':'',

'colvars':'on',
'colvarsscript':'colvars.script',
'colvarsconfig':'',


'pairlistdist':'11',
'stepspercycle':'20',
'margin':'0.0',

'timestep' : '2.0',
'rigidbonds' : 'all',
'useflexiblecell':'no',
'usesettle' : 'on', # since rigidbonds are on, faster than the SHAKE algorithm
'useconstantarea': 'no',
'pmegridspacing':'1.0',

}

charmm_params = { 
  'amber':'no',
  '_1_4scaling' : '1.0',
  'paraTypeXplor': 'on',
  'cutoff' : '10.0',
  'switchdist' : '10.0',
  'pairlistdist' : '14.0',
}

amber_params = { # based on "Using the Amber force field in NAMD" at http://ambermd.org/namd/namd_amber.html by G. Giambasu and D. Case
  'amber':'yes',
  'readexclusions': 'yes',
  'exclude':'scaled1-4',
  '_1_4scaling' : '0.833333',
  'scnb':'2.0',
  'switching':'off',
  #'switchdist':'12', # if you ever turn on switching, you'll need to define these values
  #'pairlistdist':'11',
  'cutoff':'8',
  'fullelectfrequency' : '1',
  'stepspercycle':'10',
  'ljcorrection':'on',
}

def make_input(holo, ff, stage, temperature, write_freq, receptor_type='globular', ensemble='nvt', base='namd input', get_cell=False, fixed=False, constraints=False, settings={}):
	'''creates NAMD input files based on the settings specified:
	arguments:
	ff: can be one of 'amber' or 'charmm'
	stage: can be 'min', 'equil', 'ens', 'prod' # might not actually do this one
	temperatures: a list of temperatures to create simulation files for. Useful for temperature equilibration
	write_freq: the frequency that trajectories, energies, etc. are written
	receptor_type: can be 'membrane' or 'globular'
	ensemble: can be 'nve', 'nvt', 'npt'
	base: a descriptive string added to the input file header
	get_cell: whether the structure is parsed to find waterbox dimensions to populate the cellcenter and cellbasisvector params of the input file
	settings: additional namd parameters that will be substituted into the namd input file
	'''
	params = {}
	params.update(default_namd_input_params)
	
	# forcefield
	if ff == 'amber':
		params.update(amber_params)
	elif ff == 'charmm':
		params.update(charmm_params)
	else:
		raise Exception, "%s is not a valid ff option. Must be 'amber' or 'charmm'."


	# receptor_type
	if receptor_type == "globular":
		params.update(globular_params)
	elif receptor_type == "membrane":
		params.update(membrane_params)
	else:
		raise Exception, "%s is not a valid receptor_type option. Must be 'globular', or 'membrane'."
	
	# SEEKR stage
	if stage == 'min':
		params.update(min_params)
	if stage in ['temp_equil', 'equil']:
		params.update(temperature_equil_params)
	if stage == 'ens_equil':
		params.update(ens_equil_params)
	if stage in ['prod','reverse','forward','fwd_rev']:
		params.update(prod_params)

	params.update(write_freq_params(write_freq))
	if get_cell: params.update(cell_params(holo, get_cell)) # if the cell box dimensions needs to be specified, then call it
	if fixed: params.update(fixed_params) # if the cell box dimensions needs to be specified, then call it
	if constraints: params.update(constraint_params) # if the cell box dimensions needs to be specified, then call it

	params['base'] = base
	
	temperature = str(temperature)
	params['temperature'] = temperature
	# ensemble
	params.update(ensemble_params(ensemble,temperature))
	# user-defined
	params.update(settings) # make sure to update with the user-specified parameters last, so they take precedence
	#print "params after ensemble:", params
	namd = Namd_template(namd_input_template_location, params) # generate the namd input class
	return (namd, params) # return the input file
	