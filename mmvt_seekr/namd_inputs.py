#!/usr/bin/python

'''
Part of the SEEKR kinetic rate estimation package

contains all the messy information about NAMD min,equil,ensemble, and production simulation parameters

'''
import datetime, os
from adv_template import *
from math import sqrt
import pdb2 as pdb

self_path = os.path.dirname(os.path.realpath(__file__)) # get the path to this script

namd_input_template_location = os.path.join(self_path, 'namd_input.template')
colvar_input_template_location = os.path.join(self_path, 'colvar.template')
colvar_script_name = 'colvars.script'

class Namd_template(Adv_template):
	def __init__(self, template_filename, params):
		self.template_string = ''.join(open(template_filename, 'r').readlines()) # first load the template file and make it a string separated by newlines
		self.params = params
		self.template_string = self.fix_vars()
		rawoutput = self.parse_commands()
		template_string = string.Template(rawoutput) # create the template
		self.output = template_string.safe_substitute(params) # will substitute every indicated $variable in params. Safely in case template file contains extraneous ($$$) dollar signs
		return

	def input_gen(self, filename):
		'''generates input file called (filename) and fills the parameters from the dictionary params'''
		out_file = open(filename, 'w') # open output namd input file
		out_file.write(self.output)
		out_file.close()
		return

def _parse_pdb(filename):
	parser = pdb.Big_PDBParser() 
	struct=parser.get_structure(struct_name, filename, pqr=False, conventional=False) # load the file
	return struct

def cell_params(struct_name, shape="box"):
  '''returns the cellorigin and cellbasisvector parameters for a namd input file'''
  params = {}
  struct= _parse_pdb(struct_name)
  boxdims = pdb.minmax_width(struct) # measures the width of the waterbox
  boxctr =  pdb.center(struct)
  if shape == "box":
    params['cellbasisvector1'] = "%8f 0.0000000 0.0000000" % boxdims[0]
    params['cellbasisvector2'] = "0.0000000 %8f 0.0000000" % boxdims[1] # writes these values to 8 digits
    params['cellbasisvector3'] = "0.0000000 0.0000000 %8f" % boxdims[2]
  elif shape == "oct":
    d = boxdims[0] # assuming that the octahedron is aligned to the x-axis
    params['cellbasisvector1'] = "%8f 0.0000000 0.0000000" % d
    params['cellbasisvector2'] = "%8f %8f 0.0000000" % ((-1.0/3.0)*d, (2.0/3.0)*sqrt(2.0)*d) # writes these values to 8 digits
    params['cellbasisvector3'] = "%8f %8f %8f" % ((-1.0/3.0)*d, (-1.0/3.0)*sqrt(2.0)*d, (-1.0/3.0)*sqrt(6.0)*d)
  else:
    raise Exception("%s is not a valid ensemble option. Must be 'box' or 'oct'." % (shape,))
   
  params['cellorigin'] = '%8f %8f %8f' % boxctr
  return params
  
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

equil_params = {
	'caption':'SEEKR Anchor Equilibration', # the string at the top of the file
	'pme':'yes',
	#'fullelectfrequency':'1', # frequency (in timesteps) that electrostatics are evaluated
	#'veldcdfile':'${outname}vel.dcd',
	#'veldcdfreq':'1000',
}

prod_params = {
	'caption':'SEEKR MMVT Production', # the string at the top of the file
	'temperature':'',
	'pme':'yes',
}

def ensemble_params(ensemble, temp):
	'''populates the ensemble-relevant variables with ensemble information'''
	params = {}
	if ensemble == 'npt' or ensemble == 'nvt':
		params['langevin'] = 'on'
		params['langevintemp'] = temp
		params['langevindamping'] = '5'
		params['langevinhydrogen'] = 'no'

	if ensemble == 'npt':
		params['usegrouppressure'] = 'no'
		params['langevinpiston'] = 'on'
		params['langevinpistontarget'] = '1.01325'
		params['langevinpistonperiod'] = '100'
		params['langevinpistondecay'] = '50'
		params['langevinpistontemp'] = temp
	elif ensemble == 'nvt':
		params['langevinpiston'] = 'no'
		params['usegrouppressure'] = 'no'

	elif ensemble == 'nve':
		params['langevin'] = 'off'
		params['langevinpiston'] = 'off'
		params['usegrouppressure'] = 'no'
	else:
		raise Exception("%s is not a valid ensemble option. Must be 'npt', 'nvt', or 'nve'." % (ensemble,))
	return params

def write_freq_params(freq):
	'''specifies the dcd, xst, ... write frequency'''
	params = {}
	params['xstfreq'] = freq
	params['dcdfreq'] = freq
	params['restartfreq'] = freq
	params['outputenergies'] = freq
	params['outputtiming'] = freq
	return params

def _make_input(holo, stage, ff, write_freq, ensemble='nvt', get_cell=False, settings={}):
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
	#print(settings)
	params = {}
	params.update(default_namd_input_params)
	#ff = settings['ff']
	#holo= settings['system_pdb_filename']

	
	# forcefield
	if ff == 'amber':
		params.update(amber_params)
	elif ff == 'charmm':
		params.update(charmm_params)
	else:
		raise Exception("%s is not a valid ff option. Must be 'amber' or 'charmm'.")

	# SEEKR stage
	#if stage == 'min':
	#	params.update(min_params)
	#if stage in ['temp_equil', 'equil']:
	#	params.update(temperature_equil_params)
	if stage == 'equil':
		#write_freq = settings['equil_settings']['write_freq']
		#ensemble = settings['equil_settings']['ensemble']
		params.update(equil_params)
	if stage == 'prod':
		write_freq = settings['prod_settings']['write_freq']
		ensemble = settings['equil_settings']
		params.update(prod_params)

	params.update(write_freq_params(write_freq))
	if get_cell: params.update(cell_params(holo, get_cell)) # if the cell box dimensions needs to be specified, then call it
	#if fixed: params.update(fixed_params) # if the cell box dimensions needs to be specified, then call it
	#if constraints: params.update(constraint_params) # if the cell box dimensions needs to be specified, then call it

	#params['base'] = base
	
	temperature = str(settings['master_temperature'])
	params['temperature'] = temperature
	# ensemble
	params.update(ensemble_params(ensemble,temperature))
	# user-defined
	params.update(settings) # make sure to update with the user-specified parameters last, so they take precedence
	#print("params :", params)
	namd = Namd_template(namd_input_template_location, params) # generate the namd input class
	return (namd, params) # return the input file
	