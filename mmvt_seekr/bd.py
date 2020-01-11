#!/usr/bin/python

'''
bd.py

creates the necessary files to run a BD simulation using Browndye


'''

import sys, os, math, shutil, subprocess, glob, random #, make_fxd
import numpy as np
import pdb2 as pdb
import copy  # needed to keep track of separate structure objects
import unittest
import re
import xml.etree.cElementTree as ET # for writing xml files
from xml.dom import minidom
from adv_template import Adv_template, File_template
import apbs

RXN_FILENAME = 'rxns.xml'
INPUT_FILENAME = 'input.xml'
DEFAULT_TEMP = 298.0
empty_pqrxml = "./empty.pqrxml"

default_browndye_params = {
	'root':{
		'solvent':{
			'dielectric':'78',
			'debye-length':'7.86',
			'kT':'1',
		},
		'output':'results.xml',
		'start-at-site':'true',
		'trajectory-file':'traj',
		'include-desolvation-forces':'true',
		'n-trajectories':'10000',
		'n-threads':'1',
		'molecule0': {
			'prefix':'prot0',
			'atoms':'prot0.pqrxml',
			'apbs-grids': {
				'grid':'prot0.dx',
			},
			'solute-dielectric':'2.0',
		},
		'molecule1': {
			'prefix':'prot1',
			'atoms':'prot1.pqrxml',
			'all-in-surface':'false',
			'apbs-grids': {
				'grid':'prot1.dx'
			},
			'solute-dielectric':'2.0',
		},
		'time-step-tolerances': {
			'minimum-dx':'0.20',
		},
		'reactions':'rxns.xml',
		'seed':'111113',
		'n-trajectories-per-output':'1000',
		'n-copies':'200',
		'n-bin-copies':'200',
		'n-steps':'1000000',
		'n-steps-per-output':'1000',
		'max-n-steps':'1000000',
		'n-threads':'1',
		#'min-rxn-coord-file':'min_coord',
	},
}

default_browndye_molecule_block = { # for each diffusing molecule, creates this block
      'prefix':'prot0',
      'atoms':'prot0.pqrxml',
      'apbs-grids': {
        'grid':'prot0.dx',
      }
}

def prettify(elem):
	'''return a pretty-printed xml string for the Element'''
	rough_string = ET.tostring(elem, 'utf-8')
	reparsed = minidom.parseString(rough_string)
	return reparsed.toprettyxml(indent="  ")


def _write_browndye_input(pqrs,settings,criteria,work_dir='.',browndye_bin='', start_at_site='true',):
	'''
	generates a Browndye input file
	'''
	counter = 0
	input_xml = default_browndye_params
	debyes = []
	apbs_settings = settings['apbs_settings']
	inputgen_settings = settings['inputgen_settings']
	rxn_criteria = make_rxn_criteria(criteria,pqrs) # modifies the molecules to contain ghost atoms (for BrownDye rxn criteria)
	pqrxmls = []
	for pqr in pqrs: # for each molecule in pqr format
		prefix = pqr.struct_id # name of the molecule
		#print "PREFIX", prefix
		pqrfile = os.path.join(work_dir, prefix+'.pqr')
		pqr.save(pqrfile,pqr=True,endmdl=False)
		#print "pqrfile:", pqrfile
		dxfile, debye = apbs.main(pqrfile, inputgen_settings=inputgen_settings, apbs_settings=apbs_settings,) # get the electrostatic grid and debye length for the molecule
		debyes.append(debye)
		pqrxmlfile = pqr2xml(pqrfile, pqr2xml_program=os.path.join(browndye_bin, 'pqr2xml')) # call the pqrxml program using the Browndye software suite
		pqrxmls.append(pqrxmlfile)
		molecule_xml = copy.deepcopy(default_browndye_molecule_block) # create a copy of the molecule block, keep the default solute dielectric
		molecule_xml['apbs-grids']['grid']=os.path.basename(dxfile)
		molecule_xml['prefix'] = prefix
		molecule_xml['atoms'] = os.path.basename(pqrxmlfile)
		molecule_xml['solute-dielectric'] = '2.0'
		input_xml['root']['molecule%d'%counter] = molecule_xml # copy the molecule tree over to the input xml

		counter += 1

	#for pqr in pqrs:
		#pqr.save() # save the pqr files containing the ghost molecules
	rxn_file = open(os.path.join(work_dir,RXN_FILENAME), 'w')
	rxn_file.write(rxn_criteria)
	rxn_file.close()
	input_xml['root']['seed'] = int(random.random() * 10000) # give each one a random seed to make the simulation rounds different
	input_xml['root']['solvent']['debye-length'] = str(max(map(int,map(float,debyes)))) # set the debye length
	input_xml['root']['reactions'] = RXN_FILENAME
	input_xml['root']['solvent']['kT'] = float(settings['temperature']) / DEFAULT_TEMP
	input_xml['root']['n-threads'] = settings['threads']
	input_xml['root']['n-trajectories'] = settings['prods_per_anchor']
	input_xml['root']['start-at-site'] = start_at_site
	# if fhpd_mode:
	# 	input_xml['root']['seed'] = "$RANDOM"
	# 	input_xml['root']['include-desolvation-forces'] = "true"
	# 	del input_xml['root']['molecule1']['apbs-grids'] # remove the apbs-grid because we have no desolvation forces
	# 	input_xml['root']['molecule0']['prefix'] = "$REC"
	# 	input_xml['root']['molecule0']['atoms'] = "$REC.pqrxml"
	# 	input_xml['root']['molecule0']['apbs-grids']['grid'] = "$RECDX"
	# 	input_xml['root']['molecule1']['prefix'] = "$LIG"
	# 	input_xml['root']['molecule1']['atoms'] = "$LIG.pqrxml"
	# 	input_xml['root']['reactions'] = "$RXN"
	# 	input_xml['root']['n-trajectories'] = "$NTRAJ"


	# set number of steps,copies,trajectory outputs, etc
	# print "input_xml:", input_xml
	input_text = dict2xml(input_xml).text() # generate xml text for the file
	input_xml = {}
	input_file = open(os.path.join(work_dir,INPUT_FILENAME), 'w')
	input_file.write(input_text) # write an xml file for the input to bd
	input_file.close()

	return pqrxmls

class dict2xml(object):


    def __init__(self, structure):
        if len(structure) == 1:
            self.doc     = Document()
            rootName    = str(structure.keys()[0])
            self.root   = self.doc.createElement(rootName)

            self.doc.appendChild(self.root)
            self.build(self.root, structure[rootName])

    def build(self, father, structure):
        if type(structure) == dict:
            for k in structure:
                tag = self.doc.createElement(k)
                father.appendChild(tag)
                self.build(tag, structure[k])

        elif type(structure) == list:
            grandFather = father.parentNode
            tagName     = father.tagName
            grandFather.removeChild(father)
            for l in structure:
                tag = self.doc.createElement(tagName)
                self.build(tag, l)
                grandFather.appendChild(tag)

        else:
            data    = str(structure)
            tag     = self.doc.createTextNode(data)
            father.appendChild(tag)

    def display(self):
        print self.doc.toprettyxml(indent="  ")

    def text(self):
        return self.doc.toprettyxml(indent="  ")

def pqr2xml(pqrfile, pqr2xml_program='pqr2xml'):
  '''simply runs the pqr2xml program that comes with Browndye.'''
  #print "pqr2xml_program:", pqr2xml_program
  no_ext = os.path.splitext(pqrfile)[0] # get everything but the extension
  xmlfile = no_ext + '.pqrxml'
  command = '%s < %s > %s' % (pqr2xml_program, pqrfile, xmlfile)
  #if verbose: print("now running command:", command)
  result = os.system(command) # run the pqr2xml program
  if result != 0: raise Exception("There was a problem running pqr2xml")
  return xmlfile

def make_rxn_criteria(criteria,pqrs):
	''' takes the rxn criteria pairs and makes ghost atoms in each structure, then populates a rxn file

	criteria format:
		= [[(coords of pqr1), (coords of pqr2), distance], next...]

	example:
		= [[(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), 7.0],]
	'''
	# create the xml file
	roottag = ET.Element("roottag")
	first_state = ET.SubElement(roottag, "first-state")
	first_state.text = "start"
	reactions = ET.SubElement(roottag, "reactions")
	reaction_list = []
	#pqrs = map(copy.deepcopy, raw_pqrs)
	ghost0_coords_already_used = []
	ghost0_ids = {}
	ghost1_coords_already_used = []
	ghost1_ids = {}
	for i in range(len(criteria)):
		#coord0 = criteria[i][0]
		#coord1 = criteria[i][1]
		radius = criteria[i]['radius']
		# We don't have any duplicate ghost atoms
		if (criteria[i]['centerx'], criteria[i]['centery'], criteria[i]['centerz']) in ghost0_coords_already_used: # check if there is already a ghost atom in the molecule
			ghost0_id = ghost0_ids[(criteria[i]['centerx'], criteria[i]['centery'], criteria[i]['centerz'])]
		else: # then the ghost atom doesn't already exist
			ghost0_id = create_ghost_atom_in_pqr(pqrs[0], criteria[i]['centerx'], criteria[i]['centery'], criteria[i]['centerz']) # the star spreads the coords out into individual arguments
			ghost0_coords_already_used.append((criteria[i]['centerx'], criteria[i]['centery'], criteria[i]['centerz']))
			ghost0_ids[(criteria[i]['centerx'], criteria[i]['centery'], criteria[i]['centerz'])] = ghost0_id

		# check the ligand for already existing ghost atoms
		if (criteria[i]['ligx'], criteria[i]['ligy'], criteria[i]['ligz']) in ghost1_coords_already_used: # check if there is already a ghost atom in the molecule
			ghost1_id = ghost1_ids[(criteria[i]['ligx'], criteria[i]['ligy'], criteria[i]['ligz'])]
		else: # then the ghost atom doesn't already exist
			ghost1_id = create_ghost_atom_in_pqr(pqrs[1], criteria[i]['ligx'], criteria[i]['ligy'], criteria[i]['ligz'])
			ghost1_coords_already_used.append((criteria[i]['ligx'], criteria[i]['ligy'], criteria[i]['ligz']))
			ghost1_ids[(criteria[i]['ligx'], criteria[i]['ligy'], criteria[i]['ligz'])] = ghost1_id


		siteid = criteria[i]['siteid']
		index = criteria[i]['index']
		# now add to the xml file
		reaction_list.append(ET.SubElement(reactions, "reaction"))
		name=ET.SubElement(reaction_list[-1],"name")
		name.text = "%s_%s" % (siteid,index)
		state_before = ET.SubElement(reaction_list[-1], "state-before")
		state_before.text = "start"
		state_after = ET.SubElement(reaction_list[-1], "state-after")
		state_after.text = "end"
		criterion = ET.SubElement(reaction_list[-1], "criterion")
		n_needed = ET.SubElement(criterion, "n-needed")
		n_needed.text = "1"
		pair = ET.SubElement(criterion, "pair")
		atoms = ET.SubElement(pair, "atoms")
		atoms.text = "%d %d" % (ghost0_id, ghost1_id)
		distance = ET.SubElement(pair, "distance")
		distance.text = "%f" % (radius,)

	criteria_xml = prettify(roottag)
	return criteria_xml




def main(settings):

	rec_struct = settings['rec_struct']
	lig_struct = settings['lig_struct']
	browndye_bin = settings['browndye_bin_dir']
	empty_pqrxml = os.path.abspath(settings['empty_pqrxml_path'])

	pqrs = [copy.deepcopy(rec_struct), copy.deepcopy(lig_struct)]
	pqrs[1].struct_id='bd_ligand'
	b_surface_path = settings['b_surface_path']
	if not os.path.exists(b_surface_path): os.mkdir(b_surface_path)
	b_surface_criteria = []

	for site in settings['b_surface_ending_surfaces']:
		b_surface_criteria.append({'centerx':site['x'], 'centery':site['y'], 'centerz':site['z'], 'ligx':lig_center[0], 
			'ligy':lig_center[1], 'ligz':lig_center[2], 'radius':site['radius'], 'index':site['index'], 'siteid':site['siteid']}) # add every site to the criteria list
	print("bsurface_criteria:", b_surface_criteria)
	b_surface_pqrxmls = _write_browndye_input(pqrs, settings, b_surface_criteria, work_dir=b_surface_path,	browndye_bin=browndye_bin, start_at_site='false',) # write input for this part



	return


