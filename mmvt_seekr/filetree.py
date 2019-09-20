#!/usr/bin/python

'''
filetree.py

Contains functions to be called by seekr.py that generate the necessary filetree to execute simulations.

'''
import sys, os, math, shutil, subprocess #, make_fxd
import unittest
import numpy as np
#import pdb2 as pdb

verbose = True
mdtree = {'building':{}, 'equil':{}, 'prod':{}} #'test':{'dir1':{'dir1_1':{}, 'dir1_2':{}}, 'dir2':{}}})
bdtree = {'prod':{}}


class Filetree():
	'''defines a file tree: a framework of directories to be populated with files'''
	def __init__(self, selflist):
		self.tree = selflist
		return

	def make_tree(self, rootdir, branch={}):
		'''will construct the file tree given a root directory, populating
it with all subsequent branches and leaves in the file tree'''
		#print 'rootdir', rootdir
		assert os.path.isdir(rootdir), "rootdir argument must be a real directory"
		if not branch: branch = self.tree
		for subbranch in branch.keys():
			# first create each subbranch
			subbranchpath = os.path.join(rootdir,subbranch)
			if not os.path.exists(subbranchpath):
				os.mkdir(subbranchpath)
			if not branch[subbranch]: # if its an empty list, then we have a leaf
				continue
			else: # then we can descend further
				self.make_tree(subbranchpath, branch=branch[subbranch])
		return




def md_filetree(settings):
	rootdir = settings['rootdir']
	#wet_configs=settings['wet_configs']
	anchor_list = settings['anchor_list']
	#milestone_pos_rot_list = settings['milestone_pos_rot_list']
	#dry_configs=settings['dry_configs']

	anchor_dirlist = []
	md_file_paths = []
	#bd_file_paths = []
	anchor_counter = 0
	
	for i in range(len(anchor_list)): # for each configuration, it gets its own directory
		anchor = anchor_list[i]
		anchor_name = "anchor_%s" % (i)
		anchor_filetree = Filetree({anchor_name:{}})
		if verbose: print("Creating directory:", anchor_name)
		anchor_filetree.make_tree(rootdir) # create this anchor's directory
		#milestone_pos_rot_list[i][0].directory = anchor_name # update this milestones directory information
		anchor_dir = os.path.join(rootdir, anchor_name)
		# MD filetree
		md_file_path={}
		md_dir=os.path.join(anchor_dir,'md') # directory for MD
		md_filetree=Filetree({'md':mdtree})
		md_filetree.make_tree(anchor_dir) # make the MD filetree
		for key in mdtree.keys():
			md_file_path[key] = os.path.join(md_dir,key)
		#wet_holo_filename=os.path.join(md_dir,'holo_wet.pdb')
		#wet_config.save(wet_holo_filename, amber=True, standard=False) # write the holo structure into the md directory
		#md_file_path['wet_holo'] = wet_holo_filename
		


		anchor_dirlist.append(anchor_name)
		md_file_paths.append(md_file_path)
		#if bd_dir: bd_file_paths.append(bd_dir) # bd_file_path
		anchor_counter+=1 # increment the loop counter


	return anchor_dirlist, md_file_paths,