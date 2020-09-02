#!/usr/bin/env python

#===============================================
# MODULE DOCSTRING
#=================================================

"""
MMVT SEEKR module for building a milestoning model

Extracts transition statistics from each anchor and generates a milestoning 
model 
"""

#=================================================
#GLOBAL IMPORTS
#=================================================

import os, sys, re, random, unittest, glob, tempfile
from string import Template
import xml.etree.ElementTree as ET


#===========================
#GLOBAL VARIABLES
#============================

k_boltz = 1.3806488e-23
R_GAS_CONSTANT = 0.0019872 # in kcal/mol*K
DEFAULT_TEMP = 0.0 # temperature to assign if not in the xml file
DEFAULT_MD_TIME_FACTOR = 2.0e-15
DEFAULT_BD_TIME_FACTOR = 1000.0
DEFAULT_ROOTDIR = "./"
DEFAULT_BDDIR = "./"
DEFAULT_SHAPE = "sphere"
GRID_TRANSITION = "SEEKR: Milestone Transition: "
GRID_TRANSITION_COMPILE = re.compile(GRID_TRANSITION)
INITIAL_TRANSITION = "source: none,"
INITIAL_TRANSITION_COMPILE = re.compile(INITIAL_TRANSITION)
VT_COLLISION = "SEEKR: Cell Collision: "
VT_COLLISION_COMPILE= re.compile(VT_COLLISION)
INITIAL_COLLISION = "SEEKR: Cell Collision: current: none"
INITIAL_COLLISION_COMPILE = re.compile(INITIAL_COLLISION)
FORWARD_OUTPUT_GLOB = "vt_milestoning_*.out.results"

#===============================================================
#Classes
#==============================================================

class Model():
   ''' Object representing the entire milestoning model.

   Contains information about all milestones and cell anchors, their connectivity, etc.


   Attributes
   ----------
   sites : list
      a list of all site classe contained in the model
   temperature : float
      temperature at which the model is constructed
   md_time_factor : float
      lenth of a single MD step in the model (seconds)
   bd_time_factor : float
      length of a single BD timestep (seconds)
   num_sites: int
      number of milestone sites contained in the model
   '''


   def __init__(self): # initialize all variables used in the model
      self.sites = []
      self.temperature = 0.0
      self.browndye_bin_dir = ''
      self.md_time_factor = 0.0
      self.bd_time_factor = 0.0
      self.num_sites = 0
      self.num_milestones = 0
      self.positional_milestones = []
      self.rotational_milestones = []
      self.b_surface_milestone = None
      #self.index_dict = {}
      #...

   def _add_site(self, site): 
      """ append a new milestone site to this model"""
      self.sites.append(site)
      self.num_sites += 1
      #self.num_anchors += site.num_milestones
      #...

   def _make_index_dict(self):
      pass

   def _make_directories(self):
      """ determines which milestones have data stored in which directories"""
      i = 0
      for site in self.sites:
        for milestone in site.milestones:
         if milestone.shape == "rotational":
           self.rotational_milestones.append(milestone) # keep track of all rotational milestones for iteration
         else:
           self.positional_milestones.append(milestone)
         i += 1
      if len(self.rotational_milestones) == 0: # if there are no rotational milestones (from, say, an old milestoning xml file)
        self.rotational_milestones.append(milestone(i, 'rotational', False, False, "rotation_%d" % i, anchor="[1, 0, 0, 0]")) # then create one, we need at least one rotational milestone
      for pos_milestone in self.positional_milestones: # for each of the positional milestones
        rot_i = 0
        for rot_milestone in self.rotational_milestones:
         #id_str = '%d_%d' % (site.index, milestone.index)
         directory_name = "anchor_%s_%d" % (pos_milestone.fullname, rot_i)
         print("directory_name:", directory_name)

class Site():
   ''' object representing a milestone site: a collection of related milestones that share, say, the same origin.

   Attributes
   ----------

   anchors : list 
      list of each milestone cell/anchor contained in the milestoning model
   num_anchors : int
      number of milestone anchors/cells
   index : int
      the index associated with this particular site in the model
   name : str
      the name associated with this particular site
   '''

   def __init__(self, index, name):
      self.anchors = []
      self.num_anchors = 0
      self.index = index
      self.name = name
      #self.center ?
      #..
  
   def _add_anchor(self, anchor):
      """ add a new anchor to the site class

      Parameters
      ----------

      anchor: model.Anchor
         an instance of the Anchor class with information for a particular Voronoi cell
      """

      self.anchors.append(anchor)
      self.num_anchors +=1

class Anchor():
   '''Object representing a single Voronoi cell
   
   Attributes
   ----------

   milestones : list
      list of milestone class instances for this anchor/cell 
   num_milestones : int
      number of milestones attributed to this anchor/cell
   index : int
      numerical index of the anchor/cell
   fullname : str
      complete name of this anchor
   directory : str
      directory name containing this anchor/cell data
   bd : boolean
      whether this anchor contains BD simulation data
   md : boolean
      whether this anchor contains MD simulation data
   site : int
      index of the site which contains this anchor
   sitename : str
      name of the site containing this anchor
   total_steps : int
      total number of MD or BD steps simulated for this anchor/cell
   offsets : list
      a list of the number of steps by which each individual results file should be offset to 
      determine the total simulation time
   '''
   def __init__(self, index, md, bd, fullname, directory, siteindex, sitename, coord=None):
      self.milestones = []
      self.num_milestones = 0
      self.index = index
      self.fullname = fullname
      self.directory = directory
      #self.coord = coord
      self.bd = bd
      self.md = md
      self.site = siteindex
      self.sitename = sitename
      self.total_steps = 0
      self.offsets = {}

   def _add_milestone(self, milestone):
      self.milestones.append(milestone)
      self.num_milestones += 1

   def _parse_md_transitions(self, ):
      """Finds all production simulation output files and parses for transition event strings.
         Creates Transition and Collision objects.

      Returns
      -----------
      self.total_steps : int
      the total number of steps taken for this anchor 

      """
      forward_dir_glob = os.path.join(self.directory,FORWARD_OUTPUT_GLOB)
      forward_output_filenames = sorted(glob.glob(forward_dir_glob))
      # read files and sift out the transition lines
      unsorted_transitions = []
      unsorted_collisions = []
      replicate = 1
      self.offsets[replicate] = 0
      for filename in forward_output_filenames:
        transition_lines = []
        vt_collisions = []
        file_max_steps = 0  
        print('parsing transitions from file:', filename)
        for line in open(filename,'r'):
         if re.match(GRID_TRANSITION_COMPILE, line):
           if re.search(INITIAL_TRANSITION_COMPILE, line): continue #we don't want to count a transition for the first milestone touched by the simulation  
           else:
            transition_lines.append(line)
         if re.match(VT_COLLISION_COMPILE, line):
           if re.match(INITIAL_COLLISION_COMPILE, line): continue
           else:
            vt_collisions.append(line)
        #pprint(transition_lines)
      # feed the lines into the Transition object
        for line in transition_lines:
         transition = Transition(line)
         transition.replicate = replicate
         unsorted_transitions.append(transition)

        for line in vt_collisions:
         collision = Collision(line)
         collision.replicate = replicate
         unsorted_collisions.append(collision)
        file_max_steps = collision.step
        
        self.total_steps += file_max_steps
        self.offsets[replicate+1] = self.offsets[replicate]+file_max_steps 
        replicate += 1
      self.transitions = unsorted_transitions
      self.collisions = unsorted_collisions
      print('anchor', self.index, self.offsets)
      return self.total_steps

   def _get_md_transition_statistics(self, md_time_factor=DEFAULT_MD_TIME_FACTOR, max_step=None, ):
      """Parse the transition data to obtain transition counts and times.

      Parameters
      ---------
      self : object
         the Anchor class object
      md_time_factor : float, optional
         time factor corresponding to length of a single MD timestep
      max_step : int
         maximum step number for parsing transitions, anything longer is neglected

      Returns
      ---------
      counts : dict
         dictionary containing all transition statistics for each anchor
      total_counts : dict
         dictionary containing the total transition counts for each anchor 
      total_times : dict
         dictionary containing the total time spent in each anchor/cell
      avg_times : dict
         dictionary containing the average time spent in each anchor/cell

      """
      counts = {} # all the sources and their destinations
      #cell_counts = {} # keeps track of transitions between voronoi cells
      #cell_times = {} #keeps track of total time simulated ina voronoi cell
      total_counts = {} # keeps track of all counts to any destination
      total_times = {} # the total time of all counts, to be averaged later
      avg_times = {} # the times to transition out of each source
      site_index = self.site
      site_name = self.sitename
      
      for transition in self.transitions:
    
        source = transition.src
        dest = transition.dest
        time = transition.time
        anchor = transition.anchor
        stepnum = transition.cur_step
        src_key = source
        dest_key = dest

        if max_step != None and int(stepnum + self.offsets[transition.replicate]) > max_step:
         break
        if self.index in list(counts.keys()):
         if src_key in list(counts[self.index].keys()):
           if dest_key in list(counts[self.index][src_key].keys()):
            counts[self.index][src_key][dest_key] += 1
           else:
            counts[self.index][src_key][dest_key] = 1
           total_times[self.index][src_key] += time * md_time_factor
           total_counts[src_key] += 1
         else:
           counts[self.index][src_key] = {dest_key:1}
           total_counts[src_key] = 1
           total_times[self.index][src_key] = (time * md_time_factor)
        else:
         counts[self.index]= {src_key: {dest_key: 1}}
         total_times[self.index] = {src_key : (time * md_time_factor)}
         total_counts[src_key] =1

      return counts, total_counts, total_times, avg_times

   def _get_md_vt_collisions(self, md_time_factor=DEFAULT_MD_TIME_FACTOR, max_step=None,):
    """Parse the transition data to obtain all collisions with cell boundaries.

    Parameters
    ---------
    self : object
       the Anchor class object
    md_time_factor : float, optional
         time factor corresponding to length of a single MD timestep
    max_step : int
        maximum step number for parsing transitions, anything longer is neglected

    Returns
    ---------
    cell_counts : dict
        dictionary containing all cell collision statistics for the anchor
    total_time : float
      total time spent in the anchor
    """
    cell_counts = {}
    assert(len(self.collisions) > 0), \
      "No collisions detected in anchor %s" % self.index
    for collision in self.collisions:
      if max_step != None and int(collision.step + self.offsets[collision.replicate]) > max_step:
       break

      curr_cell = collision.curr_cell
      new_cell = collision.new_cell   
      if curr_cell in list(cell_counts.keys()):
       if new_cell in list(cell_counts[curr_cell].keys()):
         cell_counts[curr_cell][new_cell] += 1
       else:
         cell_counts[curr_cell][new_cell] = 1
      else:
       cell_counts[curr_cell] = {new_cell:1} 
    print("cell_counts", cell_counts)
    total_time = collision.step *md_time_factor + self.offsets[collision.replicate] * md_time_factor 

    return cell_counts, total_time

   def get_bd_transition_statistics(self, results_filename="results.xml", bd_time=0.0):
    '''read the BD results.xml file for the anchors to determine transition statistics and 
    times for the BD stage

      Parameters
      ---------
      self : object
         the Anchor class object
      results_filename: string
          name of file containing BD transition results
      bd_time: float
         time factor by which a bd timestep is scaled

      Returns
      ---------
      counts : dict
         dictionary containing BD milestone transition statistics
      total_counts : dict
         dictionary containing the total BD transition counts  
      total_times : dict
         dictionary containing the total BD simulation time
      avg_times : dict
         dictionary containing the average BD simulation time 


    '''
    bd_results_filename = os.path.join(self.directory, results_filename)
    counts = {}
    total_counts = {}
    total_times = {}
    avg_times = {}
    #src_key = '%s_%s' % (self.sitename, self.index)
    src_key = self.index
    counts[src_key], total_times[src_key] = _parse_bd_results(bd_results_filename)
    total_counts[src_key] = 0
    for dest_key in counts[src_key].keys():
      total_counts[src_key] += counts[src_key][dest_key]
    if total_times[src_key] == None:
      avg_times[src_key] = bd_time
      total_times[src_key] = bd_time
    else:
      avg_times[src_key] = total_times[src_key] / total_counts[src_key]
    
    return counts, total_counts, total_times, avg_times

class Milestone():
   """ Object representing a single milestone

   Attributes
   ----------
   id : int 
      numberical index of the milestone
   shape : str
      shape of the milestone (spherical,etc.)
   end : bool
      whether this milestone is an end state
   radius : float
      radial distance of the milestone (priomarily for spherical milestones)
   transitions : list
      list containing all transition statistics associated with this milestone


   """
   def __init__(self, id, group, value=None, end="false" ):
      self.id = id # most of this information is read from the provided milestone.xml file
      self.group = group
      self.end = end.lower()
      #self.fullname = fullname
      #self.directory = directory
      #self.anchor = anchor
      #self.normal = normal
      self.value = value
      #self.bd = bd
      #self.md = md
      #self.site = siteindex
      #self.sitename = sitename
      self.transitions = [] # all the transition statistics associated with this milestone


class Collision():
   """Object representing a collision with a voronoi cell edge
   """
   def __init__(self,line):
      line = line.strip() # we have to parse the line    
      self.line = line
      linetail = line[len(VT_COLLISION)-1:] # just take the last bit of the line, the important part, but not the endline
      linelist = linetail.split(',') # split line into a list of elements
      dictlist = [a.strip().split(': ') for a in linelist] # map the line to a list of lists for dictionary conversion
      linedict = dict(dictlist) # convert the list of lists into a dictionary
      #self.src = int(linedict['current'].strip())
      self.curr_cell = int(linedict['current'].strip())
      self.new_cell = int(linedict['new'].strip())
      self.step = int(linedict['stepnum'].strip())

class Transition():
   """ Object representing a transition event between two milestones. 

   Used to construct transition statistics for analysis.

   """
   def __init__(self, line):
      line = line.strip() # we have to parse the line
      self.line = line
      linetail = line[len(GRID_TRANSITION)-1:] # just take the last bit of the line, the important part, but not the endline
      linelist = linetail.split(',') # split line into a list of elements
      dictlist = [a.strip().split(': ') for a in linelist] # map the line to a list of lists for dictionary conversion
      linedict = dict(dictlist) # convert the list of lists into a dictionary
      self.src = int(linedict['source'].strip())
      self.dest = int(linedict['destination'].strip())
      self.cur_step = float(linedict['stepnum'].strip())
      self.time = float(linedict['incubation time'].strip().split()[0])
      self.anchor = int(linedict['anchor'].strip())
      #self.curr_cell = int(linedict['curr_cell'].strip())
      #self.new_cell = int(linedict['new_cell'].strip())
      #self.ligand_com = linedict['ligand COM'].strip().split()
      #self.receptor_com = linedict['receptor COM'].strip().split()
      #self.receptor_start_com = linedict['receptor start COM'].strip().split()
      #self.umbrella_step = int(linedict['ID'].strip().split()[0])
      #self.velocity_step = int(linedict['VEL_ID'].strip().split()[0])

   def print_status(self):
      print("src:", self.src)
      print("dest:", self.dest)
      print("cur_step:", self.cur_step)
      print("time:", self.time)
      #print "ligand_com:", self.ligand_com
      #print "receptor_com:", self.receptor_com
      #print "receptor_start_com:", self.receptor_start_com

#======================================================
# Functions
#======================================================

def _parse_milestoning_file(milestoning_filename):
   """given a milestoning file, will parse the XML and generate a model object.

   The model object contains all cell anchors
   
   Parameters
   --------
   milestoning_filename : str
      string corresponding to the name of the milestones XML file

   Returns
   ---------
   model : obj
      Model object to contain all milestone and cell transition information

   """
   tree = ET.parse(milestoning_filename) # now the entire XML file is in a dictionary
   root = tree.getroot() # object at bottom of the tree
   model = Model() # create the model object

   # Temperature
   xml_temp = root.find('temperature') # temperature tag
   if xml_temp != None: # make sure it exists
      model.temperature = float(xml_temp.text.strip())
   else: # some old milestones.xml files may not have this tag
      model.temperature = DEFAULT_TEMP
      
  # Rootdir
   xml_rootdir = root.find('rootdir') # temperature tag
   if xml_rootdir!= None: # make sure it exists
      rootdir = xml_rootdir.text.strip()
   else: # some old milestones.xml files may not have this tag
      rootdir = DEFAULT_ROOTDIR
      
  # BDdir
   xml_bddir = root.find('browndye_bin_dir') # temperature tag
   if xml_bddir!= None: # make sure it exists
      bddir = xml_bddir.text.strip()
   else: # some old milestones.xml files may not have this tag
      bddir = DEFAULT_BDDIR

   # MD Time Factor
   xml_md_time_factor = root.find('md_time_factor') # temperature tag
   if xml_md_time_factor != None: # make sure it exists
      model.md_time_factor = float(xml_md_time_factor.text.strip())
   else: # some old milestones.xml files may not have this tag
      model.xml_md_time_factor = DEFAULT_MD_TIME_FACTOR

   # MD Time Factor
   xml_bd_time_factor = root.find('bd_time_factor') # temperature tag
   if xml_bd_time_factor != None: # make sure it exists
      model.bd_time_factor = xml_bd_time_factor.text.strip()
   else: # some old milestones.xml files may not have this tag
      model.xml_bd_time_factor = DEFAULT_BD_TIME_FACTOR


   site_counter = 0 
   milestone_id_list = []
   for branch in root:
      if branch.tag != "site": continue # make sure that we are parsing a site
      site_name = branch.find('name').text.strip()
      print("Now parsing milestones from site:", site_name, "in XML file.")
      site_obj = Site(site_counter, site_name) # create the site object
      for anchor in branch: #iterate through each anchor in the site
         if anchor.tag != "anchor": continue #ensure we are reading voronoi anchors
         index = anchor.find('index').text.strip()
         #coord = anchor.find('coord').text.strip()
         fullname = anchor.find('name')
         if fullname != None:
            fullname = fullname.text.strip() # parse all the attributes of the milestone
         else: # then just name it by its anchor
            fullname = str(anchor)
         directory_text = anchor.find('directory').text
         if directory_text:
            directory = directory_text.strip() # directory where this milestone is located
         else:
            directory = None
         bd =boolean(anchor.find('bd').text.strip())
         md =boolean(anchor.find('md').text.strip())
         anchor_obj = Anchor(index, md, bd, fullname, directory, site_counter, site_name,) 
      

         for milestone in anchor: #read all of the milestones in this anchor
            if milestone.tag != "milestone": continue
            id = milestone.find('id').text
            radius = None
            normal = None
            group_xml = milestone.find('group')
            if group_xml != None:
               group = group_xml.text.strip()
            #else:
            #   shape = DEFAULT_SHAPE # by default
            #if shape == "plane": # get other information based on the milestone shape
               #  anchor = milestone.find('anchor').text.strip()
            #   normal = milestone.find('normal').text.strip()
            #elif shape == "sphere":
               #  anchor = milestone.find('anchor').text.strip()
            #   radius = milestone.find('radius').text.strip()
            #elif shape == "rotational":
            #   pass
            end_xml = milestone.find('end')
            if end_xml != None:
              end = milestone.find('end').text.strip()
            else:
              end = "false"
            value_xml = milestone.find("value")
            value = value_xml.text.strip()
            milestone_obj = Milestone(id, group, value, end)
            if milestone_obj.id not in milestone_id_list: milestone_id_list.append(milestone_obj.id)
            anchor_obj._add_milestone(milestone_obj)
         site_obj._add_anchor(anchor_obj)
     
      model._add_site(site_obj)
      site_counter += 1 


      bd_index = len(milestone_id_list) -1
      b_surf_index = len(milestone_id_list)
      model.b_surface = Anchor(index=b_surf_index,  md=False, bd=True, 
                               fullname="b_surface", 
                               directory=os.path.join(rootdir, "b_surface"), 
                               siteindex=0, sitename="b_surface")

      #print("bd index", bd_index)
      model.bd_milestone = Anchor(index=bd_index, md=False, bd=True, 
        fullname="bd_milestone", 
        directory=os.path.join(rootdir, "bd_milestone"), 
        siteindex=0, sitename="bd_milestone") 
      
      model.browndye_bin_dir = bddir
     

   return model 

def add_dictionaries(dict1, dict2):
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

def _parse_bd_results(bd_results_filename):
  ''' given a BD results file name, will open the file and extract information 
  about state transitions

    Parameters
    --------
    bd_results_filename : str
      string corresponding to the name of the bd results file

    Returns
    ---------
    bd_dict: dictionary
      dictionary containing transition statistics for bd simulations
    bd_time: float
      total BD simulation time (if recorded by Browndye)
    

  '''
  
  bd_dict = {}
  bd_time = None
  tree = ET.parse(bd_results_filename)
  root = tree.getroot()
  for tag in root:
    if tag.tag == "reactions":
      reactions = tag
      for tag2 in reactions:
        i = 0
        if tag2.tag == "escaped":
          bd_dict['inf'] = int(tag2.text)
        elif tag2.tag == "completed":
          text = tag2[0].text.strip()
          milestone = text.split('_')[1] #get just the milestone number
          #print "site:", site
          n = tag2[1].text
          #print "n:", n
          #name = outer_state[i] + '_' + str(site)
          bd_dict[milestone] = int(n)
          i += 1
        elif tag2.tag == "time":
          bd_time = float(tag2.text)

  return bd_dict, bd_time

# def parse_bound_state_args(bound_args):
#    """Parses a user-defined string corresponding to all bound states

#    Parameters
#    -----------
#    bound_args : str
#       string of bound state values
#       can use ',' or ':' for a range

#    Returns
#    ----------
#    bound_dict : dict
#       dictionary of bound state indices

#    """

#    bound_dict = {}
#    bound_pairs = bound_args.split(',')
#    for pair in bound_pairs:
#       # print 'PAIR'+ pair
#       site_index = pair.split(':')
#       # print site_index
#       if len(site_index) == 1:
#          site = 'all'
#          index = site_index[0]
#       elif len(site_index) == 2:
#          site = site_index[0]
#          index = site_index[1]
#       if site not in bound_dict:
#          bound_dict[site] = [index]
#       else:
#          bound_dict[site].append(index)
#          # print bound_dict
#    return bound_dict


def _read_transition_statistics_from_files(model, verbose):
   """Parses the transitions statistics from the simulation output files for later analysis

   Parameters
   -------
   model : obj
      object containing all anchor and milestone information

   Returns
   ---------

   total_steps : int
      total number of MD steps taken in all simulations in all anchors

   """
   total_steps = 0
   for site in model.sites:
      for anchor in site.anchors:
         if anchor.md == True and anchor.directory:
            #if verbose: print 'parsing md transitions for:Anchor', milestone.fullname
            #print info['max_steps']
            print('parsing md transitions for:Anchor', anchor.fullname)
            max_steps = anchor._parse_md_transitions()
            print(max_steps, total_steps)
            if max_steps > total_steps:
              total_steps = max_steps
  
   return total_steps


def make_model( milestone_filename="milestones.xml", verbose=False):
   """ Top level function used to parse inputs and create the milestoning model
   
   Parameters
   ------
   milestone_filename : str , optional
      name of the XML file containing milestone and anchor information
   bound_states : str, optional
      string of indices of all bound states for the model

   Returns
   ---------
   model: obj
      the complete milestoning model with extracted transition statistics
   bound_dict : dict
      dictionary of bound state indices
   max_steps : int
      maximum/total number of steps from the anchor with the longest simulation time

   """

   #bound_dict = parse_bound_state_args(bound_states)
   model = _parse_milestoning_file(milestone_filename)
   max_steps = _read_transition_statistics_from_files(model, verbose)

   return model, max_steps


def boolean(arg):
  if str(arg).lower() in ("false","", "0", "0.0", "[]", "()", "{}"):
   return False
  else:
   return True



