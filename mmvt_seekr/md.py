def main(settings):
  rootdir = settings['rootdir']
  wet_configs=settings['wet_configs']
  raw_milestone_list = settings['raw_milestone_list']
  milestone_pos_rot_list = settings['milestone_pos_rot_list']
  #dry_configs=settings['dry_configs']

  config_dirlist = []
  md_file_paths = []
  bd_file_paths = []
  config_counter = 0
  
  for i in range(len(wet_configs)): # for each configuration, it gets its own directory
    wet_config = wet_configs[i]
    anchor_name = "anchor_%s" % (wet_config.struct_id,)
    anchor_filetree = Filetree({anchor_name:{}})
    if verbose: print "Creating directory:", anchor_name
    anchor_filetree.make_tree(rootdir) # create this anchor's directory
    milestone_pos_rot_list[i][0].directory = anchor_name # update this milestones directory information
    anchor_dir = os.path.join(rootdir, anchor_name)
    # MD filetree
    md_file_path={}
    if settings['md'] == True and milestone_pos_rot_list[i][0].md == True: # then prep this anchor for an MD simulation
      md_dir=os.path.join(anchor_dir,'md') # directory for MD
      md_filetree=Filetree({'md':mdtree})
      md_filetree.make_tree(anchor_dir) # make the MD filetree
      for key in mdtree.keys():
        md_file_path[key] = os.path.join(md_dir,key)
      wet_holo_filename=os.path.join(md_dir,'holo_wet.pdb')
      wet_config.save(wet_holo_filename, amber=True, standard=False) # write the holo structure into the md directory
      md_file_path['wet_holo'] = wet_holo_filename
    # BD filetree
    bd_file_path={}
    bd_dir = None
    if settings['bd'] == True and milestone_pos_rot_list[i][0].bd == True: # then prep this anchor for BD simulation
      dry_config= wet_config #pdb.dry(wet_config) # since we didn't iterate through this list, we must use the index
      bd_dir=os.path.join(anchor_dir,'bd') # directory for MD
      bd_filetree=Filetree({'bd':bdtree})
      bd_filetree.make_tree(anchor_dir) # make the MD filetree
      for key in bdtree.keys():
        bd_file_path[key] = os.path.join(bd_dir,key)
      dry_holo_filename=os.path.join(bd_dir,'holo_dry.pdb')
      dry_config.save(dry_holo_filename, amber=True, standard=False) # write the holo structure into the md directory
      md_file_path['dry_holo'] = dry_holo_filename


    config_dirlist.append(anchor_name)
    if md_file_path: md_file_paths.append(md_file_path)
    if bd_dir: bd_file_paths.append(bd_dir) # bd_file_path
    config_counter+=1 # increment the loop counter


  return config_dirlist, md_file_paths, bd_file_paths, raw_milestone_list