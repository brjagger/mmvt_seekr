# MMVT SEEKR Tutorial

SEEKR is an open source package for calculation protein-ligand binding kinetics using MD, BD, and milestoning. This tutorial will demonstrate the basics needed to setup, run, and analyze a SEEKR calculation.


## System Preapration

The preparation portion SEEKR is handled by the python script 'seekr.py', which reads all the necessary setup configurations from an input file. For this tutorial, a completed input file 'test.seekr' is provided in the 'tests/tryp_inp_test' directory. 

Take a look at this file using your favorite text editor.

The input file is divided into various sub-sections with headers describing the parameters being set in each section.

The **General System Info** section lets us set the 'project_name'. The 'rootdir' should be set to the path where you would like to create the SEEKR filetree.

The **System Files** section specifies the path to each of the coordinate/parameter files that SEEKR will use to set up the calculation. For this tutorial, the paths all correspond to files in the 'mmvt_seekr/data' directory

The **Milestone Definition** block is where the user can specify the collective variable that will be discretized by milestones and all necessary parameters. For this example, a single distance based collective variable, 'milestone1' is defined. The 'colvar_type' is set to distance, and the two 'group' variables contain lists of atom indices, the distance will be measured between the centers of mass of these two atom groups.

The 'milestone_values' variable allows the user to specify at which collective variable values a milestone will be places (and the corresponding simulation files generated). Finally the 'equil_rest_force' variable is used to set the strength of a harmonic restraint (in kcal/mol) that pulls the system to the corresponding milestone value. (This approach is not recommended for all systems) 

The **Anchor Equilibration** section contains parameters set in the NAMD input files for the simulations that pull the ligand to the desired milestone value. This includes the simulation ensemble, length of simulation, and trajectory write frequency. These simulation files will only be generated if the 'equil' parameter is set to 'true'.

Finally, the **Production** block contains parameters needed to generate the NAMD simulation files for the MMVT production runs such as the simulation ensemble ('prod_ensemble') and the length of the simulation ('prod_steps').

The SEEKR filetree is created by running the following on the command line from inside the 'tryp_inp_test' directory
```bash
python ../../seekr.py test.seekr
```

Once this runs, you should find a new filetree in the location you specified with the 'rootdir' variable in the .seekr input file. Navigate to this location.

Inside the top-level directory, you should finde multiple directories named 'anchor_#'. Each of these corresponds to one simulation cell whose edges are milestones we specified in the input file. You will also see a 'milestones.xml' file that containes information about each of the milestones and will be useful for analysis later.

Navigate into one of the anchor directories and you will notice another directory called 'md'. Some anchors may also have a 'bd' directory.
Inside the 'md' directory are three more directories:

 'building' contains the coordinates and force field parameter files for the system

 'equil' contains a NAMD input file and colective variables module script used to run simulations to pull the ligand to the appropriate milestone distance

 'prod' contains the NAMD input file 'prod_1.namd' used to run the MMVT simulations. If you look at this file, you will notice some code at the bottom of the file which facilitates the boundary checking and velocity reversal criteria needed for MMVT SEEKR.


## Running Simulations

Note: MMVT simulations are relatively expensive, meaning you probably do not want to run the complete set for this tutorial. Sample outputs can be found in the './mmvt_seekr/data' directory

Equilibration (if used) and production simulations can both be run using standard NAMD command line executables.

## Analysis

The analysis portion of SEEKR is designed to be conducted in a jupyter notebook by importing the necessary SEEKR libraries, 'analyze.py', 'model.py', and 'plots.py'.

A sample jupyter notebook 'tutorial.ipynb' can be found in the 'data' directory and walks through the analysis portion of the tutorial.



