## gEDES -generalized Ensemble Docking with Enhanced sampling of pocket Shape-

This folder contains all the scripts and the tools to prepare and analyze gEDES simulations, and to perform docking calculations with Autodock.
There are 3 folders numbered from 1 to 3. A brief explanation of the content of each folder is reported below (note that there are several README.txt files within each folder, 
containing more exhaustive information on specific tasks and scripts):


*1_setup_and_run_MD_simulations:* Folder containing the scripts used to setup the input files for AMBER and GROMACS simulations starting from the pdb file of the protein of interest.

*2_clustering_scripts:* Folder containing the scripts to extract cluster conformations (as pdb files) from the COLVAR files and then from the corresponding trajectory. 
Note that in the case of standard MD simulations performed with AMBER, it is sufficient to rerun the trajectory with GROMACS/PLUMED to record and save the values of the CVs in a COLVAR file.

*3_Autodock_scripts:* Autodock scripts to perform docking, relaxation, and rescoring of the poses.

Note that a tutorial explaining the use of a previous version of gEDES is available at the following link: 
- https://molmod.dsf.unica.it/enhanced-sampling-of-binding-pocket-volume-and-shape-tutorial


gEDES paper (Feb, 2025)
- https://doi.org/10.1021/acs.jcim.4c01810
