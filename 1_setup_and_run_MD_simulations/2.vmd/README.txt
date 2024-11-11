#the script "divide_into_list.tcl" takes a binding site and defines 3 pairs of groups of atoms on opposite sides of 3 orthogonal intertia planes.
#the groups of atoms are identified by index numbers that define the groups of atoms of a COORDINATION variable in PLUMED.
#to exectute the script, just load your protein pdb in VMD and load the script from the tkconsole (source divide_into_list.tcl) after editing the list of residues lining the (putative) BS you want to enhance the sampling of.
#
#four files will be printed:
#1. serials_RoG --> serial numbers defining the heavy atoms of the BS, which will enter into the CV GYRATION
#2. serials_xy --> serial numbers defining the two groups of atoms (A - Z1 and B - Z2) entering the definition of the first COORDINATION CV
#3. serials_xz --> serial numbers defining the two groups of atoms (A - Y1 and B - Y2) entering the definition of the second COORDINATION CV
#4. serials_yz --> serial numbers defining the two groups of atoms (A - X1 and B - X2) entering the definition of the third COORDINATION CV
#
#NOTE: there is also a basic description of the usage of this script here: https://molmod.dsf.unica.it/enhanced-sampling-of-binding-pocket-volume-and-shape-tutorial-part-2-4

