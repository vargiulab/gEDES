For clustering, one need to have cumulative trajectory in any format (dcd is preferred), and the holo parameters and structure (parm7 and rst7 formats, respectively) for the rmsd calculation using cpptraj (this is just for comparison purposes).

You need to run the scripts below in this order. Please modify them to define the names/paths for your system: 
1. ./gen_colvar_meta.sh 
2. ./R_clustering_8CVs.sh

