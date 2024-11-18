#################################
# Part 0: loading the trajectory
#################################

set outfile_ROG [open "RoG_trajectory.dat" w]
set outfile_RMSD [open "RMSD_fromholo.dat" w]

# reference structure, molid: 0
mol new ../holo_processed.pdb type pdb
set mol [mol new frames.pdb type pdb waitfor all]
# set mol [mol new meta_clusters.pdb type pdb waitfor all]

#expr bs
set reference [atomselect 0 "protein and noh and resid 8 9 10 11 12 13 15 31 32 35 36 57 58 59 64 84 85 86 88 92 119 120 123 134 137 156 158 167 200 201 202 205"]
set compare [atomselect 1 "[$reference text]"]
set pdbselection [atomselect 1 "protein"]

# If one want to use frame 0 of the trajectory as the reference;
# set reference [atomselect 0 "protein" frame 0]
# the frame being compared
# set compare [atomselect 0 "protein"]

# Loop performing the work for each frame of the trajectory; alignment is performed, the RMSD from the reference structure is calculated, and the corresponding frame is saved.

                set num_steps [molinfo $mol get numframes]
                for {set frame 0} {$frame < $num_steps} {incr frame} {
                        $compare frame $frame
                        animate goto $frame
                        set trans_mat [measure fit $compare $reference]
                        $pdbselection move $trans_mat
		        animate write pdb ./aligned_$frame.pdb beg $frame end $frame sel [atomselect top "protein"] 
                        set rmsd [measure rmsd $compare $reference]
		        set rog [measure rgyr $compare]
			puts $outfile_ROG  "$frame $rog"
			puts $outfile_RMSD "$frame $rmsd"
               }
       close $outfile_ROG
       close $outfile_RMSD
quit
