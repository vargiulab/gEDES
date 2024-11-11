#!/bin/bash

# Run gmx_mpi command here
mpirun -np 8 gmx_mpi mdrun -gpu_id 0123 -multidir rep0 rep1 rep2 rep3 rep4 rep5 rep6 rep7 -s topol -deffnm BE_AK -cpi BE_AK.cpt -plumed plumed-common.dat -dlb yes -replex 50000 -nsteps 275000000 

# Check the exit status of gmx_mpi in a loop
while true
do
    if ! ps ax | grep -v grep | grep gmx_mpi > /dev/null
    then
        # Launch relaunch-meta.sh when gmx_mpi stops
	mpirun -np 8 --hostfile host gmx_mpi mdrun -gpu_id 0123 -multidir rep0 rep1 rep2 rep3 rep4 rep5 rep6 rep7 -s topol -deffnm BE_AK -append -cpi BE_AK.cpt -plumed plumed-common.dat -dlb yes -replex 50000 -nsteps 275000000
       # break
    fi
    sleep 30
done

