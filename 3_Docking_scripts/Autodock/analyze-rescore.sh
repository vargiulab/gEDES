#!/bin/bash
p=$1
l=$2
bs=$3

root=..
scridir=./accessory_scripts
ligprepdir=./amber-bcc-ligs

#cycle over trajectories of rescored poses
for i in $(ls -d */)
do
    cd ${i}/rescore
    cp ${scridir}/generate-rescore-tcl.sh ${scridir}/generate-symmetry-corrected-RMSD-plots.sh ${scridir}/fit-and-measure-rescore.tcl ${scridir}/compare_Xray*vmd .
    ./generate-rescore-tcl.sh && ./generate-symmetry-corrected-RMSD-plots.sh ${p} ${l} ${bs}
    cd ../..
done
