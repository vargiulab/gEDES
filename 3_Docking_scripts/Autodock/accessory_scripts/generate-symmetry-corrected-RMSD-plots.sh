#!/bin/bash
source /etc/modules.sh
module load amber20
source $AMBERHOME/amber.sh

p=$1 #e.g. AK
l=$2 #e.g. ADP
bs=$3 #e.g. adp_exp
TOP=${ligprepdir}/holo_complex_${p}-amber.parm7
REF=${ligprepdir}/holo_complex_${p}-amber.rst7
for i in $(ls compare*vmd)
do
    j=$(basename $i .vmd | sed 's/compare_Xray-//')
    vmd -dispdev none -e $i <<EOF
animate write dcd traj-${j}.dcd waitfor all sel [atomselect top "all"]
EOF
    #BSCPPTRJ=$(grep "_${bs} " ~/.vmdrc | grep -v \# | awk '{for (i=1;i<=NF;i++) {if ($i=="resid") {count=1} ; if (count==1) {printf"%d,",$i+1}}}END{print ""}' | sed -e 's/resid\,/:/' -e 's/)}\,//')
    #this is needed because protein residues in docking complexes start from 2...
    BSCPPTRJ=$(grep "_${bs} " ~/.vmdrc | grep -v \# | awk '{for (i=1;i<=NF;i++) {if ($i=="resid") {count=1;count2=0;printf":"} ; if (count==1) {count2++; if (count2>1) {printf"%d,",$i+1}}}}END{print ""}' | sed -e 's/,$//')
    cat<<EOF>symm-corr-RMSD-${j}.cpptrj
parm $TOP
trajin traj-${j}.dcd
reference $REF [REF]
symmrmsd RMSDPROT @CA,C,O,N @CA,C,O,N out RMSD_prot-${j}.rmsd remap ref [REF]
symmrmsd RMSDBS (${BSCPPTRJ})&!@H= (${BSCPPTRJ})&!@H= out RMSD_bs-${j}.rmsd remap ref [REF]
symmrmsd RMSDLIG :${l}&!@H= out RMSD_lig-${j}.rmsd nofit ref [REF]
EOF
    cpptraj.cuda -i symm-corr-RMSD-${j}.cpptrj
    paste RMSD_bs-${j}.rmsd RMSD_lig-${j}.rmsd | awk '{if (NR==1) {print "#BS LIG"} else {print $2, $4}}' > scatter-bs-lig-${j}.rmsd
    head -n 11 scatter-bs-lig-${j}.rmsd > scatter-bs-lig-${j}-top10.rmsd
done
