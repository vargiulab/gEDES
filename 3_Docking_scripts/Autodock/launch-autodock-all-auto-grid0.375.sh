#!/bin/bash
#load program and modules
source /etc/modules.sh
module load vmd193 amber20
source ${AMBERHOME}/amber.sh
#set directories and system names
LIG=$1
p=$2
bs=$3
boxplus=0
boxscale=1.0
grid=0.375
nargs=${#@}
if [ $nargs != 3 ]; then
    echo "usage: ./ LIG (e.g. UDP, NEO, CMG) protein (e.g. GLUCO, RICIN, CDK2), bs (eg gluco, ricina, cdk2)"
    exit
fi
root=..
scridir=./accessory_scripts
mgldir=/home/vargiu/local/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24
dirconf=.
n=$(($(ls ${dirconf}/${p}*pdbqt | wc -l)-1))
parfile=AD4.1_bound.dat
ln -s ${root}/${LIG}.pdbqt .
ln -s ${root}/${parfile} .
#set the BS (defined as macro in .vmdrc) and the grid
sed -e "s/XXX/${bs}/" -e "s/GGG/${grid}/" ${scridir}/bs_dim_center.tcl > bs_dim_center.tcl
#run docking
for i in $(seq 0 $n)
do
    cut -c 1-54 ${p}_${i}.pdbqt > tmp.pdb
    vmd -dispdev none tmp.pdb -e bs_dim_center.tcl
    #pts=$(head -n 1 bs_dimension.dat)
    #pts=$(head -n 1 bs_dimension.dat | awk 'BEGIN{FS=",";OFS=","}{print $1+'${boxplus}',$2+'${boxplus}',$3+'${boxplus}'}')
    pts=$(head -n 1 bs_dimension.dat | awk 'BEGIN{FS=",";OFS=","}{printf "%d,%d,%d\n", $1*'${boxscale}',$2*'${boxscale}',$3*'${boxscale}'}')
    cnt=$(tail -n 1 bs_dimension.dat)
    ${mgldir}/prepare_gpf4.py -l ${LIG}.pdbqt -r ${p}_${i}.pdbqt -p npts="${pts}" -p gridcenter="${cnt}" -p spacing="${grid}" -o ${p}_${i}.gpf
    autogrid4 -p ${p}_${i}.gpf -l ${p}_${i}.glg
    ${mgldir}/prepare_dpf4.py -l ${LIG}.pdbqt -r ${p}_${i}.pdbqt -p ga_num_evals=25000000 -p ga_run=1 -p parameter_file=${parfile} -o tmp.dpf
    sed 's/#/ #/g' tmp.dpf > ${p}_${i}.dpf
    #sed "s/holo/${p}_${i}/g" ../${p}-template.dpf > ${p}_${i}.dpf
    autodock4 -p ${p}_${i}.dpf -l ${p}_${i}.dlg
    rm ${p}_${i}.*map ${p}_${i}.*fld ${p}_${i}.*xyz
done
rm bs_dimension.dat tmp.pdb
