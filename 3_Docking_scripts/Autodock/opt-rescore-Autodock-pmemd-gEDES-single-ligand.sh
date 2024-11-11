#!/bin/bash
#load program and modules
source /etc/modules.sh
module load vmd193 amber20
source $AMBERHOME/amber.sh
#set directories and system names
root=..
#root2=${root}/$root
scridir=./accessory_scripts
mgldir=~/local/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24
ligprepdir=./amber-bcc-ligs
LIG=$1
bs=$2
box=$3
grid=$4
#add box when rescoring so that the ligand is inside!
boxplus=15
#
nargs=${#@}
if [ $nargs != 4 ]; then
    echo "usage: ./opt-rescore-Autodock-pmemd.sh LIG (e.g. UDP, NEO, CMG) bs (eg gluco, ricina, cdk2), box of simulation (eg 40), and Autodock grid (0.25)"
    exit
fi

num=$(wc -l list-ordered-energy.tcl | awk '{print $1}')
tail -n $(($num-1)) list-ordered-energy.tcl | awk '{print $1, $2, $3}' > load-all.tcl

mkdir rescore${grid}
cd rescore${grid}
#cycle over all complexes
cat <<EOF>measure-rmsd-lig.tcl
set num [expr {[molinfo num] - 1}]
for {set id 0} {\$id < \$num} {incr id} {
    set name [molinfo \$id get name]
    set ref [atomselect top "noh not protein"]
    set cpd [atomselect \$id "noh not protein"]
    set all [atomselect \$id all]
    set trans_mat [measure fit \$cpd \$ref]
    \$all move \$trans_mat
    set rmsd [measure rmsd \$cpd \$ref]
    puts \$rmsd
    #if { \$rmsd < 0.5} { 
    set out [open "measure-rmsd-lig.dat" w]
    puts \$out "\$name \$rmsd"
    close \$out
    quit
    #}
}
quit
EOF
#select the appropriate prep file
for i in $(cat ../load-all.tcl | awk '{print $NF}')
do
    j=$(basename $i .pdb)
    mkdir $j
    echo $j
    for kk in $(ls ${root}/${root}/${ligprepdir}/${LIG}*.mol2)
    do
	vmd -dispdev none -m ${kk} ${root}/${root}/${i} -e measure-rmsd-lig.tcl
    done
    #
    mv measure-rmsd-lig.dat $j
    cd $j
    lj=$(basename $(cat measure-rmsd-lig.dat | awk '{print $1}') .mol2)
    ln -s ${root}/${root}/${i} .
    ln -s ${root}/${root}/${scridir}/*tcl .
    ln -s ${root}/${root}/${root}/${scridir}/opt-restr?.pmemd .
    ln -s ${root}/${root}/AD4.1_bound.dat .
    #copy prep and frcmod files
    ln -s ${root}/${root}/${ligprepdir}/${lj}.prep .
    ln -s ${root}/${root}/${ligprepdir}/${LIG}.frcmod .
    vmd -dispdev none ${root}/${root}/$i -e save-noh.tcl
    #ligands are named ${LIG} in the pdbs...
    awk 'BEGIN{a=0}{if (($1 == "ATOM") && ($4 !~ /'${LIG}'/)) {if (a==0) {print "TER"}; print $0; a=1} else {print $0}}' tmp.pdb | cut -c 1-54 > ${j}_noh.pdb
    #generate topology
    cat<< EOF > create-top.leap
loadamberprep ${lj}.prep
loadamberparams ${LIG}.frcmod
com=loadpdb ${j}_noh.pdb
setbox com vdw ${box}
saveamberparm com com.parm7 com.rst7
quit
EOF
    tleap -s -f leaprc.protein.ff14SB -f leaprc.gaff -f create-top.leap
    ncyc=1
    for kk in $(seq 1 $ncyc)
    do
	source ./opt-restr${kk}.pmemd
	#$LIG
    done
    ~/scripts/cpptraj/nc_to_pdb.cpptrj com.parm7 rest1.rst7 ${j}_opt1.pdb
    rm cpptraj.in
    cd ..
    ln -s ${j}/${j}_opt1.pdb .
done
    #
    #now rescore...
    #
for i in *_opt1.pdb
do
    ii=$(basename $i .pdb)
    #prefix=$(echo ${ii} | awk 'BEGIN{FS="_";OFS="_"}{print $1,$2,$3}')
    #prefix=$(echo ${ii} | awk 'BEGIN{FS="_";OFS="_"}{print $1,$2}')
    awk -f ./accessory_scripts/reorder_pdb_auto-with-TERs.awk ${i} | grep -v TER > tmp2.pdb ; mv tmp2.pdb ${i} 
    #name=$(echo ${i} | awk 'BEGIN{FS="_-"}{print $1}')
    #ligands are named ${LIG} in the pdbs...
    grep "${LIG}" ${i} > ${ii}_${LIG}.pdb
    ${mgldir}/prepare_ligand4.py -l ${ii}_${LIG}.pdb -o ${ii}_${LIG}.pdbqt -U nphs -Z
    #ligands are named ${LIG} in the pdbs...
    grep -v "${LIG}" $i > ${ii}_rec.pdb
    ${mgldir}/prepare_receptor4.py -r ${ii}_rec.pdb -o ${ii}_rec.pdbqt -U nphs
    sed -e "s/XXX/${bs}/" -e "s/GGG/${grid}/" ${root}/${root}/${scridir}/bs_dim_center.tcl > bs_dim_center.tcl
    cut -c 1-54 ${ii}_rec.pdbqt > tmp.pdb
    vmd -dispdev none tmp.pdb -e bs_dim_center.tcl
    pts=$(head -n 1 bs_dimension.dat | awk 'BEGIN{FS=",";OFS=","}{print $1+'${boxplus}',$2+'${boxplus}',$3+'${boxplus}'}')
    cnt=$(tail -n 1 bs_dimension.dat)
    ${mgldir}/prepare_gpf4.py -l ${ii}_${LIG}.pdbqt -r ${ii}_rec.pdbqt -p npts="${pts}" -p spacing="${grid}" -o tmp.gpf -y
    sed 's/#/ #/g' tmp.gpf > ${ii}-rescore.gpf
    autogrid4 -p ${ii}-rescore.gpf -l ${ii}-rescore.glg
    ${mgldir}/prepare_dpf4.py -l ${ii}_${LIG}.pdbqt -r ${ii}_rec.pdbqt -o tmp.dpf
    awk '{if ($1 ~ /^about/) {print "epdb" ;exit} else print $0}' tmp.dpf | sed 's/#/ #/g' > ${ii}-rescore.dpf
    autodock4 -p ${ii}-rescore.dpf -l ${ii}-rescore.dlg
    rm *.*map *.*fld *.*xyz
    newaff=$(grep "Estimated Free Energy" ${ii}-rescore.dlg | awk '{print $9}')
    mv ${i} ${ii}-rescored_${newaff}.pdb
    rm ${ii}_rec.pdbqt ${ii}_rec.pdb ${ii}_${LIG}.pdb ${ii}_${LIG}.pdbqt ${ii}-rescore.glg ${ii}-rescore.gpf 
    rm bs_dimension.dat tmp.pdb tmp.dpf tmp.gpf
done
