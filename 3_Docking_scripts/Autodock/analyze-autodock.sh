#!/bin/sh
#script to analyze the poses and generate a summary table.
#script directory
scridir=./accessory_scripts
#parameters
#only average the first "avcut" poses within each cluster (HADDOCK-style)
#topcut=1
#avcut=4
#cutoff for final analysis (last command of this script)
#cutstat=2
#minimum cutoff for ligand clustering
#cutmin=1.5
#usage
nargs=${#@}
if [ $nargs != 3 ]; then
    echo "usage: ./analyze-full.sh LIG (e.g. UDP, NEO, CMG) protein (e.g. GLUCO, RICIN, CDK2) bs (eg gluco, ricina, cdk2)"
    exit
fi
#create folders
if [ ! -e frames-ali ]
    then
	mkdir frames-ali
fi
if [ ! -e clusters ]
    then
	mkdir clusters
	mkdir clusters/ali-bm
	mkdir clusters/ali-bs
fi
#
#load modules and set the names of the variables
module load amber20 && source $AMBERHOME/amber.sh
l=$1
p=$2
bs=$3
##0. copy the files if missing
#for i in fit-bs-and-deviation-multi-cutoff.tcl fit-bm-and-deviation-multi-cutoff.tcl fit-bm-and-deviation-traj.tcl fit-bs-and-deviation-traj.tcl
#do
#    if [ ! -f $i ]
#    then
#	ln -s ${scridir}/$i .
#    fi
#done

#1. first generate the complexes
n=$((1+$(grep ATOM ../${l}.pdbqt | wc -l)))
#n=$(awk '{if (($1=="ATOM") && ($NF !~ /H/)) print $0}' ../${l}.pdbqt | wc -l)
for i in $(ls *.dlg)
do
    pf=$(basename $i .dlg)
    e=$(grep -m 1 "DOCKED: USER    Estimated Free Energy of Binding" $i | awk '{print $9}')
    grep -m $((${n}+1)) -e "DOCKED: USER    Estimated Free Energy of Binding" -e "DOCKED: ATOM" -e "DOCKED: ENDROOT" $i \
	| sed -e 's/DOCKED: USER/REMARK/g' -e 's/DOCKED: ATOM/ATOM/g' -e 's/DOCKED: \ENDROOT/TER/g' -e 's/ENDROOT/TER/g' -e "s/${l} \{3,5\}[0-9]\{1,3\}/${l}     0/g" > l.pdb
    cat l.pdb $pf.pdbqt | awk '{if ($1 ~/^ATOM/) {print substr($0,0,54)} else {print $0}}' > ${pf}_complex_${e}.pdb
done
rm l.pdb
echo "mol new ../holo_complex_${p}.pdbqt type pdb" > list-ordered-energy.tcl
ls -l *complex*.pdb | sort -t "_" -gk 4 | awk '{if (NR==1) {printf "mol new     %s type pdb\n", $NF} else {printf "mol addfile %s type pdb\n", $NF}}' >> list-ordered-energy.tcl
numcom=$(ls -l *complex*.pdb | wc -l)
##2. then perform cluster analysis (and generate the trajectory ordered by energy)
#TOP=$(ls ${p}_0_complex_*.pdb)
#traj=traj-ordered-energy.dcd
#tail -n ${numcom} list-ordered-energy.tcl > load-and-save-traj-ordered.energy.tcl
##
#tail -n ${numcom} list-ordered-energy.tcl > load-fit-bm-and-save-traj-ordered.energy.tcl
#sed -e "s/XXX/${l}/"\
#    -e "s/YYY/${p}/"\
#    -e "s/ZZZ/${bs}/" ${scridir}/fit-bm-traj-first-frame-ori.tcl \
#    | awk '{if ($1 !~ /^#/) {print $0}}'\
#    > fit-bm-traj-first-frame.tcl
#
#cat<<EOF>>load-fit-bm-and-save-traj-ordered.energy.tcl
#source fit-bm-traj-first-frame.tcl
#animate write dcd "${traj}" waitfor all sel [atomselect top "all"]
#quit
#EOF
#vmd -dispdev text -e load-fit-bm-and-save-traj-ordered.energy.tcl
##clustering with cutoff defined by cutclus
#cutclus=$(echo "${n}*0.05" | bc -l)
#echo "#"
#echo "the cutoff for clustering is $n * 0.075 = ${cutclus}, or ${cutmin} Angstrom if that number is lower"
#echo "#"
#if [ $(echo ${cutclus}'>'${cutmin} | bc -l) -eq 0 ]
#then
#    cutclus=${cutmin}
#fi
#cat > cpptraj.in <<EOF
#parm $TOP [top]
#trajin $traj parm [top]
##cluster :${l}&!@H= hieragglo clusters 10 complete dme[:${l}&!@H=] nofit out cluster-num-vs-time.dat summary cluster-summary.out info cluster-info.out repout cl repfmt pdb clusterout cl-traj clusterfmt pdb
#cluster :${l}&!@H= hieragglo epsilon ${cutclus} complete dme[:${l}&!@H=] nofit out cluster-num-vs-time.dat summary cluster-summary.out info cluster-info.out repout cl repfmt pdb clusterout cl-traj clusterfmt pdb
#EOF
#cpptraj -i cpptraj.in
#
##3. now generate the files with the energies of the cluster structures within each cluster
#start=0
##end=$(($(wc -l list-ordered-energy.tcl | awk '{print $1}' )-1))
#end=$(sort -gk 2 cluster-num-vs-time.dat | tail -n 1 | awk '{print $2}')
#
#awk 'BEGIN{FS="_"}{print $4}' load-and-save-traj-ordered.energy.tcl | sed 's/.pdb//g' | grep -v ^$ | awk '{print NR,$1}' > topenergy_numerow.dat
#
#for ((i=start; i<=end; i++))
#do
#    grep " ${i}\$" cluster-num-vs-time.dat > cluster${i}.dat
#    awk '{print $1}' cluster${i}.dat > cluster${i}_index.dat
#    awk '{printf"\\b%s\\b\n",$1}' cluster${i}_index.dat > cluster${i}_index_ready.dat
#    
#    while IFS='' read -r line || [[ -n "$line" ]]
#    do
#	#    echo "Text read from file: $line"
#	grep -m 1 "${line} " topenergy_numerow.dat >> cluster${i}_energies.dat
#    done < cluster${i}_index_ready.dat
#done
#
##4.then calculate distortion with respect to Xray
##
##FIRST ALIGNING THE BM AND CALCULATIN THE RMSD OF BS, LIG, AND BM
##
##4.1. FOR EACH CLUSTER SEPARATELY
#ls -1 cl-traj.c* > cluster_listed
##num_lines=$(wc -l cluster_listed | cut -d ' ' -f1)
#
##for i in $(seq 0 $num_lines)
#for i in $(seq 0 $end)
#do
#    sed -e "s/CLUSTFILE/cl-traj.c${i}/"\
#	-e "s/XXX/${l}/"\
#	-e "s/YYY/${p}/"\
#	-e "s/ZZZ/${bs}/"\
#	-e "s/HOLOX/\.\.\/holo\_complex\_${p}\.pdbqt/"\
#	-e "s/rmsdlig/rmsdlig${i}/" \
#	-e "s/rmsdbs/rmsdbs${i}/" \
#	-e "s/rmsdbmode/rmsdbmode${i}/" \
#	-e "s/POPOL/popol${i}/" ${scridir}/fit-bm-and-deviation-ori.tcl \
#	| awk '{if ($1 !~ /^#/) {print $0}}'\
#	> fit-and-deviation${i}.tcl
#    vmd -dispdev text -e fit-and-deviation${i}.tcl
#done
#
##4.2. FOR THE TRAJ MADE BY CLUSTERS
#sed -e "s/XXX/${l}/"\
#    -e "s/YYY/${p}/"\
#    -e "s/ZZZ/${bs}/" ${scridir}/fit-bm-and-deviation-traj-ori.tcl \
#    | awk '{if ($1 !~ /^#/) {print $0}}'\
#    > fit-bm-and-deviation-traj.tcl
#vmd -dispdev text -e list-ordered-energy.tcl <<EOF
#source fit-bm-and-deviation-traj.tcl
#EOF
#
##5. now perform statistics and generate the table
########
##topcut
#for i in $(seq 0 $end)
#do
#    awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++}} END {print sum/i}' rmsdbs${i}.dat > avbs${i}-top.dat
#    awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++}} END {print sum/i}' rmsdlig${i}.dat > avlig${i}-top.dat
#    awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++}} END {print sum/i}' rmsdbmode${i}.dat > avbmode${i}-top.dat
#    awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' rmsdbs${i}.dat > sdbs${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' rmsdlig${i}.dat > sdlig${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' rmsdbmode${i}.dat > sdbmode${i}.dat
#done
##
#rm Av_BS-top.dat Av_LIG-top.dat Av_BM-top.dat Av_BS_st.dat Av_LIG_st.dat Av_BM_st.dat Population.dat
#for i in $(seq 0 $end)
#do
#    cat avbs${i}-top.dat >> Av_BS-top.dat
#    cat avlig${i}-top.dat >> Av_LIG-top.dat
#    cat avbmode${i}-top.dat >> Av_BM-top.dat
#    cat sdbs${i}.dat >> Av_BS_st.dat
#    cat sdlig${i}.dat >> Av_LIG_st.dat
#    cat sdbmode${i}.dat >> Av_BM_st.dat
#    cat popol${i}.dat >> Population.dat
#done
##
#rm top_energies_clusters.dat Av_EN-top.dat ST_EN.dat
#for i in $(seq 0 $end)
#do
#    sort -g -k2 cluster${i}_energies.dat | awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++}} END {print sum/i}' > aven${i}-top.dat 
#    sort -g -k2 cluster${i}_energies.dat | awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' > sden${i}.dat
#    cat aven${i}-top.dat >> Av_EN-top.dat
#    cat sden${i}.dat >> ST_EN.dat
#    sort -g -k2 cluster${i}_energies.dat |awk 'NR==1{print $2}' >> top_energies_clusters.dat
#done
##
#paste Population.dat Av_BS-top.dat Av_BS_st.dat Av_LIG-top.dat Av_LIG_st.dat Av_BM-top.dat Av_BM_st.dat Av_EN-top.dat ST_EN.dat top_energies_clusters.dat > ok.dat
#sed 's/,/./g' ok.dat > ok1.dat
#sed 's/[:[:blank:]:]/;/g' ok1.dat > summary-commas.dat
#echo "#Pop      Av_BS Av_BS_st Av_LIG Av_LIG_st Av_BM Av_BM_st Av_EN ST_EN top_energies_clusters" > summary.dat
#sed 's/;/ /g' summary-commas.dat >> summary.dat
#awk 'BEGIN{print"#Pop    <BS>    BSsd    <LIG>   LIGsd   <BM>    BMsd    <EN>    ENsd    ENtop"}{if (NR>1) {printf "%3d\t", $1; for (i=2; i<=NF; i++) {printf"%.2f\t", $i}; print ""}}' summary.dat > summary-tabs.dat
#mv summary.dat summary-ali-bm-top_clusters.dat
#mv summary-tabs.dat summary-tabs-ali-bm-top_clusters.dat
#mv summary-commas.dat summary-commas-ali-bm-top_clusters.dat
#######
##avcut
#for i in $(seq 0 $end)
#do
#    awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++}} END {print sum/i}' rmsdbs${i}.dat > avbs${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++}} END {print sum/i}' rmsdlig${i}.dat > avlig${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++}} END {print sum/i}' rmsdbmode${i}.dat > avbmode${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' rmsdbs${i}.dat > sdbs${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' rmsdlig${i}.dat > sdlig${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' rmsdbmode${i}.dat > sdbmode${i}.dat
#done
##
#rm Av_BS.dat Av_LIG.dat Av_BM.dat Av_BS_st.dat Av_LIG_st.dat Av_BM_st.dat Population.dat
#for i in $(seq 0 $end)
#do
#    cat avbs${i}.dat  >> Av_BS.dat
#    cat avlig${i}.dat  >> Av_LIG.dat
#    cat avbmode${i}.dat >> Av_BM.dat
#    cat sdbs${i}.dat >> Av_BS_st.dat
#    cat sdlig${i}.dat >> Av_LIG_st.dat
#    cat sdbmode${i}.dat >> Av_BM_st.dat
#    cat popol${i}.dat >> Population.dat
#done
##
#rm top_energies_clusters.dat Av_EN.dat ST_EN.dat
#for i in $(seq 0 $end)
#do
#    sort -g -k2 cluster${i}_energies.dat | awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++}} END {print sum/i}' > aven${i}.dat 
#    sort -g -k2 cluster${i}_energies.dat | awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' > sden${i}.dat
#    cat aven${i}.dat  >> Av_EN.dat
#    cat sden${i}.dat  >> ST_EN.dat
#    sort -g -k2 cluster${i}_energies.dat |awk 'NR==1{print $2}' >> top_energies_clusters.dat
#done
##
#paste Population.dat Av_BS.dat Av_BS_st.dat Av_LIG.dat Av_LIG_st.dat Av_BM.dat Av_BM_st.dat Av_EN.dat ST_EN.dat top_energies_clusters.dat > ok.dat
#sed 's/,/./g' ok.dat > ok1.dat
#sed 's/[:[:blank:]:]/;/g' ok1.dat > summary-commas.dat
#echo "#Pop      Av_BS Av_BS_st Av_LIG Av_LIG_st Av_BM Av_BM_st Av_EN ST_EN top_energies_clusters" > summary.dat
#sed 's/;/ /g' summary-commas.dat >> summary.dat
#awk 'BEGIN{print"#Pop    <BS>    BSsd    <LIG>   LIGsd   <BM>    BMsd    <EN>    ENsd    ENtop"}{if (NR>1) {printf "%3d\t", $1; for (i=2; i<=NF; i++) {printf"%.2f\t", $i}; print ""}}' summary.dat > summary-tabs.dat
#mv summary.dat summary-ali-bm.dat
#mv summary-tabs.dat summary-tabs-ali-bm.dat
#mv summary-commas.dat summary-commas-ali-bm.dat
#######
##mv files to clusters directory
#mv av*.dat rmsdbs*.dat rmsdlig*.dat rmsdbmode*.dat clusters/ali-bm
##
##NOW ALIGNING THE BS ONLY, AND CALCULATIN THE RMSD OF BS, LIG, AND BM
##
##5.1. FOR EACH CLUSTER SEPARATELY
#for i in $(seq 0 $end)
#do
#    sed -e "s/CLUSTFILE/cl-traj.c${i}/"\
#	-e "s/XXX/${l}/"\
#	-e "s/YYY/${p}/"\
#	-e "s/ZZZ/${bs}/"\
#	-e "s/HOLOX/\.\.\/holo\_complex\_${p}\.pdbqt/"\
#	-e "s/rmsdlig/rmsdlig${i}/" \
#	-e "s/rmsdbs/rmsdbs${i}/" \
#	-e "s/rmsdbmode/rmsdbmode${i}/" \
#	-e "s/POPOL/popol${i}/" ${scridir}/fit-bs-and-deviation-ori.tcl \
#	| awk '{if ($1 !~ /^#/) {print $0}}'\
#	> fit-and-deviation${i}.tcl
#    vmd -dispdev text -e fit-and-deviation${i}.tcl
#done
#
##5.2. FOR THE TRAJ MADE BY CLUSTERS
#sed -e "s/XXX/${l}/"\
#    -e "s/YYY/${p}/"\
#    -e "s/ZZZ/${bs}/" ${scridir}/fit-bs-and-deviation-traj-ori.tcl \
#    | awk '{if ($1 !~ /^#/) {print $0}}'\
#    > fit-bs-and-deviation-traj.tcl
#vmd -dispdev text -e list-ordered-energy.tcl <<EOF
#source fit-bs-and-deviation-traj.tcl
#EOF
##mv aligned pdbs in folder
#mv frame*ali*pdb frames-ali
##mv clusters in folder
#mv cl.c*pdb cl-traj* clusters
#
##5. now perform statistics and generate the table
########
##topcut
#for i in $(seq 0 $end)
#do
#    awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++}} END {print sum/i}' rmsdbs${i}.dat > avbs${i}-top.dat
#    awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++}} END {print sum/i}' rmsdlig${i}.dat > avlig${i}-top.dat
#    awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++}} END {print sum/i}' rmsdbmode${i}.dat > avbmode${i}-top.dat
#    awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' rmsdbs${i}.dat > sdbs${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' rmsdlig${i}.dat > sdlig${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' rmsdbmode${i}.dat > sdbmode${i}.dat
#done
##
#rm Av_BS-top.dat Av_LIG-top.dat Av_BM-top.dat Av_BS_st.dat Av_LIG_st.dat Av_BM_st.dat Population.dat
#for i in $(seq 0 $end)
#do
#    cat avbs${i}-top.dat >> Av_BS-top.dat
#    cat avlig${i}-top.dat >> Av_LIG-top.dat
#    cat avbmode${i}-top.dat >> Av_BM-top.dat
#    cat sdbs${i}.dat >> Av_BS_st.dat
#    cat sdlig${i}.dat >> Av_LIG_st.dat
#    cat sdbmode${i}.dat >> Av_BM_st.dat
#    cat popol${i}.dat >> Population.dat
#done
##
#rm top_energies_clusters.dat Av_EN-top.dat ST_EN.dat
#for i in $(seq 0 $end)
#do
#    sort -g -k2 cluster${i}_energies.dat | awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++}} END {print sum/i}' > aven${i}-top.dat 
#    sort -g -k2 cluster${i}_energies.dat | awk 'BEGIN{i=0}{if (NR <='${topcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' > sden${i}.dat
#    cat aven${i}-top.dat >> Av_EN-top.dat
#    cat sden${i}.dat >> ST_EN.dat
#    sort -g -k2 cluster${i}_energies.dat |awk 'NR==1{print $2}' >> top_energies_clusters.dat
#done
##
#paste Population.dat Av_BS-top.dat Av_BS_st.dat Av_LIG-top.dat Av_LIG_st.dat Av_BM-top.dat Av_BM_st.dat Av_EN-top.dat ST_EN.dat top_energies_clusters.dat > ok.dat
#sed 's/,/./g' ok.dat > ok1.dat
#sed 's/[:[:blank:]:]/;/g' ok1.dat > summary-commas.dat
#echo "#Pop       Av_BS Av_BS_st Av_LIG Av_LIG_st Av_BM Av_BM_st Av_EN ST_EN top_energies_clusters" > summary.dat
#sed 's/;/ /g' summary-commas.dat >> summary.dat
#awk 'BEGIN{print"#Pop     <BS>    BSsd    <LIG>   LIGsd   <BM>    BMsd    <EN>    ENsd    ENtop"}{if (NR>1) {printf "%3d\t", $1; for (i=2; i<=NF; i++) {printf"%.2f\t", $i}; print ""}}' summary.dat > summary-tabs.dat
#mv summary.dat summary-ali-bs-top_clusters.dat
#mv summary-tabs.dat summary-tabs-ali-bs-top_clusters.dat
#mv summary-commas.dat summary-commas-ali-bs-top_clusters.dat
#######
##avcut
#for i in $(seq 0 $end)
#do
#    awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++}} END {print sum/i}' rmsdbs${i}.dat > avbs${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++}} END {print sum/i}' rmsdlig${i}.dat > avlig${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++}} END {print sum/i}' rmsdbmode${i}.dat > avbmode${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' rmsdbs${i}.dat > sdbs${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' rmsdlig${i}.dat > sdlig${i}.dat
#    awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' rmsdbmode${i}.dat > sdbmode${i}.dat
#done
##
#rm Av_BS.dat Av_LIG.dat Av_BM.dat Av_BS_st.dat Av_LIG_st.dat Av_BM_st.dat Population.dat
#for i in $(seq 0 $end)
#do
#    cat avbs${i}.dat  >> Av_BS.dat
#    cat avlig${i}.dat  >> Av_LIG.dat
#    cat avbmode${i}.dat >> Av_BM.dat
#    cat sdbs${i}.dat >> Av_BS_st.dat
#    cat sdlig${i}.dat >> Av_LIG_st.dat
#    cat sdbmode${i}.dat >> Av_BM_st.dat
#    cat popol${i}.dat >> Population.dat
#done
##
#rm top_energies_clusters.dat Av_EN.dat ST_EN.dat
#for i in $(seq 0 $end)
#do
#    sort -g -k2 cluster${i}_energies.dat | awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++}} END {print sum/i}' > aven${i}.dat 
#    sort -g -k2 cluster${i}_energies.dat | awk 'BEGIN{i=0}{if (NR <='${avcut}') {sum+=$2;i++;array[NR]=$2}} END {for(x=1;x<=i;x++){sumsq+=((array[x]-(sum/i))**2);}print sqrt(sumsq/i)}' > sden${i}.dat
#    cat aven${i}.dat  >> Av_EN.dat
#    cat sden${i}.dat  >> ST_EN.dat
#    sort -g -k2 cluster${i}_energies.dat |awk 'NR==1{print $2}' >> top_energies_clusters.dat
#done
##
#paste Population.dat Av_BS.dat Av_BS_st.dat Av_LIG.dat Av_LIG_st.dat Av_BM.dat Av_BM_st.dat Av_EN.dat ST_EN.dat top_energies_clusters.dat > ok.dat
#sed 's/,/./g' ok.dat > ok1.dat
#sed 's/[:[:blank:]:]/;/g' ok1.dat > summary-commas.dat
#echo "#Pop       Av_BS Av_BS_st Av_LIG Av_LIG_st Av_BM Av_BM_st Av_EN ST_EN top_energies_clusters" > summary.dat
#sed 's/;/ /g' summary-commas.dat >> summary.dat
#awk 'BEGIN{print"#Pop     <BS>    BSsd    <LIG>   LIGsd   <BM>    BMsd    <EN>    ENsd    ENtop"}{if (NR>1) {printf "%3d\t", $1; for (i=2; i<=NF; i++) {printf"%.2f\t", $i}; print ""}}' summary.dat > summary-tabs.dat
#mv summary.dat summary-ali-bs.dat
#mv summary-tabs.dat summary-tabs-ali-bs.dat
#mv summary-commas.dat summary-commas-ali-bs.dat
#######
##mv files to clusters directory
#mv av*.dat rmsdbs*.dat rmsdlig*.dat rmsdbmode*.dat clusters/ali-bs
##
#mv cluster*_energies.dat popol*.dat clusters
#
##6. now clean
##for i in {${start}..${end}}
#for i in $(seq ${start} ${end})
#do
#    rm cluster${i}_index_ready.dat\
#       cluster${i}_index.dat\
#       cluster${i}.dat\
#       fit-and-deviation${i}.tcl\
#       avbs${i}.dat\
#       avlig${i}.dat\
#       avbmode${i}.dat\
#       sdbs${i}.dat\
#       sdlig${i}.dat\
#       sdbmode${i}.dat\
#       aven${i}.dat\
#       sden${i}.dat\
#       fit-and-deviation${i}.tcl
#done
#rm Av_BS.dat Av_LIG.dat Av_BM.dat Av_BS_st.dat Av_LIG_st.dat Av_BM_st.dat Population.dat top_energies_clusters.dat Av_EN.dat ST_EN.dat cluster_listed cpptraj.in ok.dat ok1.dat topenergy_numerow.dat
##now perform final analysis to get file statistics.dat
#${scridir}/get-numbers-docking.sh $cutstat $numcom
