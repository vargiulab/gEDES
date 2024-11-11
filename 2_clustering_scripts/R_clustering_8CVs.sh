#!/bin/bash
#@@@@@@@@@@ Script to perform a multistep cluster analysis with hierarchical clustering and k-means. @@@@@@
#@@@@@@@@@@ made by Andrea Basciu and refined by Han Kurt, Mohd Athar, and Attilio V. Vargiu @@@@@@@@@@@@@@

# The first clustering gives a predefined number of structures serving as starting point for k-means.
# Moreover, multistep clustering is performed so that a minimum number of structures is found for each
# the Nw "windows" in which we divide the first CVs (generally the RoG of the BS).

# $$$$$$$$$$$$$$$$$$$$ USAGE $$$$$$$$$$$$$$$$ #
# First, run the script without arguments, it will provide optimal value for the width of each window for the first CV
# Next, use this value to run the clustering:
# parameters N, Center, Sigma must be provided by the user when launching the script
#

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REMEMBER TO SET THE CORRECT BS IN trajanalyze_templ.tcl !!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "REMEMBER TO SET THE CORRECT BS IN trajanalyze_templ.tcl"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
sleep 5s

# grab RoG range from index_RoGBS.dat
RoG_min=$(sort -gk2 index_RoGBS.dat | head -1 | awk '{print $2}')
RoG_max=$(sort -gk2 index_RoGBS.dat | tail -1 | awk '{print $2}')
RoG_difference=$(echo "$RoG_max - $RoG_min" | bc)

# to calculate the number of windows, we se the maximum between between 10 and the ratio between the CV range (max-min) and a \"sensible\" variation of the same CV, here set to 0.2.
sigmasens=0.2
Nw=$(echo "$RoG_difference/$sigmasens"| bc)
if [ $Nw -lt 10 ]
then
    Nw=10
fi
    
if [ $# -le 1 ]; then
echo "";
echo "Warning: The script requires three input parameters"
echo "[N] = Number of desired clusters (for example, N=500)"
#echo "[Center] = Center of the distribution in terms of RMSD from the apo";
echo "[Sigma] = Width of each window"
echo "[minclwin] = Minimum number of clusters in each window (should be such that the number of required clusters is equal or larger than number of windows times minclwin"
echo "Provide in input: [N] [Sigma] [minclwin]"
echo "The number of windows is suggested by default to be the maximum between 10 and the ratio between the CV range (max-min) and a \"sensible\" variation of the same CV, here set to 0.2. In this case the value is: $Nw (set parameter $Nw). Depending on their width (sigma), they";
echo "will contain more or fewer structures. They will be built starting from the minimum RoG upwards.";
echo "";
echo "From the analysis of the file related to RoG, it appears that RoG_max - RoG_min = $RoG_difference";
echo "To cover $RoG_difference in $Nw windows, it is recommended to use the following value for the width of the windows:";
suggested_sigma=$(echo "$RoG_difference/$Nw"| bc -l)
printf "%.3f" $suggested_sigma
echo "";
echo "";
exit 1;
fi

N=$1       # No. of clusters you want. The number of total clusters may differ slightly from this number.
#center=$2   # Center of the distribution (assumed approximately Gaussian)
sigma=$2    # Dev. distribution standards
#step=0.1     # Width (in terms of sigma) of the (half) windows. Note: the total number of windows remains fixed
minclwin=$3    # Minimum number of clusters in each window

#set the maximum number of clusters for each window to N / (num windows) [ rounded...]
maxclwin=$(echo "(${N}/${Nw}+0.5)/1" | bc)

# Note: with sigma = 0.4 and step 0.2, each window is 0.4*sigma wide
#
##################################################
# Define the ends of the windows: upper (top) and lower (bottom)
##################################################
##### Walls for Window 1 #########
top[1]=$(echo "$RoG_min+$sigma"|bc -l)
btm[1]=$(echo "$RoG_min"|bc -l)
echo "Le finestre sono quindi: Window 1 = [${btm[1]}; ${top[1]}]"
for i in `seq 2 $Nw`; do
    iminus1=$(echo "$i -1"|bc)
#    topvalue=$(echo "$btm[$i]+$sigma"|bc -l)
#    echo "top[1] = ${top[$iminus1]}"
#    echo "topvalue = $topvalue"
    declare btm[$i]=${top[$((${i}-1))]}
    declare top[$i]=$(echo "${btm[${i}]}+${sigma}" | bc -l)
    echo "Window $i: [${btm[${i}]};${top[${i}]}]"
done

###################################################
# Let's get the structures with the RMSDs inside the chosen windows
# It could be done with a for loop on awk, but I haven't succeeded, at least trying dynamic variables
###################################################
for i in `seq 1 $Nw`; do
awk -v top1=${top[${i}]} -v btm1=${btm[${i}]} '{if ($2 < top1 && $2 >= btm1) print $0}' index_RoGBS.dat > RoGBS-sigma_win${i}.dat
echo "## Divide in slices ##"
echo "Slice $i between ${btm[${i}]} and ${top[${i}]}"
echo "###"
done

#####################################################################################
########## Extract clusters within each window  ###################
echo ""
echo ""
printf "%s\n" "Number of clusters for each slice"
N_tot=$(wc -l index_RoGBS.dat)
N_tot_lines=$(echo $N_tot | cut -d ' ' -f1)
for i in `seq 1 $Nw`; do
    a=$(wc -l RoGBS-sigma_win$i.dat)
    num_lines=$(echo $a | cut -d ' ' -f1)
    #echo "num_lines: $num_lines"
    N_cluster=$(echo "$N*($num_lines/$N_tot_lines)"|bc -l)
    #echo "N_cluster: $N_cluster"
    N_cl_rounded=$(echo "($N_cluster+0.9)/1" | bc)
    #minclwin=$(echo "($minclwin+0.1)/1" | bc)
    echo "$N_cl_rounded ${minclwin}"
    if [ "$N_cl_rounded" -lt "${minclwin}" ]
       #if [ "$N_cl_rounded" -lt 2 ];
    then
	#cluster_num[$i]=$(echo "${minclwin}")
	cluster_num[$i]=${minclwin}
	echo "==> [Slice$i]: cluster number artificially increased to ${minclwin}!"
    else
	if [ "$N_cl_rounded" -gt "${maxclwin}" ]
	then
	    cluster_num[$i]=${maxclwin}
	else
	    #echo "N_cl_rounded: $N_cl_rounded"
	    cluster_num[$i]=${N_cl_rounded}
	    #echo "${cluster_num[$i]}"
	fi
    fi
    printf "%s" "Slice" $i ": "
    printf "%.0f"  ${cluster_num[$i]}
    printf "%s" " | Structures in the slice: " $num_lines
    printf "\n"
done
###
# Note: Due to approximation errors, it may be that the total number of clusters exceeds N.
# To solve this issue, if needed, the number of clusters in the last window could simply be calculated as N minus the clusters from the other windows.
###

echo -n "N_total_clusters=\$(echo " > tmp_source
ee=$((${Nw}-1))
for i in `seq 1 $ee`; do
    echo -n "\${cluster_num[$i]}+" >> tmp_source
done
echo "\"\${cluster_num[$Nw]}\" | bc -l)" >> tmp_source

source tmp_source
rm tmp_source
#N_total_clusters=$(echo "${cluster_num[1]}+${cluster_num[2]}+${cluster_num[3]}+${cluster_num[4]}+${cluster_num[5]}+${cluster_num[6]}+${cluster_num[7]}+${cluster_num[8]}+${cluster_num[9]}+${cluster_num[10]}+${cluster_num[11]}+${cluster_num[12]}+${cluster_num[13]}+${cluster_num[14]}+${cluster_num[15]}" | bc -l)

echo "Total Number of cluster: $N_total_clusters"

############## Parte 2 #####################
# Go from RMSD from apo, to files for cluster analysis, which includes the RoGBS-meta of all residues, in the right format
############################################
# The first part of the loop goes from the RMSD file from the apo to printing a file with index, various RMSDs, window, rmsd from the holo, from the apo, etc.
# The last command in awk instead puts only the second block of data into a file, containing the parameters for the cluster analysis
# The last part could be avoided by having files formatted differently
echo " "
echo "...I save the CVs for each window to do the clustering on.."
echo " "
for i in `seq 1 $Nw`; do
    echo "Working on the window: " $i "..."
    if [ -f lines_sigma${i}.dat ]
    then
	rm lines_sigma${i}.dat
    fi
    awk '{printf"\\b%s\\b\n",$1}' RoGBS-sigma_win$i.dat > sigma_${i}_ready.dat
    while IFS='' read -r line || [[ -n "$line" ]]
    do
	grep -m 1 "${line} " CVS_clustering.dat >> lines_sigma${i}.dat
    done < sigma_${i}_ready.dat
    awk '{print $2}' lines_sigma${i}.dat > dataclust_sigma${i}.dat
done

############## Part 3 #####################
# Now we do the cluster analysis with R: hierarchical agglomerative ward.D2 method.
# We start from an R template, adapt it to the newly created files and launch the cluster analysis
############################################
echo "Preparing files for R..."
for i in `seq 1 $Nw`; do
    sed "s/FILENAME/dataclust_sigma${i}.dat/" hclust_8cvs_PR_template.R > Rhclust_win${i}_tmp.R
    sed "s/NCLUSTER/${cluster_num[i]}/" Rhclust_win${i}_tmp.R > Rhclust_win${i}_tmp2.R
    sed "s/OUTPUTNAME/clusterselected_s${i}.dat/" Rhclust_win${i}_tmp2.R > Rhclust_win${i}.R
    rm Rhclust_win${i}_tmp.R Rhclust_win${i}_tmp2.R
    R < Rhclust_win${i}.R --no-save
done
##############

####!!!!for i in `seq 1 11`; do
####!!!!    grep -v "x" clusterselected_s$i.dat > tmp1
####!!!!    sed 's/"//g' tmp1 > tmp2
####!!!!    awk '{printf"\\b%s\\b\n",$1}' tmp2 > tmp3
####!!!!    while IFS='' read -r line || [[ -n "$line" ]]
####!!!!    do
####!!!!	grep -m 1 "${line} " holoRMSD.dat >> cluster_s${i}_RMSDfromholo.dat
####!!!!	awk '{printf "%s", $1}' tmp3 > index_cluster_s${i}
####!!!!	
####!!!!    done < tmp2
####!!!!    rm tmp1 tmp2 tmp3
####!!!!   done

##################### SAVING THE PDBS FROM HIERARCHICAL CLUSTERING ###########################

#mkdir frames_hierarchical
#for i in `seq 1 15`; do
 #  frame=$(grep -v "x" clusterselected_s${i}.dat | sed 's/\"//g' | awk '{printf "%s ", $0}')
 # sed "s/FRAMESTOSAVE/$frame/" lines_template.sh > lines${i}.sh
# chmod +x lines${i}.sh
# ./lines${i}.sh
#done

############################################################################
####################### NOW THE K-MEANS PART ###############################
############################################################################

awk '{print $2}' CVS_clustering.dat > CVS_clustering_for_kmeans.dat
for i in `seq 1 $Nw`;
do grep -v "V1" clusterselected_s${i}.datmatrix > tmp;
   sed 's/"/ /g' tmp > tmp${i};
done
cat tmp? tmp?? > tmp_total

awk '{print $2" "$3" "$4" "$5" "$6" "$7" "$8" "$9}' tmp_total > tmp_total_ok
rm tmp_total
awk '!visited[$0]++' tmp_total_ok > tmp_clean
rm tmp_total_ok

awk '{printf "%s,",$1}' tmp_clean > CV0_center
awk '{printf "%s,",$2}' tmp_clean > CV1_center
awk '{printf "%s,",$3}' tmp_clean > CV2_center
awk '{printf "%s,",$4}' tmp_clean > CV3_center 
awk '{printf "%s,",$5}' tmp_clean > CV4_center
awk '{printf "%s,",$6}' tmp_clean > CV5_center
awk '{printf "%s,",$7}' tmp_clean > CV6_center
awk '{printf "%s,",$8}' tmp_clean > CV7_center 

Nclustkmeans=$(wc -l tmp_clean)
Ntotlineskmeans=$(echo $Nclustkmeans | cut -d ' ' -f1)

cat CV?_center | sed 's/\,$//' > CV_ALL_center
cat k-means-head.R CV_ALL_center k-means_8CVs-tail.R > k-means_tmp.R
sed "s/CENTERS/$Ntotlineskmeans/g" k-means_tmp.R > k-means_ready.R
R < k-means_ready.R --no-save

############################################################################
##################### SAVING THE PDBS FROM HIERARCHICAL-KMEANS CLUSTERING ###########################

mkdir frames_hierarchkmeans
frame=$(grep -v "x" matrix_k-means-init-centers.dat | sed 's/\"//g' | awk '{printf "%s ", $0}')
sed "s/FRAMESTOSAVE/$frame/" lines_templatekmeans.sh > lines_kmeans.sh
chmod +x lines_kmeans.sh
./lines_kmeans.sh

########################## CONCATENATE THE FRAMES ################################

cd frames_hierarchkmeans
if [ -f frames.pdb ]
then
    echo "frames.pdb already exists, moving to frames.pdb.bak. CHECK..."
    mv frames.pdb frames.pdb.bak
fi
cat *.pdb > frames.pdb
cd ..

################# ANALYSIS OF RMSD WITH RESPECT TO THE VARIOUS HOLO MODELS  #######################

cd frames_hierarchkmeans
mkdir clusters_aligned
#cp ../trajanalyze_templ.tcl .
ln -sf ../rmsdsymm.cpptrj
./rmsdsymm.cpptrj
mv aligned* clusters_aligned
cd ..
mv frames_hierarchkmeans frames_hierarchkmeans_${N}cl_${sigma}sigma_${minclwin}mincl_ward

#rm tmp trajanalyze_templ.tcl
#done
