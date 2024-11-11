module load plumed
plumed driver --mf_dcd meta_concat.dcd --plumed plumed-common.dat
awk 'BEGIN{i=0}(NR>1){i++;print i,$2}' COLVAR_apoMD > index_RoGBS.dat
awk 'BEGIN{i=0}{if($1 !~ /#/) {i++; printf"%s %s,%s,%s,%s,%s,%s\n", i,$2,$3,$4,$5,$6,$7}}' COLVAR_apoMD > CVS_clustering.dat
awk '{print $2}' CVS_clustering.dat > CVS_clustering_for_kmeans.dat
