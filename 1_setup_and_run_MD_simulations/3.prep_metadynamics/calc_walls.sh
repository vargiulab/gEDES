#!/bin/bash

R_min1=$1
R_max1=$2
#R_min2=$3
#R_max2=$4

R_diff1=$(echo "$R_max1 - $R_min1" | bc -l)
R_step1=$(echo "$R_diff1/389" | bc -l)

#R_diff2=$(echo "$R_max2 - $R_min2" | bc -l)
#R_step2=$(echo "$R_diff2/800" | bc -l)

echo "R_diff1=$R_diff1"
#echo "R_diff2=$R_diff2"
echo "R_step1=$R_step1"
#echo "R_step2=$R_step2"

index1=5
index2=1

#while [ $index1 -lt 255000000 ]
while [ $index1 -lt 389 ] || [ $index2 -lt 389 ];
do

    for i in `seq 55250000 500000 250000000`; do
	
	r_1=$(echo "scale=3; ($R_max1 - $index2*($R_step1))/1"|bc )
#	r_2=$(echo "scale=3; ($R_max2 - $index2*($R_step2))/1"|bc )
	
	echo "STEP$index1=$i  AT$index1=$r_1"
	
	index1=$(echo "$index1 +1" | bc)
	index2=$(echo "$index2 +1" | bc)

#	index1= `expr $index1 + 1`
	#index2= `expr $index2 + 1` 

    done

done



#####R_min1=$1
#####R_max1=$2
#####R_min2=$3
#####R_max2=$4
#####
#####R_diff1=$(echo "$R_max1 - $R_min1" | bc -l)
#####R_step1=$(echo "$R_diff1/200000000" | bc -l)
#####
#####R_diff2=$(echo "$R_max2 - $R_min2" | bc -l)
#####R_step2=$(echo "$R_diff2/200000000" | bc -l)
#####
#####echo "R_diff1=$R_diff1"
#####echo "R_diff2=$R_diff2"
#####echo "R_step1=$R_step1"
#####echo "R_step2=$R_step2"
#####
#####index1=5
#####index2=1
#####
######while [ $index1 -lt 255000000 ]
#####while [ $index1 -lt 255000000] || [ $index2 -lt 255000000 ];
#####do
#####
#####    for i in `seq 55000000 250000 255000000`; do
#####	
#####	r_1=$(echo "$R_max1 - $index2*($R_step1)"|bc -l)
#####	r_2=$(echo "$R_max2 - $index2*($R_step2)"|bc -l)
#####	
#####	echo "STEP$index1=$i  AT$index1=$r_1,$r_2"
#####	
#####	index1=$(echo "$index1 +1" | bc)
#####	index2=$(echo "$index2 +1" | bc)
#####
######	index1= `expr $index1 + 1`
#####	#index2= `expr $index2 + 1` 
#####
#####    done
#####
#####done


