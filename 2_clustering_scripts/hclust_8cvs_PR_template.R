#############################################
#CV Clustering Analysis Script: Hierarchical Clustering 
##############################################
#### install.packages("randomForest", repos="http://cran.cnr.berkeley.edu")
###

# Install R packages
#install.packages("fpc")
#install.packages("cluster")
#install.packages("rgl")

# Recall packages
#setwd("./")
working_dir <- getwd()

library(fpc)
library(cluster)
#library(flashClust)
#library(rgl)

# Read the data and put it as a string
CV <- read.csv("FILENAME")
str(CV)
# Scale the data to allow the clustering | compulsory for hclust()
data <- scale(CV)

#### Choose the number of clusters 
k <- NCLUSTER

# Extract filename from working directory and "FILENAME"
filename <- basename(paste(working_dir, "/", "FILENAME", sep=""))

# Check if k is valid (between 1 and the maximum possible number of clusters)
max_clusters <- nrow(CV) 

if (k > max_clusters) {
  k <- max_clusters
  # to print only file name
  warning("Desired number of clusters (", NCLUSTER, ") exceeds the number of data points. Using maximum possible clusters (", max_clusters, ") for the file: ", filename)
}

cat("Using maximum possible clusters (", max_clusters, ") in window of size: ", max_clusters, "\n")

####
# Perform the clustering and cut the tree into k clusters
# 
d <- dist(data[,1:8],method = "euclidean")
#d <- dist(data[,5:5],method = "euclidean")

# d <- dist(data[,1:2],method = "euclidean") # you can select which columns to use from here, in case you do not they will all be used automatically.
# Standardize the formula for distance calculation to the actual collective variables used. For the rest, the size of the representative matrix can also be kept at the maximum size corresponding to the number of CVs in the file. It is not necessary to reduce it if a subset of CVs is used. If it is not reduced, all the CVs related to the clusters will be printed in the final file, but what determines the cluster analysis are the columns used in data and the formula to calculate the representatives. If errors like "index out of bounds" or errors in vector/matrix size occur, they are almost always attributable to errors of this kind. To perform the analysis with a single column, it can be selected in data[,1:1].

#hclust.fit1 <- flashClust(d, method="complete")  ## Change the method from here
#hclust.fit1 <- flashClust(d, method="ward")  ## Change the method from here
hclust.fit1 <- hclust(d, method="ward.D")  # Use hclust from stats package
groups1 <- cutree(hclust.fit1, k) # cut tree into k clusters

###### Extract the cluster representative for each cluster
index_min <- c()
#index <- matrix(nrow=nrow(cluster1), ncol=3)
matrice_dei_rappresentativi <- matrix(nrow=k, ncol=9)  # N° CVs +1 (row number) +2: window number + index
indice_finale <- c()

for(i in 1:k) {                                            # For each of the k clusters
cluster1 <- CV[groups1 ==i,]                               # Select this cluster
distanceMatrix <- matrix(nrow=nrow(cluster1), ncol=1)      # Build a matrix with a single column and the number of rows equal to the number of items in the cluster
#distanceMatrix <- matrix(nrow=length(cluster1), ncol=1)    # If we are working with a single collective variable

for(j in 1:nrow(cluster1)){                                # For each row (item) j of a cluster write in the line j of the distance matrix the distance between this point and the center of the cluster. The n-th coordinate of the center is calculated as the mean among all the n-th coordinates of the entries

#for(j in 1:length(cluster1)){

#distanceMatrix[j] = sqrt((cluster1[j]-mean(cluster1))^2)

#distanceMatrix[j] = sqrt((cluster1[j,1]-mean(cluster1[[1]]))^2 + (cluster1[j,2]-mean(cluster1[[2]]))^2 + (cluster1[j,3]-mean(cluster1[[3]]))^2)
distanceMatrix[j] = sqrt((cluster1[j,1]-mean(cluster1[[1]]))^2 + (cluster1[j,2]-mean(cluster1[[2]]))^2 + (cluster1[j,3]-mean(cluster1[[3]]))^2 + (cluster1[j,4]-mean(cluster1[[4]]))^2 + (cluster1[j,5]-mean(cluster1[[5]]))^2 + (cluster1[j,6]-mean(cluster1[[6]]))^2 + (cluster1[j,7]-mean(cluster1[[7]]))^2 + (cluster1[j,8]-mean(cluster1[[8]]))^2)

#distanceMatrix[j] = sqrt((cluster1[j,1]-mean(cluster1[[1]]))^2+(cluster1[j,2]-mean(cluster1[[2]]))^2)

}

index_min[i] <- which.min(distanceMatrix)   # It writes into the index_min vector the row number of each cluster corresponding to the object with the minimum "spatial" distance (considering all coordinates). However, this index is NOT referring to the initial item rows (in the CV matrix), but it refers to the row within each cluster.

matrice_dei_rappresentativi[i,1] <- row.names(cluster1[index_min[i],]) # I insert the row number from the original data of the chosen item.
matrice_dei_rappresentativi[i,2] <- cluster1[index_min[i],1] # and the eight collective variables.
matrice_dei_rappresentativi[i,3] <- cluster1[index_min[i],2]
matrice_dei_rappresentativi[i,4] <- cluster1[index_min[i],3]
matrice_dei_rappresentativi[i,5] <- cluster1[index_min[i],4]
matrice_dei_rappresentativi[i,6] <- cluster1[index_min[i],5]
matrice_dei_rappresentativi[i,7] <- cluster1[index_min[i],6]
matrice_dei_rappresentativi[i,8] <- cluster1[index_min[i],7]
matrice_dei_rappresentativi[i,9] <- cluster1[index_min[i],8]

indice_finale[i] <- row.names(cluster1[index_min[i],])  # Inside this vector, instead, I only put the row indices of the clusters.

}
write.table(indice_finale, file ="./OUTPUTNAME",row.names=FALSE)
write.table(matrice_dei_rappresentativi, file ="./OUTPUTNAMEmatrix",row.names=FALSE)

##########################################################
