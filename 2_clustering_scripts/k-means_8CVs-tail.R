), nrow=CENTERS, ncol=8) 

index <- c()
 model <- kmeans(CV[,1:8],start,iter.max=10000)
 
 #calculate indices of closest instance to centroid
 for (i in 1:k){
   rowsum <- rowSums(abs(CV[which(model$cluster==i),1:8] - model$centers[i,]))
     index[i] <- as.numeric(names(which.min(rowsum)))
     }
 index
#model$cluster
a <- model$centers

write.table(index, file ="./matrix_k-means-init-centers.dat",row.names=FALSE)
#write.table(index, file ="./matrix_k-means-init-random.dat",row.names=FALSE)
