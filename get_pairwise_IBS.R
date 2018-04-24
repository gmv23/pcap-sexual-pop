#!/usr/bin/env Rscript

# This script takes a VCF file and creates a pairwise IBS similarity matrix 
# for use in identifying and filtering out clones

getwd()

snps <- read.table("./capPF_filtered.012.pos")
head(snps)
geno <- read.table("./capPF_filtered.012")
print("geno loaded")
geno <- geno[,-1]
geno <- t(geno)
indvs <- read.table("./capPF_filtered.012.indv", stringsAsFactors = F)
indvs <-unlist(indvs$V1)
head(indvs)
geno[geno==-1] <- NA
print('finished replacing NAs')

ibs <- function(x,y){
  
  if(length(x) != length(y)){
    stop("Not good")}
  
  total <- 0
  number_sites <- length(x)
  
  for(i in 1:length(x)){
    if(is.na(x[i]) | is.na(y[i])){
      number_sites <- number_sites - 1
      next
    }else if(x[i]==y[i]){
      add <- 2
    }else if(abs(x[i] - y[i]) == 1){
      add <- 1
    }else if(abs(x[i] - y[i]) == 2){
      add <- 0
    }
    
    total <- total + add
  }
  return(total/(2*number_sites))
}


d <- ncol(geno)

IBS_matrix <- matrix(nrow=d, ncol=d)

print("got to loop")

for(i in 1:(d-1)){
	for (j in (i +1):d){
		IBS_matrix[i,j] <- ibs(geno[,i], geno[,j])
	}
	print(i)
}

rownames(IBS_matrix) <- indvs
colnames(IBS_matrix) <- indvs

write.csv(IBS_matrix, "IBS_matrix.csv")
