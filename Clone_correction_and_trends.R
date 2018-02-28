#!/usr/bin/env Rscript

# This script takes a pairwise IBS matrix, identifies clonal groups, and
# filters genotype files in order to just keep one individual from each clonal group
# Current version is modified from R script written on local machine 1/15/17

getwd()
library(igraph)
print("Done with igraph")

################################# READ IBS MATRIX AND GENOTYPE FILES ##############################

#Load IBS matrix and make individuals row and column names
IBS_matrix <- read.csv("IBS_matrix.csv", header = T)
row.names(IBS_matrix) <- IBS_matrix$X
IBS_matrix$X <- NULL
IBS_matrix <- as.matrix(IBS_matrix)
print("Loaded IBS")

#Load marker and individual info
snps <- read.table("capPF_filtered.012.pos")
indvs <- read.table("capPF_filtered.012.indv", stringsAsFactors = F)
indvs <-unlist(indvs$V1)
print("indvs")

#Load genos and change -1s to NAs
geno <- read.table("capPF_filtered.012")
print("Loaded geno")
geno <- geno[,-1]
print("column done")
geno <- as.matrix(geno)
geno[geno==-1] <- NA
print("NAs done")

############################### LOOK AT IBS BETWEEN PARENTAL REPLICATES ############################

#Pull out entries of parental replicates
IBS_664 <- IBS_matrix[grep("664", indvs), grep("664", indvs)]
IBS_6180 <- IBS_matrix[grep("6180", indvs), grep("6180", indvs)]

IBS_664_values <- as.vector(IBS_664[!is.na(IBS_664)])
IBS_6180_values <- as.vector(IBS_6180[!is.na(IBS_6180)])

#Plot histogram of pairwise IBS between parental replicates
pdf("plots/IBS_between_parents.pdf")
old.par <- par(no.readonly=T)
par(mar=c(5,5,4,2))
hist(c(IBS_6180_values, IBS_664_values), xlim = c(0.95, .99),
     main = "IBS between parental replicates",
     xlab = "IBS",
     cex.main = 2, cex.lab = 2, cex.axis =2,
     col = "lightblue")
par(old.par)
dev.off()

############################### ASSIGN INDIVIDUALS TO CLONAL GROUPS ################################

#Turn high IBS cells of matrix to 1
modify_matrix <- function(x){
  if(is.na(x) | x<.95){
    return(0)
  }else{
    return(1)
  }
}
clone_or_not <- structure(sapply(IBS_matrix, modify_matrix), dim=dim(IBS_matrix))

# Create network -> Each isolate is a node and there is an edge between isolates that are clones
g <- graph_from_adjacency_matrix(clone_or_not, "undirected")

# Clusters are isolates that only have edges between themselves and not the rest of the network (ie clones)
g.clusters <- clusters(graph = g)

#### Create table of clonal group assignments

# Make list of cluster size corresponding to each member of network (used later)
cluster_sizes <- rep(NA, length(indvs))
for(i in 1:length(cluster_sizes)){
  member <- g.clusters$membership[i]
  size <- sum(g.clusters$membership == member)
  cluster_sizes[i] <- size
}

# Prepare table and variables for loop
clonal_groups <- 1:(g.clusters$no)
clone_assignments <- matrix(ncol=2)
colnames(clone_assignments) <- c("Sample", "Clonal_group")
counter <- 0

# Assign individuals to clonal groups starting with largest group
for(i in 1:length(unique(g.clusters$csize))){ #loop through all unique cluster sizes
  # Start with largest cluster size
  current_size <- sort(unique(g.clusters$csize), decreasing=T)[i] 
  # how many groups of this size are there
  same_size_clonal_groups <- unique(g.clusters$membership[cluster_sizes == current_size]) 
  #loop through groups of that size
  for(j in 1:length(same_size_clonal_groups)){ 
    counter <- counter +1
    old_clonal_group_id <- same_size_clonal_groups[j] #Assignment to group from g.clusters$membership
    new_clonal_group_assignment <- clonal_groups[counter] #New assignment going from largest to smallest
    clone_assignments <- rbind(clone_assignments, cbind(
      indvs[which(g.clusters$membership == old_clonal_group_id)],
      new_clonal_group_assignment))
  }
}
clone_assignments <- clone_assignments[-1,]
clone_assignments <- as.data.frame(clone_assignments, stringsAsFactors = F)
clone_assignments$Clonal_group <- as.integer(clone_assignments$Clonal_group)

## Get rid of 0664 reintroductions from geno, indvs, IBS_matrix, and clone_assignments

group_664 <- unique(clone_assignments$Clonal_group[grep("664", clone_assignments$Sample)])
reintros <- grep("664",clone_assignments$Sample[clone_assignments$Clonal_group==group_664], invert=T, value=T)
reintros_pos <- which(indvs %in% reintros)

geno_short <- geno[-reintros_pos,]
indvs_short <- indvs[-reintros_pos]
IBS_matrix_short <- IBS_matrix[-reintros_pos, -reintros_pos]
clone_assignments_short <- clone_assignments[-which(clone_assignments$Sample %in% reintros),]

################################ CREATE CONSENSUS PARENTAL GENOTYPES ################################

#### Get consensus 0664 and 06180 genotypes

# Function works on a vector containing genotypes 
# (single marker across multiple individuals in -1,01,2 format)

get_consensus <- function(x){
  l <- length(x)
  # NA if over half are NAs
  if(sum(is.na(x)) >= floor(l/2)){
    return(NA)
  # Majority otherwise. If tie, NA.
  }else{
    counts <- as.data.frame(table(x))
    if(max(counts$Freq) <= ceiling(sum(counts$Freq)/2)){
      return(NA)
    }else{
      return(as.integer(as.character(counts$x[which.max(counts$Freq)])))
    }
  }
}

#Get parental consensus
cons_0664 <- apply(geno_short[grep("664", indvs_short),], 2, get_consensus)
cons_06180 <- apply(geno_short[grep("6180", indvs_short),], 2, get_consensus)

#Remove all 664 and 6180 isolates from files
parental_pos <- grep("664|6180", indvs_short)
parental_names <- grep("664|6180", indvs_short, value=T)

geno_prog <- geno_short[-parental_pos,]
indvs_prog <- indvs_short[-parental_pos]
IBS_matrix_prog <- IBS_matrix_short[-parental_pos, -parental_pos]
clone_assignments_prog <- clone_assignments_short[-which(clone_assignments_short$Sample %in% parental_names),]

# Make object with consensus parental genotypes
parental_geno <- rbind(cons_0664, cons_06180)

############################## SAMPLE ONE ISOLATE FROM EACH CLONAL GROUP ##############################

#Rename clone_assignments_prog to cap to make it easier to handle
cap <- clone_assignments_prog

#Sample all but one from each clonal group to remove
clones_to_remove <- c()
for(i in unique(cap$Clonal_group)){
  clones <- cap$Sample[cap$Clonal_group==i]
  clones_toss <- sample(clones, size = length(clones)-1, replace=F)
  clones_to_remove <- c(clones_to_remove, clones_toss)
}

#Generate clone corrected (cc) data sets
clone_assignments_cc <- cap[-which(cap$Sample %in% clones_to_remove),]
clones_to_remove_pos <- which(indvs_prog %in% clones_to_remove)
indvs_cc <- indvs_prog[-clones_to_remove_pos]
geno_cc <- geno_prog[-clones_to_remove_pos,]
IBS_matrix_cc <- IBS_matrix_prog[-clones_to_remove_pos, -clones_to_remove_pos]

#######################################   WRITE NEW FILES   ##########################################

setwd("./filtered_data/")

# To make same 012 format that vcftools spits out, need to turn NAs back into -1s
# and need to make column of row numbers
prep_to_print <- function(x){
  x[is.na(x)] <- -1
  x <- cbind(0:(nrow(x)-1), x)
  return(x)
}

# Clonal group assignments -- both before removing reintros and parents and after
write.csv(clone_assignments, "Clonal_groups_including_0664_reintros.csv", row.names = F)
write.csv(clone_assignments_prog, "Clonal_groups_progeny_no_reintros.csv", row.names = F)

#Progeny geno and indvs files
write.table(prep_to_print(geno_prog), "capPF_progeny.012", sep="\t", col.names = F, row.names=F)
write.table(indvs_prog, "capPF_progeny.012.indv", quote=F, col.names=F, row.names=F)
write.csv(IBS_matrix_prog, "IBS_matrix_progeny.csv")

#Consensus parental genotypes
write.table(prep_to_print(parental_geno), "capPF_parental_consensus.012", sep="\t", col.names=F, row.names=F)

#Clone corrected files
write.table(prep_to_print(geno_cc), "capPF_cc.012", sep="\t", col.names = F, row.names=F)
write.table(indvs_cc, "capPF_cc.012.indv", quote=F, col.names=F, row.names=F)
write.csv(IBS_matrix_cc, "IBS_matrix_cc.csv")

setwd("..")

######################################## LOOK AT TRENDS IN CLONALITY ####################################

setwd("./plots/")

### How many genotyped isolates passed filters per year and how many of them are genotypically unique

#Extract year from sample name
clone_assignments_prog$Year <- sapply(clone_assignments_prog$Sample, FUN = function(x) 
  return(unlist(strsplit(x, "PF"))[1]))

#Sum number of genotypes passing filters to this point by year
isolate_sums <- as.data.frame(table(clone_assignments_prog$Year))
colnames(isolate_sums) <- c("Year", "Isolates")
isolate_sums$Clonal_types <- NA

#How many unique clonal groups are there in each year
for(year in isolate_sums$Year){
  isolate_sums$Clonal_types[isolate_sums$Year==year] <- 
    length(unique(clone_assignments_prog$Clonal_group[clone_assignments_prog$Year==year]))
}

## Total isolates and number clonal types per year

pdf("isolate_sums.pdf")
colrs <- c("skyblue", "dodgerblue4") #Can change me if desired
isolate_sums.bp <- t(isolate_sums[,-1])

#Save x coordinates from barplot
isolate_sums.x <- barplot(isolate_sums.bp,
                          beside=T,
                          space=c(0,1.5),
                          plot=F)

#Make empty plotting window with axes and labels
plot(0, type="n",
     xlim = c(0, max(isolate_sums.x)*1.10),
     ylim = c(0, max(isolate_sums$Isolates*1.10)), 
     xaxt="n", xlab = "Year", ylab = "Number")

#Add horizontal grey lines
for(i in seq(0,max(isolate_sums$Isolates)*1.10, by=5)){
  abline(h=i, 
         col = adjustcolor("gray", alpha.f=0.3))
}

#Add bars
barplot(isolate_sums.bp, beside=T,
        names.arg = isolate_sums$Year,
        col = colrs,
        space = c(0,1.5),
        add=T)

#Add legend
legend("topright", legend = c("Total isolates", "Unique genotypes"),
      fill = colrs, bg="white")

dev.off()

## Distribution of clonal_group sizes by year 
## This code is bare bones and graph looks very bad

pdf("clonal_sizes_by_year.pdf", height = 10, width=7)

old.par <- par(no.readonly = T)
par(mfcol = c(length(isolate_sums$Year), 1))

#Get size distribution across all years and max group size
size_distribution_overall <- as.data.frame(table(table(clone_assignments_prog$Clonal_group)))
colnames(size_distribution_overall) <- c("Size", "Groups")
max_group_size <- max(as.integer(as.character(size_distribution_overall$Size)))

#Plot distribution year by year
for(year in isolate_sums$Year){
  # Get year size distribution
  cap_year <- clone_assignments_prog[clone_assignments_prog$Year==year,]
  year_distribution_short <- as.data.frame(table(table(cap_year$Clonal_group)))
  # Pad with 0s
  year_distribution <- as.data.frame(cbind("Size" = 1:max_group_size, "Number" = 0))
  year_distribution$Number[year_distribution$Size %in% year_distribution_short$Var1] <-
    year_distribution_short$Freq
  # Normalize so that variable is percentage of isolates represented by group of x size
  year_distribution$Percent <- year_distribution$Number*year_distribution$Size/sum(year_distribution$Size*year_distribution$Number)
  # Plot year distribution
  barplot(height=year_distribution$Percent, names=year_distribution$Size)
}
par(old.par)
dev.off()


## Clonal group size distribution

size_dist <- as.data.frame(cbind("Size" = 1:max_group_size, "Number" = 0))
size_dist$Number[size_dist$Size %in% size_distribution_overall$Size] <-
  size_distribution_overall$Groups

pdf("clonal_size_distribution.pdf", width=10)
old.par <- par(no.readonly = T)
par(mar=c(5,5,4,2))
barplot(height = size_dist$Number, names.arg=size_dist$Size,
        main = "Distribution of clonal group sizes",
        xlab = "Number of samples per unique genotype",
        ylab = "Frequency",
        cex.main =2, cex.lab =1.8, cex.names = 1.5,
        col = "lightblue")
par(old.par)
dev.off()

setwd("..")
