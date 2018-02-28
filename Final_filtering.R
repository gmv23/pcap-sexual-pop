####################################### LOAD AND CLEAN DATA ##############################
snps <- read.table("capPF_cc_parental_reps_dip.012.pos")
colnames(snps) <- c("chrom", "pos")

indvs <- read.table("capPF_cc_parental_reps_dip.012.indv", stringsAsFactors = F)
indvs <-  unlist(indvs[,1])

geno <- read.table("capPF_cc_parental_reps_dip.012")
geno <- geno[,-1]
geno <- as.matrix(geno)
geno[geno==-1] <- NA

depths <- read.table("cc_diploid_depth_ratios.txt", header = T)

scaff <- read.table("scaffold_lengths.txt")

###################################### FILTER OUT HIGHLY HETEROZYGOUS ####################

get_het_rate <- function(x){
  hets <- sum(x==1, na.rm=T)
  total <- sum(!is.na(x))
  return(hets/total)
}
het_rates <- apply(geno, 2, get_het_rate)

pdf("Het_rates.pdf")
hist(het_rates,
     main = "Histogram of site heterozygosity",
     xlab = "Site heterozygosity")
dev.off()

high_hets <- which(het_rates>0.9)

print(paste(length(high_hets),"sites are removed due to heterozygosity > 0.9"))
snps <- snps[-high_hets,]
geno <- geno[,-high_hets]
depths <- depths[-high_hets,]

####################################### FILTER OUT WEIRD DEPTH RATIOS ##############################

depth_ratio_avg <- apply(depths[,3:ncol(depths)], 1, mean, na.rm=T)

pdf("Avg_depth_ratios.pdf")
hist(depth_ratio_avg,
     main = "Histogram of avg depth ratios",
     xlab = "Site avg depth ratios")
dev.off()

weird_depths <- which(depth_ratio_avg < 0.2 | depth_ratio_avg > 0.8)

print(paste(length(weird_depths),"sites are avg read depth ratio > 0.8 or < 0.2"))
snps <- snps[-weird_depths,]
geno <- geno[,-weird_depths]
depths <- depths[-weird_depths,]

###################################### FILTER OUT SMALL SCAFFOLDS #########################
colnames(scaff) <- c("scaffold", "length")
#Figure out how many markers are on each scaffold
markers_per_scaff <- as.data.frame(table(depths$CHROM))
colnames(markers_per_scaff) <- c("chrom", "num_markers")
markers_per_scaff$chrom <- as.integer(as.character(markers_per_scaff$chrom))
scaff$num_markers <- NA
scaff$num_markers[match(markers_per_scaff$chrom, scaff$scaffold)] <- markers_per_scaff$num_markers

pdf("scaffold_markers.pdf")
plot(scaff$length, scaff$num_markers,
     xlab = "Scaffold length",
     ylab = "Markers per scaffold",
     main = "Scaffold length vs markers per scaffold")
abline(v=300000)
dev.off()

print("Correlation between scaffold length and number markers per scaffold:")
cor(scaff$length, scaff$num_markers, use = "complete.obs")

large_scaffs <- scaff$scaffold[scaff$length > 300000]
large_markers <- which(depths$CHROM %in% large_scaffs)
n_markers_large <- length(large_markers)

print(paste(nrow(depths)-n_markers_large,"markers on scaffolds < 300 kb are removed"))
snps <- snps[large_markers,]
geno <- geno[,large_markers]
depths <- depths[large_markers,]

######################################## PRUNE #######################################

#Prune 1 marker for each 1 kb window

sites_to_remove <- c()
for(i in unique(depths$CHROM)){
  window <- c(1,1000)
  end <- F
  remove_from_scaffold <- 0
  while(end == F){
    window_snps <- which((snps$chrom == i) & (snps$pos >= window[1]) & (snps$pos <= window[2]))
    if(length(window_snps) != 0){
      to_remove <- sample(window_snps, length(window_snps)-1,F)
      sites_to_remove <- c(sites_to_remove, to_remove)
      remove_from_scaffold <- remove_from_scaffold + length(to_remove)
    }
    if(window[2] >= max(snps$pos[snps$chrom==i])){
      end <- T
    }
    window <- window + 1000
  }
  print(paste(remove_from_scaffold, "of", sum(snps$chrom==i), "markers removed from scaffold", i))
}

print(paste("Total of", length(sites_to_remove), "markers removed by pruning"))
print(paste("Final data set of", nrow(snps)-length(sites_to_remove), "markers"))

snps_prune <- snps[-sites_to_remove,]
geno_prune <- geno[,-sites_to_remove]
depths_prune <- depths[-sites_to_remove,]
########################################## MAKE PARENTAL CONSENSUS ####################

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
cons_0664 <- apply(geno[grep("664", indvs),], 2, get_consensus)
cons_06180 <- apply(geno[grep("6180", indvs),], 2, get_consensus)

cons_0664_prune <- apply(geno_prune[grep("664", indvs),], 2, get_consensus)
cons_06180_prune <- apply(geno_prune[grep("6180", indvs),], 2, get_consensus) 

#Remove all 664 and 6180 isolates from files
parental_pos <- grep("664|6180", indvs)
parental_names <- grep("664|6180", indvs, value=T)

geno <- geno[-parental_pos,]
geno_prune <- geno_prune[-parental_pos,]
indvs <- indvs[-parental_pos]

# Make object with consensus parental genotypes
final_geno <- rbind(cons_0664, cons_06180, geno)
total_indvs <- c("cons_0664", "cons_06180", indvs)

final_geno_prune <- rbind(cons_0664_prune, cons_06180_prune, geno_prune)
 
########################################## SAVE NEW FILES #################################

# To make same 012 format that vcftools spits out, need to turn NAs back into -1s
# and need to make column of row numbers
prep_to_print <- function(x){
  x[is.na(x)] <- -1
  x <- cbind(0:(nrow(x)-1), x)
  return(x)
}

#Final 012 files
write.table(prep_to_print(final_geno_prune), "capPF_final.012", sep="\t", col.names = F, row.names=F)
write.table(total_indvs, "capPF_final.012.indv", quote=F, col.names=F, row.names=F)
write.table(snps_prune, "capPF_final.012.pos", quote=F, col.names=F, row.names=F)

#Unpruned
write.table(prep_to_print(final_geno), "capPF_final_unpruned.012", sep="\t", col.names = F, row.names=F)
write.table(snps, "capPF_final_unpruned.012.pos", quote=F, col.names=F, row.names=F)


