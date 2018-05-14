###########################################################################################

# This script analyzes population structure in the PF and lab F1 population
# Principle components analysis
# Fst
# Separation of generations by inbreeding coefficient and Mendelian error proportion
# Distribution of MAF and marker heterozygosity by year

###############################################   LOAD AND CLEAN DATA     ###########################################

setwd("~/Documents/work/Smart_lab/P_capsici/Sexual_pop/analysis/structure_all/")

library(RColorBrewer)
library(pcaMethods)
library(StAMPP)
library(vioplot)

#Load SNPs data frame
snps <- read.table("capPF_final.012.pos")
colnames(snps) <- c("chrom", "pos")

#Load individuals data frame
indvs <- read.table("capPF_final.012.indv", stringsAsFactors = F)
indvs <-  unlist(indvs[,1])

#Load and clean geno data frame
geno <- read.table("capPF_final.012")
geno <- geno[,-1]
geno <- as.matrix(geno)
geno[geno==-1] <- NA

###############################################   POPULATE ISOLATES DATA FRAME    ###########################################

#Get name, year, plate, and lane info from individual names
isolates <- cbind("isolate" = sapply(indvs, function(x) return(unlist(strsplit(x,":"))[1])), 
                  "year" = substr(sapply(indvs, function(x) return(unlist(strsplit(x,":"))[1])),1,2),
                  "plate" = sapply(indvs, function(x) return(unlist(strsplit(x,":"))[2])),
                  "lane" = sapply(indvs, function(x) return(unlist(strsplit(x,":"))[3])))

rownames(isolates) <- NULL
isolates <- as.data.frame(isolates, stringsAsFactors=F)

isolates$year[isolates$year=="co"] <- "P"
isolates$year[isolates$year=="68"] <- "Lab"
isolates$year <- as.factor(isolates$year)

year_colors <- brewer.pal(9, "Set1")
isolates$color <- year_colors[match(isolates$year, sort(unique(isolates$year)))]

###############################################    PCA     ###########################################

### PCA including all isolate and parents

# PCA (scaled and centered matrix)
PF_pca <- pca(geno, nPcs=4, method='nipals', 
              scale='uv', center=TRUE)

#Percent variance explained
PF_pve <- round(PF_pca@R2*100,2)

#Make plots using first four PCs
old.par <- par(no.readonly=T)
par(mar=c(5,5,4,2), xpd=NA)

for(i in 1:(PF_pca@nPcs-1)){
  for(j in (i+1):PF_pca@nPcs){
    
    pdf(paste("plots/PF_pca",i,j,".pdf",sep=""))

    plot(PF_pca@scores[,i],
         PF_pca@scores[,j],
         col = isolates$color,
         pch = 16,
         xlab = paste("PC ", i, ": ", PF_pve[i], "%", sep=""),
         ylab = paste("PC ", j, ": ", PF_pve[j], "%", sep=""),
         cex.lab=1.5,
         cex.axis=1.2)
    
    text(c(PF_pca@scores[1,i], PF_pca@scores[2,i]),
         c(PF_pca@scores[1,j], PF_pca@scores[2,j]),
         labels = c("0664", "06180"),
         cex=1.1, 
         pos=4, xpd = NA)
    
    legend("topright",
           legend = sort(unique(isolates$year)),
           fill = year_colors,
           cex=0.9,
           ncol=2)
    
    dev.off()
    
  }
}

par(old.par)

#Which are the isolates that cluster with 0664?
isolates$isolate[PF_pca@scores[,1] > 75]

## Now excluding parents and lab isolates

geno.field <- geno[-c(1:2, which(isolates$year == "Lab")),]
isolates.field <- isolates[-c(1:2, which(isolates$year == "Lab")),]

# PCA (scaled and centered matrix)
PF_pca <- pca(geno.field, nPcs=4, method='nipals', 
              scale='uv', center=TRUE)

#Percent variance explained
PF_pve <- round(PF_pca@R2*100,2)

#Make plots
old.par <- par(no.readonly=T)
par(mar=c(5,5,4,2), xpd=NA)

for(i in 1:(PF_pca@nPcs-1)){
  for(j in (i+1):PF_pca@nPcs){
    
    pdf(paste("plots/PF_pca_field",i,j,".pdf",sep=""))
    
    plot(PF_pca@scores[,i],
         PF_pca@scores[,j],
         col = isolates.field$color,
         pch = 16,
         xlab = paste("PC ", i, ": ", PF_pve[i], "%", sep=""),
         ylab = paste("PC ", j, ": ", PF_pve[j], "%", sep=""),
         cex.lab=1.5,
         cex.axis=1.2)
    
    legend("topright",
           legend = sort(unique(isolates$year))[-9],
           fill = year_colors[-9],
           cex=0.9,
           ncol=2)
    
    dev.off()
    
  }
}

par(old.par)

##################################################### PAIRWISE FST ###########################################

#Turn geno into allele frequency format
geno.stamp <- geno/2

#Make column names SNP names
colnames(geno.stamp) <- paste("S_",snps$chrom, "_", snps$pos, sep="")

#Add isolate, 'population' (year), ploidy, and format information

geno.stamp <- cbind("Isolate" = isolates$isolate,
                    "Year" = as.character(isolates$year),
                    "Ploidy" = 2,
                    "Format" = 'freq',
                    geno.stamp)
geno.stamp <- as.data.frame(geno.stamp, stringsAsFactors = F)

#Make object for StAMPP functions
geno.stamp <- stamppConvert(geno.stamp, "r")

#Weir and Cockerham's Fst, 100 bootstraps, 95% CI
fst.stamp <- stamppFst(geno.stamp,100,95)

#Get rid of parents
fst.mat <- fst.stamp$Fsts
fst.mat <- fst.mat[-1,-1]

#Put pairwise fst matrix in year order
year.order <- c("Lab", "09", "10", "11", "12", "13", "14", "15")

#Hacky function to reorder pairwise Fst matrix in year order
order_matrix <- function(mat, ord){
  
  mat.sort <- matrix(NA, nrow=nrow(mat), ncol=ncol(mat))
  rownames(mat.sort) <- ord
  colnames(mat.sort) <- ord
  
  for(i in 1:(nrow(mat.sort)-1)){
    ord.i <- ord[i]
    for(j in (i+1):ncol(mat.sort)){
      ord.j <- ord[j]
      mat.values <- c(mat[ord.i,ord.j], mat[ord.j,ord.i])
      mat.value <- mat.values[which(!is.na(mat.values))]
      mat.sort[i,j] <- mat.value
    }
  }
  return(mat.sort)
}

#Order pairwise Fst matrix and p-value matrix
fst.sort <- order_matrix(fst.mat, year.order)
p.sort <- order_matrix(fst.stamp$Pvalues[-1,-1], year.order)

#Get red color gradient for plotting
red_colors <- colorRampPalette(brewer.pal(9, 'Oranges'), space = "Lab")(30)[-(1:2)]

# 'Rotate' matrix to get right position for plotting
rotate <- function(x) t(apply(x,2,rev))
fst.rot <- rotate(rotate(rotate(fst.sort)))

#Make pairwise Fst plot

pdf("plots/Pairwise_fst.pdf")

image(fst.rot,
      bty = "n", xaxt = 'n', yaxt='n',
      col = red_colors,
      main = "Pairwise Fst across years")

axis.pos <- seq(0,1,length.out=nrow(fst.sort))

axis(1, axis.pos, labels = rownames(fst.rot))
axis(2, axis.pos, labels = colnames(fst.rot))

fst.plot <- rotate(rotate(fst.sort))
p.plot <- rotate(rotate(p.sort))

for(i in 2:nrow(fst.plot)){
  for(j in 1:(i-1)){
    if(p.plot[i,j] < 0.05){
      font.type = 2 #Make text bold if Fst is significant
      cex.type = 0.9
    }else{
      font.type = 1
      cex.type = 0.75
    }
    print(p.plot[i,j])
    text(axis.pos[j], rev(axis.pos)[i], 
         round(fst.plot[i,j], 3),
         cex = cex.type,
         font = font.type)
  }
}

dev.off()

# Plot pairwise Fst as a function of time

pdf("plots/fst_by_time_red.pdf")

year_diffs <- c()
fsts <- c()
fst.cols <- c()

fst.year <- fst.plot[rownames(fst.plot) != "Lab", colnames(fst.plot) != "Lab"]

for(i in 2:nrow(fst.year)){
  for(j in 1:(i-1)){
    year.i <- as.integer(rownames(fst.year))[i]
    year.j <- as.integer(colnames(fst.year))[j]
    fst <- fst.year[i,j]
    year_diff <- year.j - year.i
    year_diffs <- c(year_diffs, year_diff)
    fsts <- c(fsts,fst)
    if(year.i == 15 | year.j == 15){
      fst.col <- 'red'
    }else if(year.i == 14 | year.j == 14){
      fst.col <- 'gray'
    }else{
      fst.col <- 'gray'
    }
    fst.cols <- c(fst.cols, fst.col)
  }
}

plot(year_diffs, fsts, pch=16, col = fst.cols,
     main = "Pairwise Fst vs pairwise time difference",
     ylab = "Pairwise Fst",
     xlab = "Time difference (years)")

dev.off()

###############################################    INBREEDING OVER TIME     ###########################################

#Calculate individual heterozygosity

get_heterozygosity <- function(x){
  hets <- sum(x==1, na.rm=T)
  total <- sum(!is.na(x))
  return(hets/total)
}

het_levels <- apply(geno,1,get_heterozygosity)
boxplot(het_levels ~ isolates$year)
stripchart(het_levels ~ isolates$year, vertical=T, add=T, pch=16)

#Get expected mafs based on parental genotypes

get_maf <- function(x){
  if(any(is.na(x))){
    return(NA)
  }else{
    return(min(sum(x)/4, 1-sum(x)/4))
  }
}

expected_mafs <- apply(geno[1:2,],2,get_maf)

#Calculate inbreeding coefficient (using parental genotypes to get He)

Ffun <- function(geno, mafs){
  if (length(geno) != length(mafs)){
    return('ERROR!')
  } else {
    maf_na <- which(is.na(mafs))
    geno <- geno[-maf_na]
    mafs <- mafs[-maf_na]
    miss = which(is.na(geno))
    if (length(miss) == 0){
      Ho = sum(geno == 1)/length(geno)
      He = mean(2*mafs*(1-mafs))
    } else {
      Ho = sum(geno[-miss] == 1,na.rm=TRUE)/length(na.omit(geno))
      He = mean(2*mafs[-miss]*(1-mafs[-miss]))
    }
    return(1-Ho/He)
  }
}

Fin <- apply(geno,1,Ffun, mafs=expected_mafs)

### Make the plot

groups <- c("Lab", "09", "10", "11", "12", "13", "14", "15")

pdf("plots/Inbreeding.pdf")

plot("n",
     xlim = c(0.75,length(groups)+0.25),
     ylim = c(min(Fin),max(Fin)),
     xlab = "Year",
     ylab = "Inbreeding coefficient",
     main = "Inbreeding coefficients over time",
     xaxt='n')

rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
     col = adjustcolor("gray",alpha.f=0.8))

for(i in seq(-1,1,by=0.2)){
  abline(h=i,
         col = "white",
         lty=2)
}

for(i in 1:length(groups)){
  vioplot(Fin[which(isolates$year==groups[i])],
          add=T,
          at = i,
          col=year_colors[i],
          names=groups[i],
          drawRect=F)
  stripchart(Fin[which(isolates$year==groups[i])],
          add=T,
          at=i,
          vertical=T,
          pch=20,
          method = "jitter",
          jitter = 0.025)
}

box("plot", col="white", lty=1, lwd=5)
axis(1, at=1:length(groups), labels=c("Lab", "2009", "2010", "2011", "2012", "2013", "2014", "2015"))

dev.off()

#Which are the more heterozygous lab isolates?
cbind(isolates$isolate[isolates$year == "Lab" & Fin > 0],
      Fin[isolates$year == "Lab" & Fin > 0])

############################################ MAF AND HET DISTRIBUTION   ###########################################

find_maf <- function(x){
  total <- 2*sum(!is.na(x))
  if(total == 0){
    maf <- NA
  }else{
    a1_f <- sum(x, na.rm = T)
    if(a1_f/total <= 0.5){
      maf <- a1_f/total
    }else{
      maf <- 1 - a1_f/total
    }  
  }
  return(maf)
}

# Plot MAF distributions

pdf("plots/maf_distributions.pdf")

plot(0, type="n", xlim = c(0,0.55), ylim=c(0,5.5),
     main = "MAF distributions by year",
     xlab = "MAF",
     ylab = "Density")
counter <- 0
for(group in groups){ 
  counter <- counter + 1
  marker_mafs <- apply(geno[isolates$year==group,], 2, find_maf)
  lines(density(marker_mafs, na.rm=T, adjust=1.50),
        lwd = 2,
        main = group,
        col = year_colors[counter])
}
legend(0.05,5.7,
       legend = groups,
       fill = year_colors,
       bty = "n",
       ncol = 4,
       cex= 1,
       x.intersp = 0.3, 
       text.width = 0.1,
       y.intersp = 1.25)

dev.off()

# Plot heterozygosity distributions

pdf("plots/het_distributions.pdf")

plot(0, type="n", xlim = c(0,1), ylim=c(0,4),
     main = "Marker heterozygosity distributions by year",
     xlab = "Marker heterozygosity",
     ylab = "Density")
counter <- 0
for(group in groups){ 
  counter <- counter + 1
  marker_hets <- apply(geno[isolates$year==group,], 2, get_heterozygosity)
  lines(density(marker_hets, na.rm=T, adjust=1.75),
        lwd = 2,
        main = group,
        col = year_colors[counter])
}
legend(0.05,4,
       legend = groups,
       fill = year_colors,
       bty = "n",
       ncol = 4,
       cex= 1,
       x.intersp = 0.3, 
       text.width = 0.2,
       y.intersp = 1.25)

dev.off()

#####

#PCA of just 2015 isolates

# PCA (scaled and centered matrix)
PF_pca <- pca(geno[isolates$year=="15",], nPcs=4, method='nipals', 
              scale='uv', center=TRUE)

#Percent variance explained
PF_pve <- round(PF_pca@R2*100,2)

Fin_15 <- Fin[which(isolates$year=="15")]
colors_15 <- rep('gray', length(Fin_15))
colors_15[Fin_15 < -0.2] <- "red"

plot(PF_pca@scores[,1], PF_pca@scores[,2], col=colors_15, pch=16)

###############################################    MENDELIAN ERRORS   ###########################################

#Make matrix where 1=ME and 0=non ME
#Rows are individuals excluding parents
#Columns are markers

n <- nrow(geno)
p <- ncol(geno)
me_matrix <- matrix(NA,nrow=n, ncol=p)

for(i in 1:n){
  for(j in 1:p){
    g <- geno[i,j]
    p1 <- geno[1,j]
    p2 <- geno[2,j]
    if(is.na(p1) | is.na(p2) | is.na(g) |
       (p1==1 & p2==1)){
      me_matrix[i,j] <- NA
    }else if(
    (p1==0 & p2==0 & g!=0) |
    (p1==2 & p2==2 & g!=2) |
    (p1==2 & p2==0 & g!=1) |
    (p1==0 & p2==2 & g!=1) |
    (p1==1 & p2==0 & g==2) |
    (p1==0 & p2==1 & g==2) |
    (p1==1 & p2==2 & g==0) |
    (p1==2 & p2==1 & g==0)){
      me_matrix[i,j] <- 1
    }else{
      me_matrix[i,j] <- 0
    }   
  }
}

ind_mes <- apply(me_matrix,1,function(x)sum(x,na.rm=T)/sum(!is.na(x)))
site_mes <- apply(me_matrix[3:n,],2,function(x)sum(x,na.rm=T)/sum(!is.na(x)))

hist(ind_mes)
hist(site_mes)
boxplot(ind_mes ~ isolates$year)

# Find which SNPS have high tendency of MEs in putative F1s and remove them

# Which are (conservatively) putative F1s --- just based on lab F1, 2009, and 2010
f1_coords <- which(isolates$year %in% c("Lab", "09", "10") & Fin < 0)

# Get rid of any markers with more than 10% MEs in putative F1s
f1_site_mes <- apply(me_matrix[f1_coords,],2,function(x)sum(x,na.rm=T)/sum(!is.na(x)))
hist(f1_site_mes, breaks=50)
me_enriched_snps <- which(f1_site_mes > 0.10)
length(me_enriched_snps)
me_snps <- snps[me_enriched_snps,]
# Write file with info on ME enriched SNPs for later analysis
me_snps.csv <- write.csv(me_snps, "me_enriched_snps.csv", quote = F, row.names = F)

# Get corrected ME statistics
ind_mes_corrected <- apply(me_matrix[,-me_enriched_snps],1,function(x) sum(x,na.rm=T)/sum(!is.na(x)))
site_mes_corrected <- apply(me_matrix[,-me_enriched_snps],2,function(x) sum(x,na.rm=T)/sum(!is.na(x)))

#Make ME plot
pdf("plots/Mendelian_errors.pdf")

plot("n",
     xlim = c(0.75,length(groups)+0.25),
     ylim = c(0,max(ind_mes_corrected, na.rm=T)),
     xlab = "Year",
     ylab = "Proportion Mendelian errors",
     main = "Proportion of Mendelian errors per isolate",
     xaxt='n')

rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
     col = adjustcolor("gray",alpha.f=0.8))

for(i in seq(0,1,by=0.025)){
  abline(h=i,
         col = "white",
         lty=2)
}

for(i in 1:length(groups)){
  vioplot(ind_mes_corrected[which(isolates$year==groups[i])],
          add=T,
          at = i,
          col=year_colors[i],
          names=groups[i],
          drawRect=F)
  stripchart(ind_mes_corrected[which(isolates$year==groups[i])],
             add=T,
             at=i,
             vertical=T,
             pch=20,
             method="jitter",
             jitter = 0.025
              )
}

box("plot", col="white", lty=1, lwd=5)
axis(1, at=1:length(groups), c("Lab", "2009", "2010", "2011", "2012", "2013", "2014", "2015"))

dev.off()

#Which are the lab cross individuals with many MEs?
cbind(isolates$isolate[ind_mes_corrected > 0.1 & isolates$year == "Lab"],
      ind_mes_corrected[ind_mes_corrected > 0.1 & isolates$year == "Lab"])

### Assign isolates to generation based on MEs

# Get F1 cutoff based on lab cross F1s 
lab_me_proportions <- ind_mes_corrected[isolates$year=="Lab" & ind_mes_corrected < 0.1]
lab_me_proportions.sd <- sqrt(var(lab_me_proportions))
f1_cutoff <- mean(lab_me_proportions.sd + 4*lab_me_proportions.sd)

isolates$gen[ind_mes_corrected < f1_cutoff] <- "F1"
isolates$gen[ind_mes_corrected >= f1_cutoff] <- "Inbred"
isolates$gen[isolates$year=="P"] <- "Parent"

#Write isolates table
write.csv(isolates, "isolate_info.csv", quote = F, row.names = F)

####################### 'F1s' color coded by year ###################

F1_geno <- geno[isolates$gen %in% c("F1", "Parent"),]
F1_isolates <- isolates[isolates$gen %in% c("F1", "Parent"),]

# PCA (scaled and centered matrix)
F1_pca <- pca(F1_geno, nPcs=4, method='nipals', 
              scale='uv', center=TRUE)

#Percent variance explained
F1_pve <- round(F1_pca@R2*100,2)

old.par <- par(no.readonly=T)
par(mar=c(5,5,4,2), xpd=NA)

for(i in 1:(F1_pca@nPcs-1)){
  for(j in (i+1):F1_pca@nPcs){
    
    pdf(paste("plots/PF_pca_F1",i,j,".pdf",sep=""))
    
    plot(F1_pca@scores[,i],
         F1_pca@scores[,j],
         col = F1_isolates$color,
         pch = 16,
         xlab = paste("PC ", i, ": ", F1_pve[i], "%", sep=""),
         ylab = paste("PC ", j, ": ", F1_pve[j], "%", sep=""),
         main = "PCA of putative F1s",
         cex.lab=1.5,
         cex.axis=1.2)
    
    legend("topright",
           legend = sort(unique(isolates$year)),
           fill = year_colors,
           cex=0.9,
           ncol=2)
    dev.off()
    
  }
}

par(old.par)




