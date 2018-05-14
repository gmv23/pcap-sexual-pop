# This script assigns markers to LGs based on genetic map in Lamour et al. 2012
# Then produces histograms of allele depth ratios for each isolate and for each LG within each isolate
# Then attempts to estimate a ploidy level for each linkage group
# Finally prints out some plots summarizing ploidy trends

############################################   LOAD AND CLEAN DATA   ######################################
print("Beginning to import data")
depths <- read.table("depth_ratios.txt", header = T)

indvs <- read.table("capPF_filtered.012.indv", stringsAsFactors = F)
indvs <-  unlist(indvs[,1])

lg_order <- read.csv("scaffold_assignments.csv", header=T, stringsAsFactors = F)
lg_order <- lg_order[-nrow(lg_order),] #Get rid of totals row

subscaffolds <- read.csv("subscaffold_assignments.csv", header=T, stringsAsFactors = F)

#Function to divide string like 'Sc45.2' or 'LG01.01' into first or second number
get_value <- function(x, position){
  x.number.start <- regexpr("[0-9]", x) #Find where numbers start
  x.clean <- substr(x,x.number.start,nchar(x)) #Remove beginning letters
  x.split <- strsplit(x.clean, ".", fixed = T)[[1]] #Split on "."
  return(as.integer(x.split[position])) #Return first or second number
}

#Function to get marker (2nd part) from string like "1_529343"
get_marker <- function(x){
  start_pos <- regexpr("_", x)
  marker <- as.integer(substr(x, start_pos + 1, nchar(x)))
  return(marker)
}

#Separate Scaffold_Block into Scaffold and Block
lg_order$Scaffold <- sapply(lg_order$Scaffold_Block, get_value, position = 1)
lg_order$Block <- sapply(lg_order$Scaffold_Block, get_value, position = 2)

subscaffolds$Scaffold <- sapply(subscaffolds$Scaffold_Block, get_value, position = 1)
subscaffolds$Block <- sapply(subscaffolds$Scaffold_Block, get_value, position = 2)

#Clean up marker positions
subscaffolds$Lowest_Marker <- sapply(subscaffolds$Lowest_Marker, get_marker)
subscaffolds$Highest_Marker <- sapply(subscaffolds$Highest_Marker, get_marker)

print("Finished loading and cleaning data")

##################################   COMBINE GENETIC MAP DATASETS INTO ONE DATA FRAME   #######################

print("Pulling information from genetic map tables to assign scaffolds to linkage groups")

block_assignments <- matrix(NA, nrow=nrow(lg_order), ncol=6)
colnames(block_assignments) <- c("Scaffold", "Block", "Low", "High", "LG", "LG_order")

#Loop through lg_order by row
#Take out values or match values in subscaffolds data frame to populate block_assignments data frame
for(i in 1:nrow(block_assignments)){
  row_data <- unlist(lg_order[i,]) 
  #Pull out values from each row of lg_order
  scaffold <- as.integer(row_data[5])
  block <- as.integer(row_data[6])
  lg <- as.integer(get_value(row_data[1], 1))
  lg.pos <- as.integer(get_value(row_data[1], 2))
  #If scaffold isn't divided into subscaffold blocks, "low" and "high" markers cover entire scaffold
  if(is.na(block)){ 
    low <- 1
    high <- as.integer(row_data[3])
  #Otherwise find matching marker positions of subscaffold from subscaffolds data frame
  }else{
    low <- subscaffolds$Lowest_Marker[subscaffolds$Scaffold==scaffold & subscaffolds$Block==block]
    high <- subscaffolds$Highest_Marker[subscaffolds$Scaffold==scaffold & subscaffolds$Block==block]
  }
  block_assignments[i,] <- c(scaffold, block, low, high, lg, lg.pos) #Fill in row in block_assignments df
}

block_assignments <- as.data.frame(block_assignments)

#Look at block_assignments in order by scaffold and scaffold block
head(block_assignments[order(block_assignments$Scaffold, block_assignments$Block),])

print("Finished assigning scaffolds to linkage groups")

############################################   ASSIGN POSITIONS IN VCF FILE TO LGs    ##################################

print("Assign markers to linkage groups")

# Insert LG column to depths data frame after chromosome and position
depths <- data.frame(depths[1:2], "LG"=rep(NA, nrow(depths)), depths[3:ncol(depths)])

#Loop through rows in depths data frame
#Use CHROM and POS values to find matching LG from block_assignments data frame
for(i in 1:nrow(depths)){
  chrom <- depths[i,1]
  pos <- depths[i,2]
  #How many times does scaffold appear in block_assignments
  #0 means not assigned to lg, 1 means not divided into subscaffolds, 2+ means divided into subscaffold blocks
  scaffold_block_number <- sum(block_assignments$Scaffold==chrom) 
  if(scaffold_block_number==0){ #If scaffold not assigned a LG, return NA
    lg <- NA
  }else if(scaffold_block_number==1){ #If no subdivision of scaffold, match LG
    lg <- block_assignments$LG[block_assignments$Scaffold==chrom]
  }else if(scaffold_block_number > 0){ #If subdivision of scaffold, find appropriate scaffold block
    lg <- NA #NA if POS is not in any of the scaffold blocks
    chrom_blocks <- block_assignments[block_assignments$Scaffold==chrom,] #Pull out info on scaffold blocks
    for(j in 1:nrow(chrom_blocks)){ #Loop through scaffold blocks
      if(pos >= chrom_blocks$Low[j] & pos <= chrom_blocks$High[j]){
        lg <- chrom_blocks$LG[j] #If marker position is in between low and high, pull out matching LG
      }
    }
  }
  depths$LG[i] <- lg #Populate LG column of depths data frame
}
print("Finished assignign markers to linkage groups")

###############################################   DRAW HISTOGRAMS FOR ALL ISOLATES AND LGS    ######################################

print("Beginning to create histograms")
table(depths$LG) #Total number of markers per linkage group

#Function to pull out abbreviated isolate name
get_name <- function(x){
  name <- unlist(strsplit(x, "[:.]"))[1]
  if(substr(name,1,1) == "X"){ #Remove X if from column name
    return(substr(name,2,nchar(name)))
  }else{
    return(name)
  }
}

#For each isolate, make 2 plots:
#One including all depth ratios
#One separating histograms by linkage group

for(i in 4:ncol(depths)){ #Loop through individuals
  isolate_name <- get_name(colnames(depths)[i])
  
  #First by linkage group
  pdf(paste("./plots/all_isolates_by_lg/", isolate_name, "_by_lg.pdf", sep=""), height=11.5, width=8)
  
  #Save par settings before switching to 6x3 layout
  old.par <- par(no.readonly = T)
  par(mfrow=c(6,3))
  
  for(j in 1:18){ #Loop through linkage groups
    #Count number of markers
    n <- sum(depths$LG==j & !is.na(depths[,i]), na.rm = T)
    if(n > 0){ #Don't draw any histograms for lgs with 0 markers
      #QUESTION -- plot just major/total or minor/total as well to make it symmetrical?
      hist(c(depths[depths$LG==j,i], 1-depths[depths$LG==j,i]),
           main = paste("LG:", j, "n:", n),
           xlab = paste(isolate_name, "Depth ratio"),
           col = "gray")
      #Draw lines at 0.25, 0.33, 0.5, 0.66, 0.75
      abline(v = 0.333, lty=3)
      abline(v=0.666, lty=3)
      abline(v = 0.25, lty=2)
      abline(v=0.75, lty=2)
      abline(v=0.5, lty=1)
    }
  }
  
  par(old.par) #Go back to old par settings
  dev.off()
  
  #Now do all markers for isolate
  pdf(paste("./plots/all_isolates/", isolate_name, ".pdf", sep=""))
  
  n <- sum(is.na(depths[,i]), na.rm = T)
  
  #QUESTION -- plot just major/total or minor/total as well to make it symmetrical?
  hist(c(depths[,i], 1-depths[,i]),
       main = paste(isolate_name, "n:", n),
       xlab = "Depth ratio",
       col = "gray")
  abline(v = 0.333, lty=3)
  abline(v=0.666, lty=3)
  abline(v = 0.25, lty=2)
  abline(v=0.75, lty=2)
  abline(v=0.5, lty=1)
  
  dev.off()
  
}
print("Finished creating histograms")
###############################################   ASSIGN PLOIDY LEVEL TO LGs    ######################################
print("Estimating ploidy levels of each linkage group")
# Function to roughly estimate ploidy by finding histogram peak
# Somewhat hacky ----
# Take vector of depth ratios
# Make bins of 0.25, 0.33, 0.5, 0.66, and 0.75 +/- 0.10
# Find number of depth ratios that appear in each bin
# Find which bin has greatest count and return corresponding ploidy

estimate_ploidy <- function(x){
  centers <- c(0.25,0.33,0.5,0.66,0.75)
  ploidies <- c(4,3,2,3,4)
  counts <- rep(0, length(centers))
  for(i in x[!is.na(x)]){
    for(j in 1:length(centers)){
      center <- centers[j]
      if(i > (center - 0.1) & i < (center + 0.1)){
        counts[j] <- counts[j] + 1
      }  
    }
  }
  ploidy_est <- ploidies[which.max(counts)]
  return(ploidy_est)
}

#Make file of non-diploid individuals
isolate_ploidies <- apply(depths[,4:ncol(depths)], 2, estimate_ploidy)
non_diploid_isolates <- indvs[isolate_ploidies != 2]
write.table(non_diploid_isolates, "non_diploid_isolates.txt", quote = F, row.names = F, col.names = F)

#Make matrix with rows as isolates and columns as linkage groups
#Each cell corresponds to ploidy estimate for that isolate/LG combination

LG_ploidies <- matrix(NA, nrow=(ncol(depths)-3), ncol=18)
colnames(LG_ploidies) <- paste("LG", 1:18, sep = "")
rownames(LG_ploidies) <- 1:(ncol(depths)-3)

for(i in 4:ncol(depths)){ #Loop through individuals
 
  #Get isolate name and fill in row name
  row_pos <- i-3
  isolate_name <- get_name(colnames(depths)[i])
  rownames(LG_ploidies)[row_pos] <- isolate_name
  
  #Loop through linkage groups
  for(j in 1:18){
    lg_markers <- depths[depths$LG==j,i]
    lg_markers <- lg_markers[!is.na(lg_markers)]
    if(length(lg_markers) < 50){
      cell_value <- NA
    }else{
      cell_value <- estimate_ploidy(c(lg_markers, 1-lg_markers))
    }
    LG_ploidies[row_pos,j] <- cell_value
  }
}

LG_ploidies <- as.data.frame(LG_ploidies)
print("Finished estimating ploidy levels")
#############################################   INCORPORATE CLONE INFORMATION    ######################################
print("Reading in clonal information")
clones <- read.csv("Clonal_groups_including_0664_reintros.csv",
                   header=T, stringsAsFactors = F)

clones$Sample <- sapply(clones$Sample, get_name)
clonal_assignment <- clones$Clonal_group[match(rownames(LG_ploidies), clones$Sample)]
LG_ploidies <- data.frame("Name" = rownames(LG_ploidies), 
                          "Clonal_group"= clonal_assignment, 
                          LG_ploidies)
#Order by clonal group
LG_ploidies <- LG_ploidies[order(LG_ploidies$Clonal_group),]

write.csv(LG_ploidies, "ploidy_est_non_CC.csv")
# No ploidy aberrations consistent within clonal groups
print("Finished creating csv file with clonal group info and ploidy levels")

###############################################   LOOK AT TRENDS AND PLOT    ######################################

print("Plotting trends")
### Distribution of ploidy levels by linkage group

pdf("./plots/ploidy_by_lg_bp.pdf")
old.par <- par(no.readonly = T)
par(xpd=NA)

#Get matrix of counts for ploidy levels across isolates for each linkage group
lg_ploidy_counts <- matrix(NA, ncol=18, nrow=4)
for(i in 3:20){
  lg_ploidy_counts[,(i-2)] <- tabulate(LG_ploidies[,i], nbins=4)
}
lg_ploidy_counts <- lg_ploidy_counts[-1,]
colnames(lg_ploidy_counts) <- 1:18
rownames(lg_ploidy_counts) <- 1:3

#Make barplot
bp_colors <- c("burlywood", "plum4", "cadetblue3")
barplot(lg_ploidy_counts, col=bp_colors,
        main = "Ploidy levels by linkage group",
        ylab = "Counts",
        xlab = "Linkage group",
        ylim = c(0,max(apply(lg_ploidy_counts,2,sum))*1.15))
legend("topright",
       legend = c("Diploid", "Triploid", "Tetraploid"),
       fill = bp_colors,
       ncol = 3,
       x.intersp = 0.2,
       text.width = 2,
       bty="n",
       cex=0.85)

dev.off()

### Look at distribution of non-diploid linkage groups by isolate

pdf("./plots/ploidy_frequencies.pdf")

#Count non-diploid linkage groups for each isolate
non_diploid_counts <- apply(LG_ploidies[,3:20], 1, function(x) return(sum(x!=2, na.rm=T)))

#Make table of counts
non_diploid_table <- matrix(NA, nrow=19, ncol=2)
for(i in 0:18){
  i.count <- sum(non_diploid_counts==i, na.rm=T)
  non_diploid_table[(i+1),] <- c(i, i.count)
}

#Make data frame with plot information
non_diploid_table <- as.data.frame(non_diploid_table)
colnames(non_diploid_table) <- c("Number", "Counts")

#Color information is used to mask 0 values
non_diploid_table$Color <- "plum4"
non_diploid_table$Color[non_diploid_table$Counts==0] <- "white"
non_diploid_table$Border <- "black"
non_diploid_table$Border[non_diploid_table$Counts==0] <- "white"

barplot(height = (non_diploid_table$Counts), 
        names.arg = non_diploid_table$Number,
        col = non_diploid_table$Color,
        border = non_diploid_table$Border,
        xlab = "Number of non-diploid linkage groups",
        ylab = "Counts",
        main = "Ploidy patterns across isolates")

abline(h=0, col="gray", lty=1)

dev.off()

print("Finished plotting trends")
