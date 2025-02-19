#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#	Bedmethyl.R
#
#	Purpose 
#		This code was made for calculating high confidence calls 
#
# Author comment:
#	Lisa Koole
#	lisa.koole@maastrichtuniversity.nl
#
#

#-----------------------------------------------------------------------------------------------------#
# 							1. Main settings                                     ----
#-----------------------------------------------------------------------------------------------------#

# Get all general settings 
s_ROOT_dir <<- "G:/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/"
source(paste0(s_ROOT_dir,"Scripts/.Main/Settings.R"))


# library(profvis)
# profvis({
#   data(diamonds, package = "ggplot2")
#   
#   plot(price ~ carat, data = diamonds)
#   m <- lm(price ~ carat, data = diamonds)
#   abline(m, col = "red")
# })

#-----------------------------------------------------------------------------------------------------#
# 							2. Load libraries needed                             ----
#-----------------------------------------------------------------------------------------------------#
library(dplyr)
library(data.table)
library(purrr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(readxl)
library(annotatr)
library(envnames)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(pbapply)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(vcfR)
library(corrplot)
library(Hmisc)
library(matrixStats)
library(rstatix)
library(edgeR)
library(arsenal)
library(SNPfiltR)

devtools::install_github("DevonDeRaad/SNPfiltR")
#-----------------------------------------------------------------------------------------------------#
# 							3. Process bedmethyl data                            ----
#-----------------------------------------------------------------------------------------------------#

# Create list of bedmethyl files to import
temp = list.files(s_bedfiles_folder, pattern="\\.bed$", full.names = TRUE)
cat(paste0("number of samples: " , length(temp)))

pblapply(temp, function(x) {

  # Load bedmethyl file and select relevant columns
  bedfile <- tryCatch({fread(x, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="",
                   colClasses =c("character", "integer", "integer", "character", "numeric",
                                 "character", "NULL", "NULL", "NULL", "NULL", "numeric",
                                 "NULL", "integer", "NULL", "NULL", "NULL", "NULL", "NULL"))
  }, error = function(e) {
    warning(paste("Error reading file:", x, "-", e$message))
    return(NULL)
  })
  
  # Rename columns of interest
  bedfile <- setNames(bedfile, c("chrom", "start_pos0", "end_pos0", "mod_base_code","score",
                        "strand", "percent_mod", "N_canon"))

  # Extract the base name of the file
  file_base_name <- basename(x)
 
  # Extract postions
  position <- paste0(bedfile$chrom, "_", as.character(bedfile$start_pos0))
  
  # Extract 5mC values
  rows_with_m <- bedfile$mod_base_code== "m"
  methyl <- data.frame("position" = position[rows_with_m], "percent_mod" = bedfile$percent_mod[rows_with_m])

  # Extract 5hmC values
  rows_with_h <- bedfile$mod_base_code== "h"
  hydroxymethyl <- data.frame("position" = position[rows_with_h], "percent_mod" = bedfile$percent_mod[rows_with_h])

  # Extract canonical C values
  canon_perc <- (bedfile$N_canon/bedfile$score)*100
  canonical <- unique(data.frame("position" = position, canon_perc))

  # Extract combined 5mC and 5hmC values
  modified_perc <- 100 - (bedfile$N_canon/bedfile$score)*100
  modified <- unique(data.frame("position" = position, modified_perc))

  # Extract scores (number of calls)
  scores <- unique(data.frame("position" = position, "scores" = as.numeric(bedfile$score)))

  # Save processed data

  tryCatch(saveRDS(bedfile, file = file.path(s_OUT_dir, "Bedmethyl/", paste0(file_base_name, ".rds"))), error = function(e) warning(e$message))
  tryCatch(saveRDS(methyl, file = file.path(s_OUT_dir, "Bedmethyl/", paste0(file_base_name, "_methyl.rds"))), error = function(e) warning(e$message))
  tryCatch(saveRDS(hydroxymethyl, file = file.path(s_OUT_dir, "Bedmethyl/", paste0(file_base_name, "_hydroxymethyl.rds"))), error = function(e) warning(e$message))
  tryCatch(saveRDS(canonical, file = file.path(s_OUT_dir, "Bedmethyl/", paste0(file_base_name, "_canonical.rds"))), error = function(e) warning(e$message))
  tryCatch(saveRDS(modified, file = file.path(s_OUT_dir, "Bedmethyl/", paste0(file_base_name, "_modified.rds"))), error = function(e) warning(e$message))
  tryCatch(saveRDS(scores, file = file.path(s_OUT_dir, "Bedmethyl/", paste0(file_base_name, "_scores.rds"))), error = function(e) warning(e$message))

  # Free up memory
  rm(bedfile, methyl, hydroxymethyl, canonical, modified, scores)
  gc()

})


#-----------------------------------------------------------------------------------------------------#
# 						4. Load pheno data and determine CaseIDs                          ----
#-----------------------------------------------------------------------------------------------------#

# Load pheno data
load(paste0(s_ROOT_dir,s_out_folder,"pheno/Pheno.Rdata")) 

# Create temporary list with saved bedmethyl files
temp = list.files(file.path(s_OUT_dir, "Bedmethyl/"), pattern=paste0("\\.bam.bed.rds$"), 
                  full.names = TRUE)

# Extract the file basename (Sample Name)
file_base_name <- sub(".bam.bed.*", "", basename(temp))

# Match pheno information to samples included in data set
samples <- data.frame("Barcode" = sub(".*SQK-NBD114-24_","", file_base_name),
                      "Flowcell" = sub("_filtered_SQK-NBD114-24.*", "", file_base_name))


# Check for values that do not match 
cat("Values in pheno but not in samples: \n", setdiff(pheno$Flowcell, samples$Flowcell))  # Values in pheno not in samples
cat("Values in samples but not in pheno: \n", setdiff(samples$Flowcell, pheno$Flowcell))  # Values in samples not in pheno

# Set sample names
pheno_subset <- inner_join(samples, pheno, by = c("Flowcell", "Barcode")) 
sample_names <- paste0("sample_", pheno_subset$CaseID)
pheno_subset <- cbind(pheno_subset, sample_names)

# Free up memory
rm(pheno, sample_names)
gc()

#-----------------------------------------------------------------------------------------------------#
#					  5. Create dataframe for 5mC, 5hmC and canonical values       ----
#-----------------------------------------------------------------------------------------------------#

# Function for creating dataframes with CpGs (rows) and samples (columns)
combine_perc <- function(mod_type) {

  temp = list.files(file.path(s_OUT_dir, "Bedmethyl/"), 
                    pattern=paste0("\\.bed_", mod_type, ".rds$"), 
                    full.names = TRUE)
  
  print(length(temp))
  
  bedfiles_all <- readRDS(temp[1]) # Start with first sample
  
  progress_bar <- txtProgressBar(min=0, max=length(temp)) # Add progress bar
  
  for (i in 2:length(temp)) {
    # Read BED file
    file <- readRDS(temp[i])
    file_base_name <- sub(".bam.bed.*", "", basename(temp[i]))
    
    # Combine methyl values
    bedfiles_all <- full_join(bedfiles_all, file, by = "position")
    colnames(bedfiles_all)[1] <- c("position")
    colnames(bedfiles_all)[i] <- sub(".bam.bed.*", "", basename(temp[i-1])) # existing data frame
    colnames(bedfiles_all)[i+1] <- sub(".bam.bed.*", "", basename(temp[i])) # new data frame
    
    rm(file)
    gc()
    
    setTxtProgressBar(progress_bar, value = i)
  }
  
  # Move position values to rownames
  rownames(bedfiles_all) <- bedfiles_all$position
  bedfiles_all <- bedfiles_all %>% dplyr::select(!(position))
  
  # Set sample names
  names(bedfiles_all) <- paste0("sample_", pheno_subset$CaseID)
  
  # save
  saveRDS(bedfiles_all, file = file.path(s_OUT_dir, "Bedmethyl/", paste0("bedfiles_", mod_type,"_all.rds")))
  
  gc()
}

# Use combine_perc function for all types of modifications
combine_perc("methyl")
combine_perc("hydroxymethyl")
combine_perc("canonical")
combine_perc("modified")
combine_perc("scores")


#-----------------------------------------------------------------------------------------------------#
#                       6. Record of failure during QC                     ----
#-----------------------------------------------------------------------------------------------------#

# SamplesFail will be the boolean record of which samples have already failed. 
# Create dataframe with quality metrics that can later be added to phenotype information.
QCmetrics <- pheno_subset 
SamplesFail <- as.logical(rep("FALSE", nrow(pheno_subset)))

# Entries will be changed to TRUE as samples fail 
Stepsummary <- as.data.frame(matrix(ncol=0, nrow=2)) 
rownames(Stepsummary) <- c("Failed This Step", "Total Failed")


#-----------------------------------------------------------------------------------------------------#
#					         7. Number of calls and coverage                  ----
#-----------------------------------------------------------------------------------------------------#

file = readRDS(file.path(s_OUT_dir, "Bedmethyl/bedfiles_scores_all.rds"))

file_trans <- round(log(file, base = 10), digits = 2)

par(mfrow = c(2,5))

destination <- file.path(s_ROOT_dir, s_out_folder, paste0("Plots/QC1_Coverages_boxplot.pdf"))
pdf(destination, width = 11.7, height = 8.3)

for (i in colnames(file)) {
  log.file <- log(file[,i], base = 10)
  boxplot(log.file,
          ylab = "log10(coverage)",
          xlab = colnames(file)[i],
          ylim = c(0, 3))
  rm(log.file)
  gc()
}

dev.off()

help(pdf)



  


# Function to calculate the total number of (high confidence) calls per sample
calc_number_calls <- function(scores_file, cutoff) {
  
  # Calculate total number of calls
  file = readRDS(scores_file)
  
  number_calls <- t(data.frame("total_calls" = colSums(!is.na(file))/1000000, 
                               "hc_calls" = apply(file, 2, function(x) sum(x >= cutoff, na.rm = TRUE)/1000000)))
  
  gc()
  
  # Plot number of calls 
  destination <- file.path(s_ROOT_dir, s_out_folder, paste0("Plots/QC1_NumberCalls.pdf"))
  pdf(destination)
  
  barplot(number_calls, names.arg=colnames(file),
          ylab = "Number of total calls (x 1,000,000)",
          xlab="CaseID",
          main=paste0("number of calls with >=", cutoff, "X"),
          horiz = FALSE, las=2, beside=TRUE, col = c("white", "blue4"))
  
  dev.off()
  
  # Save number of calls
  write.csv(as.data.frame(number_calls), file.path(s_ROOT_dir, s_out_folder, paste0("QC/number_calls.csv")))
  return(number_calls)
}


# Store number of calls
number_calls <- calc_number_calls(scores_file = file.path(s_OUT_dir, "Bedmethyl/bedfiles_scores_all.rds"), cutoff = 10)

# Add number of calls to QCmetrics
QCmetrics <- cbind(QCmetrics, "total_calls" =  number_calls[1,], "hc_calls" = number_calls[2,])

# Mark samples that fail the minimum of CpGs with >10X
min_cpgs <- 0.5
min_cpgs_filter <- pheno_subset$sample_names %in% colnames(number_calls)[number_calls[2,] < min_cpgs]

# Samples that do not pass threshold
QCmetrics<-cbind(QCmetrics,high_confidence_filter)
SamplesFail[which(high_confidence_filter==FALSE)]<-TRUE
Step6<-c(length(which(high_confidence_filter==FALSE)),sum(SamplesFail))
Stepsummary<-cbind(Stepsummary,Step6)
print(Stepsummary)

# free up memory 
rm(number_calls)
gc()


#-----------------------------------------------------------------------------------------------------#
#                      14. Removal of CpGs                      ----
#-----------------------------------------------------------------------------------------------------#

cutoff = 6 # set min read count

# Function to remove CpGs with overall low coverage
filter_cpgs_cpm <- function(scores_file, cutoff) {
  
  # Read file and impute zero in NA values (no reads)
  file <- readRDS(scores_file)
  file[is.na(file)] <- 0
  
  # determine positions with min CUTOFF counts per million in at least 70% of samples
  keep <- filterByExpr(file, min.count = cutoff, group = NULL, min.prop=0.7)
  positions <- rownames(file[keep,])
  print(summary(positions))
  gc()
  return(positions)
}

positions_filtered <- filter_cpgs_cpm(scores_file = file.path(paste0(s_OUT_dir, "Bedmethyl/bedfiles_scores_all.rds")), 
                                      cutoff = cutoff)

# Apply filter to all files (5mC, 5hmC, modified etc.)
temp = list.files(file.path(s_OUT_dir, "Bedmethyl/"), pattern="_all.rds", 
                  full.names = TRUE)

pblapply(temp, function(bedfile) {
  
  # Read file
  file <- readRDS(bedfile)
  cutoff = cutoff
  
  # Extract CpGs
  file_filtered <- file[rownames(file) %in% positions_filtered,]
  
  # Extract the base name of the file
  file_base_name <- sub("bedfiles_", "", sub("_all.rds*", "", basename(bedfile)))
  
  # Save processed data
  saveRDS(file_filtered, file = file.path(s_OUT_dir, "Bedmethyl/", paste0("bedfiles_", file_base_name, "_all_cmpfiltered_", cutoff, "X.rds")))
  
  # Free up memory
  rm(file)
  gc()
  
}, cl = 4)

gc()



#-----------------------------------------------------------------------------------------------------#
#					 8. Basic stats and Batch effects             ----
#-----------------------------------------------------------------------------------------------------#

# Create list of bedmethyl files
temp = list.files(file.path(s_OUT_dir, "Bedmethyl/"), pattern="bedfiles_.rds", 
                  full.names = TRUE)

# Calculate basic stats
stats_all <- pblapply(temp, function(bedfile) {
  
  file = readRDS(bedfile)
  print(dim(file))
  file_base_name <- sub("bedfiles_", "", sub("_all.rds*", "", basename(bedfile)))
  
  # Calculate median and standard deviation per sample
  median_mod <- colMedians(as.matrix(file), na.rm = TRUE)
  sd_mod <- colSds(as.matrix(file), na.rm = TRUE)
  
  stats <- data.frame(median_mod, sd_mod)
  colnames(stats) <- c(paste0("median_", file_base_name), paste0("sd_", file_base_name))
  
  # Calculate median per batch (flowcell)
  batch_median <- aggregate(median_mod, by = list(pheno_subset$Flowcell), FUN = median)
  
  # Plot each flowcell / cohort as a boxplot
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC1.BoxplotFlowcellCohort_", file_base_name, ".pdf", sep="")
  pdf(file=destination)
  
  nCol<-length(unique(pheno_subset$Flowcell)) ## assumes there is a column called Plate in your phenotype file
  boxplot(median_mod ~ pheno_subset$Flowcell, ylab = paste0("Median ",file_base_name, " (%)"), 
          xlab = "Flowcell", las = 2, horizontal = FALSE ,col = rainbow(nCol)) 
  boxplot(median_mod ~ pheno_subset$Cohort, ylab = paste0("Median ",file_base_name, " (%)"), 
          xlab = "Cohort", las = 2, horizontal = FALSE, col = rainbow(nCol)) 
  
  dev.off()
  
  # Free memory
  rm(file, median_mod, sd_mod, batch_median)
  gc()

  return(stats)
})

# Combine list of data frames into one data frame
stats_all <- bind_cols(stats_all)

write.csv(as.data.frame(stats_all), file.path(s_ROOT_dir, s_out_folder, paste0("QC/stats_all.csv")))

QCmetrics <- cbind(QCmetrics, stats_all)

# Alternatively color points in original scatterplot by flowcell 
# (relevant if more variation between samples)
destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC1.ScatterplotFlowcell.pdf", sep="")
pdf(file=destination)

nCol<-length(unique(pheno_subset$Flowcell))
plot(stats_all$median_methyl, stats_all$median_hydroxymethyl, pch = 16, xlab = "Median 5mC (%)", ylab = "Median 5hmC (%)", col = rainbow(nCol)[factor(pheno_subset$Flowcell)])
legend("topright", levels(factor(pheno_subset$Flowcell)), col = rainbow(nCol), pch = 16, cex = 0.4)

dev.off()


# #-----------------------------------------------------------------------------------------------------#
# #                       12. Filter sample outliers          -----
# #-----------------------------------------------------------------------------------------------------#

detect_outliers <- function(stats_df) {
  
  # plot boxplots
  stats_df_median <- stats_df %>% select(contains("median"))
  
  boxplot(stats_df_median, las = 2)
  hist(stats_df_median)
  
  # Determine outliers in 5-hmC and 5-mC
  out_hm <- identify_outliers(stats_df, variable = "median_hydroxymethyl")
  out_m <- identify_outliers(stats_df, variable = "median_methyl")
  
  outliers <- intersect(rownames(out_hm), rownames(out_m)) # in both 5mC and 5hmC 
  
  # outliers <- unique(rownames(out_hm), rownames(out_m)) # in either 5mC or 5hmC
  
  # plot post outlier
  boxplot(stats_df_median[!rownames(stats_df_median) %in% outliers,], las = 2)
  hist(stats_df_median[!rownames(stats_df_median) %in% outliers,])
  
  return(outliers)
}

outliers <- detect_outliers(stats_all)
dev.off()



# Mark samples that are outliers
outliers_filter <- rownames(stats_all) %in% outliers

# Samples that do not pass threshold
QCmetrics<-cbind(QCmetrics,outliers_filter)
SamplesFail[which(outliers_filter==FALSE)]<-TRUE
Step6<-c(length(which(outliers_filter==FALSE)),sum(SamplesFail))
Stepsummary<-cbind(Stepsummary,Step6)
print(Stepsummary)

# free up memory
rm (stats_all)
gc()


#-----------------------------------------------------------------------------------------------------#
#                    9. Check sex                       ----
#-----------------------------------------------------------------------------------------------------#

# Filter for sex chromosomes
filter_sexchr <- function(scores_file) {
  positions <- data.frame(rownames(readRDS(scores_file)))
  colnames(positions) <- "position"
  positions_sex <- positions %>% filter(str_detect(position, "chrX_|chrY_"))
  return(positions_sex$position)
}
  
positions_sexchr <- filter_sexchr(scores_file = file.path(paste0(s_OUT_dir, "Bedmethyl/bedfiles_scores_all.rds")))


# Load function findGenderPC
findGenderPC <- function(modifications, sex, npcs = 4, file_base_name){
  
  # Select no outlier cases 
  modifications[!rownames(modifications) %in% outliers,]
  
  # Select complete cases 
  modifications_all_com<-modifications[complete.cases(modifications),]
  
  # Select positions on the sex chromosomes
  modifications_all_com <- modifications_all_com[rownames(modifications_all_com) %in% positions_sexchr,]
  pca<-prcomp(modifications_all_com)
  
  # Correlate prinicpal components with factor Gender
  pca.cor <- rep(NA, npcs)
  
  for(i in 1:npcs){
    pca.cor[i]<-cor(pca$rotation[,i], as.numeric(as.factor(sex)), use = "complete")  
    }
  
  top<-order(abs(pca.cor), decreasing = TRUE)[1]
  second<-order(abs(pca.cor), decreasing = TRUE)[2]
  print(paste("Top correlated principal components with sex:", top, ",", second))
  
  # Plot the top principal components, color by sex 
  destination <- file.path(s_ROOT_dir, s_out_folder, paste0("Plots/QC2_PredSex_", file_base_name, ".pdf"))
  pdf(destination)
  
  plot(pca$rotation[,top], pca$rotation[,second], pch = 16, col = c("green", "darkblue")[as.factor(sex)],
    xlab = paste("PC", top), ylab = paste("PC", second))
  legend("topright", levels(as.factor(sex)), pch = 16, col = c("green", "darkblue"))
  
  dev.off()
  
  # Predict sex
  predSex <- rep(NA, length(sex))
  options.sex<-levels(as.factor(sex))
  
  if(abs(pca.cor[top]) > 0.9){ print("Top PC has r > 0.9 with sex so good enough to confirm reported sexes") } else
  {print(paste("Top PC has r =", round(abs(pca.cor[top]),2), "with sex so may not be good enough to confirm reported sexes")) }
  
  if(sign(pca.cor[top]) == 1){
    predSex[which(pca$rotation[,top] < 0)] <- options.sex[1]
    predSex[which(pca$rotation[,top] > 0)] <- options.sex[2]} 
  
  else {
   predSex[which(pca$rotation[,top] < 0)] <-options.sex[2]
   predSex[which(pca$rotation[,top] > 0)] <-options.sex[1]} 
  
  predSex <- data.frame(predSex) 
  colnames (predSex) <- paste0("predSex_", file_base_name)
  gc()
  return(predSex)
}

# Use findGenderPC function to predict sex 
ReportedSex <- pheno_subset$Gender

temp = list.files(file.path(s_OUT_dir, "Bedmethyl/"), pattern="_all.rds", 
                  full.names = TRUE)

PredictedSex <- pblapply(temp, function(bedfile) {
  # Read RDS file
  file = readRDS(bedfile)
  file_base_name <- sub("bedfiles_", "", sub("_all.rds*", "", basename(bedfile)))
  
  # Predicted sex
  predSex <- findGenderPC(file,ReportedSex, npc=4, file_base_name) # npcs = 16
  return(predSex)
})

# Add predicted sex to QCmetrics
PredictedSex <- bind_cols(PredictedSex)
QCmetrics <- cbind(QCmetrics, PredictedSex)

# Update stepssummary
SamplesFail[which(PredictedSex!=ReportedSex)] <- TRUE 
Step3<-c(length(which(PredictedSex!=ReportedSex)),sum(SamplesFail)) 
Stepsummary<-cbind(Stepsummary,Step3)
print(Stepsummary)

# Free up memory
rm(PredictedSex, positions_sexchr)
gc()

#-----------------------------------------------------------------------------------------------------#
#                         10. Check SNPs                ----
#-----------------------------------------------------------------------------------------------------#

# Create list of VCF files
temp = list.files(s_vcfs_folder, pattern="\\.vcf$", full.names = TRUE)

pblapply(temp, function(x) {

  # Load VCF file
  vcf <- read.vcfR(x)
  
  # Extract DP coverage values
  dp.matrix <- vcfR::extract.info(vcf, element = "DP", 
                                  as.numeric = TRUE)
  
  # Extract values for which the depth coverage > 50
  vcf@fix <- vcf@fix[dp.matrix > 50,]
  vcf@gt <- vcf@gt[dp.matrix > 50,]

  # Extract alternative nucleotide
  vcf_snps_num <- t(vcfR2loci(vcf, return.alleles=TRUE))

  vcf_snps_num <- data.frame("position" = rownames(vcf_snps_num), "alleles" = vcf_snps_num)

  # Extract the base name of the file
  file_base_name <- basename(x)

  # Save processed data
  saveRDS(vcf_snps_num, file = file.path(s_OUT_dir, "VCF/", paste0(file_base_name, "_alleles.rds")))

  # Free up memory
  rm(vcf, vcf_snps_num)
  gc()

}, cl = 4)


# # Example for one VCF file
# temp = list.files(s_vcfs_folder, pattern="\\.vcf$", full.names = TRUE)
# 
# vcf <- read.vcfR(temp[1])
# 
# dp.matrix <- vcfR::extract.info(vcf, element = "DP", 
#                               as.numeric = TRUE)
# 
# # Extract values for which the depth coverage > 50
# vcf@fix <- vcf@fix[dp.matrix > 50,]
# vcf@gt <- vcf@gt[dp.matrix > 50,]
# 
# 
# t(vcfR2loci(vcf, return.alleles=TRUE))

# Com# Com# Combine the VCF files
temp = list.files(file.path(s_OUT_dir, "VCF/"),
                    pattern=paste0("\\.rds$"),
                    full.names = TRUE)

vcf_all <-readRDS(temp[1]) # Start with first sample

progress_bar <- txtProgressBar(min=0, max=length(temp)) # Add progress bar

for (i in 2:length(temp)) {
  # Read BED file
  file <- data.frame(readRDS(temp[i]))
  file_base_name <- sub(".bam_igv.vcf_scores*", "", basename(temp[i]))

  # Combine methyl values
  vcf_all <- merge(vcf_all, file, by = "position", all=TRUE)
  colnames(vcf_all)[1] <- c("position")
  colnames(vcf_all)[i] <- sub(".bam_igv.vcf_alleles.rds*", "", basename(temp[i-1])) # existing data frame
  colnames(vcf_all)[i+1] <- sub(".bam_igv.vcf_alleles.rds*", "", basename(temp[i])) # new data frame

  # Free memory
  rm(file)
  gc()

  setTxtProgressBar(progress_bar, value = i)
}

# Move position values to rownames
rownames(vcf_all) <- vcf_all$position
vcf_all <- vcf_all %>% dplyr::select(!(position))

# Set sample names
# names(vcf_all) <- paste0("sample_", pheno_subset$CaseID)
# colnames(vcf_all_factor) <- paste0("sample_", seq.int(names(vcf_all)))

# select complete cases
vcf_all_complete <- vcf_all[complete.cases(vcf_all),]
head(vcf_all_complete)

# Factorize the rows based on genotypes (e.g. 0, 1,2)
vcf_all_factor <- t(apply(vcf_all_complete, 1, function(row) as.numeric(factor(row))))

# save
saveRDS(vcf_all_complete, file = file.path(s_OUT_dir, "VCF/", paste0("vcf_all_complete.rds")))
saveRDS(vcf_all_factor, file = file.path(s_OUT_dir, "VCF/", paste0("vcf_all_factor.rds")))


#-----------------------------------------------------------------------------------------------------#
#      11. Check genetically identical samples correlate across SNP probes        ----
#-----------------------------------------------------------------------------------------------------#

# Load VCF dataframes
vcf_all_complete <- readRDS(file = file.path(s_OUT_dir, "VCF/", paste0("vcf_all_complete.rds")))
vcf_all_factor <- readRDS(file = file.path(s_OUT_dir, "VCF/", paste0("vcf_all_factor.rds")))

# ~ Genetic correlations -------------------------------------------
snpCor<-cor(vcf_all_factor, use="pairwise.complete.obs")
round(snpCor,3)

rcorr(vcf_all_factor, type = c("pearson","spearman")) # with p-values

#  ~  Create plots for correlation ---------------------------------

# Plot dendrogram
hclust_res <- hclust(dist(t(vcf_all_factor)), method = "complete")
plot(hclust_res)

# Visualize correlation matrix (with numbers)
corrplot(snpCor, type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45, method = 'number', addrect = 3,  col = COL2('RdBu', 20), col.lim=c(0,1))

# Visualize correlation matrix (with circles)
corrplot(snpCor, type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45, addrect = 3,  col = COL2('RdBu', 20), col.lim=c(0,1))

# Visualize correlation matrix in heatmap
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = snpCor, col = col, symm = TRUE)

heatmap_vcfs <- function(factorized_alleles) {
  # Get overview of rows with identical values across columns
  rows_no_diff <- apply(vcf_all_factor, 1, function(row){
    all(row == row[1])
  })

  print(summary(rows_no_diff))

  # Get vcf_genotypes with differences
  vcf_dif <- vcf_all_factor[rows_no_diff == FALSE,]

  # Make heatmap
  heatmap(vcf_dif)
}

heatmap_vcfs(vcf_all_factor)


# ~ Calculate Max correlation -----------------------------
# Ignore correlations between a sample and itself
for(i in 1:ncol(vcf_all_factor)){
  snpCor[i,i]<-NA
}

corMax<-apply(snpCor, 1, max, na.rm = TRUE)
hist(corMax, xlab = "Max. correlation with all other samples", main = "")
corMax

# Add corMax to QCmetrics
QCmetrics<-cbind(QCmetrics, corMax)
SamplesFail[which(corMax>0.8)]<-TRUE
Step5<-c(sum(corMax>0.8),sum(SamplesFail))
Stepsummary<-cbind(Stepsummary,Step5)
print(Stepsummary)





#-----------------------------------------------------------------------------------------------------#
#                 13. Save QC Results of the failed and passed samples   -----
#-----------------------------------------------------------------------------------------------------#

write.csv(QCmetrics[SamplesFail,], paste0(s_ROOT_dir, s_out_folder, "QC/SamplesFailedQC.csv"), row.names = FALSE)
write.csv(QCmetrics[!SamplesFail,], paste0(s_ROOT_dir, s_out_folder, "QC/SamplesPassedQC.csv"), row.names = FALSE)





#-----------------------------------------------------------------------------------------------------#
#                     15. Create matrix of weights          ----
#-----------------------------------------------------------------------------------------------------#

# ~ weights based on all samples  --------------------------

weights_single_sample <- function(scores_file) {
  
  scores_filtered <- readRDS(scores_file)

  # Compute weights using vectorized row-wise division
  scores_weights <- sweep(scores_filtered, 1, rowSums(scores_filtered), `/`)
  print(head(scores_weights))
  saveRDS(scores_weights, file = paste0(s_ROOT_dir,s_out_folder,"QC/scores_global_weights.rds")) 
  
  rm(scores_filtered, scores_weights)
  
  gc()
}

weights_single_sample(scores_file = file.path(s_OUT_dir, "Bedmethyl/bedfiles_scores_all_filtered.rds"))

# scores_weights <- apply(scores_filtered,1, function(x) matrixStats::colWeightedMeans(matrix(x))) # does not work yet


# ~ weights based on all samples within group -----------------------------

weights_global_sample <- function(scores_file, cohorts) {
  
  # Prepare list of cohorts
  cohort_weights_list <- list()
  unique_cohorts <- cohorts 
  
  # Read scores file
  scores_filtered <- readRDS(scores_file)
  
  for (cohort in unique_cohorts) {
    # samples that belong to a specific cohort
    cohort_samples <- sample_names[pheno_subset$Cohort == cohort]
    
    # scores for these specific samples
    cohort_scores <- scores_filtered[,colnames(scores_filtered) %in% cohort_samples]
    
    # scores divided by total sum of scores
    cohort_weights <- sweep(cohort_scores, 1, rowSums(scores_filtered), `/`)
    
    # return weights
    cohort_weights_list[[cohort]] <- cohort_weights
  }
  
  # Combine cohort weights
  scores_weights <- do.call(cbind, cohort_weights_list)
  
  # remove prefix (AD, MCI, CTL)
  colnames(scores_weights) <- unlist(lapply(cohort_weights_list, function(x) colnames(x)))
  
  # Order weights based on Case ID
  scores_weights <- scores_weights[,colnames(scores_filtered)]
  print(head(scores_weights))
  
  # Save file
  saveRDS(scores_weights, file = paste0(s_ROOT_dir,s_out_folder,"QC/scores_cohort_weights2.rds")) 
  
  # Free up memory
  rm(scores_weights,cohort_scores, cohort_weights, scores_filtered)
  gc()
}

weights_global_sample(scores_file = file.path(s_OUT_dir, "Bedmethyl/bedfiles_scores_all_filtered.rds"), 
                      cohorts = unique(pheno_subset$Cohort))


#-----------------------------------------------------------------------------------------------------#
#					        16. Density plots per cohort                ----
#-----------------------------------------------------------------------------------------------------#

temp = list.files(file.path(s_OUT_dir, "Bedmethyl/"), pattern="all_filtered.rds", 
                  full.names = TRUE)

pblapply(temp, function(bedfile) {
  
  # Read RDS file
  file = readRDS(bedfile)
  file_base_name <- sub("_all_filtered.rds*", "", basename(bedfile))
  
  # Data per cohort
  dat_mod <- reshape2::melt(file)
  colnames(dat_mod) <- c("sample_names", "value")
  dat_mod <- dat_mod %>% left_join(pheno_subset[,c("sample_names","Cohort")], by = join_by(sample_names))
  gc()
  
  # Basic density plots, separate graphs per group
  images <- ggplot(dat_mod, aes(x=value, fill = Cohort)) + 
    geom_density(alpha = 0.8) + scale_fill_manual(values=c("#505567", "#E5E6EB", "blue4"))+ 
    facet_grid(Cohort ~ .) + xlab("Beta value") +
    ylab("Density") +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15, face ="bold"))
  
  ggexport(images,  width = 2000, height = 2500, res = 200,
           filename = paste0(s_ROOT_dir, s_out_folder, "Plots/BetasDensities_separate_", file_base_name,".png"))
  
  # Basic density plot; groups overlapping
  images <- ggplot(dat_mod, aes(x=value, fill = Cohort)) + 
    geom_density(alpha = 0.8) + 
    scale_fill_manual(values=c("#505567", "#E5E6EB", "blue4"), labels = c("AD", "CTL", "MCI")) +
    xlab("Beta value") +
    ylab("Density") + theme_classic() +
    theme(
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15, face ="bold")
    ) 
  
  ggexport(images, width = 3000, height = 1500, res = 200,
           filename = paste0(s_ROOT_dir, s_out_folder, "Plots/BetasDensities", file_base_name, ".png"))
  
  # Free up memory
  rm(images, dat_mod, file)
  gc()
})


#-----------------------------------------------------------------------------------------------------#
#					      17. Density plots per modification type           ------
#-----------------------------------------------------------------------------------------------------#

#### MAKE IT INTERACTIVE?

# List of samples 
temp = list.files(file.path(s_OUT_dir, "Bedmethyl/"), pattern="bam.bed.rds", 
                  full.names = TRUE)

file_base_name <- sub(".bam.bed.rds*", "", basename(temp))

# ~ Density plots per sample -------------------
destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC10.ModificationDensities.pdf", sep="")
pdf(destination)

pblapply(file_base_name, function(file) {
  
  # Retrieve dataframes
  canonical <- readRDS(file.path(s_OUT_dir, "Bedmethyl/", paste0(file, ".bam.bed_canonical.rds")))
  methyl <- readRDS(file.path(s_OUT_dir, "Bedmethyl/", paste0(file, ".bam.bed_methyl.rds")))
  hydroxymethyl <- readRDS(file.path(s_OUT_dir, "Bedmethyl/", paste0(file, ".bam.bed_hydroxymethyl.rds")))
  
  # Plot the densities 
  plot(density(canonical[,2], na.rm=TRUE), main = c("Density C modifications (%)", file, sep = ""), col = "darkgrey") 
  lines(density(methyl[,2], na.rm=TRUE), col = "darkblue")
  lines(density(hydroxymethyl[,2], na.rm=TRUE), col = "orange")
  legend("topright", legend=c("canonical C", "5mC", "5hmC"), lty=1, col=c("darkgrey","darkblue", "orange"))
  
  # Clear memory
  rm(canonical, methyl, hydroxymethyl)
  gc()
  
})

dev.off()


# ~ Density plots on all data --------------------------------------

plot_densities <- function(mod_type) {
  
  # Get colors based on cohort
  colors <- as.character(factor(pheno_subset$Cohort, labels = c("forestgreen","#1a80bb", "darkorange")))
  
  # Retrieve data
  file <- readRDS(file.path(paste0(s_OUT_dir, "Bedmethyl/", "bedfiles_", mod_type, "_all_filtered.rds")))
  
  # Calculate max density 
  max_density <- max(apply(file, 2, function(x) {
    max(density(x)$y)
    }))
  
  # Plot densities
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/ModificationDensities_", mod_type, ".pdf")
  pdf(destination)
  
  plot(density(as.matrix(file)[,1], na.rm=TRUE), 
       main = c("Density C modifications (%)", mod_type, sep = ""), col = colors[1],
       ylim=c(0, max_density)) 
  
  for (i in 2:(length(file))) { lines(density(file[,i], na.rm=TRUE), col = colors[i])}
  
  legend("topright", legend=c("CTL", "MCI", "AD"), lty=1, col=c("forestgreen","#1a80bb", "darkorange"))
  
  dev.off()
  
  # Free up memory
  rm(file, max_density)
  gc()
}

plot_densities("methyl")
plot_densities("hydroxymethyl")
plot_densities("canonical")



# ~ Density plots on high confidence data --------------------------------------

plot_densities_hc <- function(mod_type) {
  
  # Get colors based on cohort
  colors <- as.character(factor(pheno_subset$Cohort, labels = c("forestgreen","#1a80bb", "darkorange")))
  
  # Retrieve data
  file <- readRDS(file.path(paste0(s_OUT_dir, "Bedmethyl/", "bedfiles_", mod_type, "_highconfidence.rds")))

  # Calculate max density 
  max_density <- max(apply(file, 2, function(x) {
    max(density(x, na.rm = TRUE)$y)
  }))
  
  # Plot densities
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/ModificationDensities_HC_", mod_type, ".pdf")
  pdf(destination)
  
  plot(density(as.matrix(file)[,1], na.rm=TRUE), 
       main = c("Density C modifications (%)", mod_type, sep = ""), col = colors[1],
       ylim=c(0, max_density)) 
  
  for (i in 2:(length(file))) { lines(density(file[,i], na.rm=TRUE), col = colors[i])}
  
  legend("topright", legend=c("CTL", "MCI", "AD"), lty=1, col=c("forestgreen","#1a80bb", "darkorange"))
  
  dev.off()
  
  # Free up memory
  rm(file, max_density)
  gc()
}

plot_densities_hc("methyl")
plot_densities_hc("hydroxymethyl")
plot_densities_hc("canonical")


#-----------------------------------------------------------------------------------------------------#
#                       18. Check PCA correlations (first 10 comps)         ----
#-----------------------------------------------------------------------------------------------------#

pca_correlations <- function(mod, filesuffix) {
  
  mod <- as.matrix(mod)
  # Recalculating PCA as the betas were updated 
  pca <- prcomp(t(mod[complete.cases(mod),]))
  
  # Plot some 
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/", "QC11.ScreePlot_", filesuffix, ".pdf", sep="")
  pdf(file=destination) 
  
  plot(pca) # scree plot, indicating which PC has most variance. 
  dev.off()
  
  Rel_variance_expl = round(pca$sdev^2 / sum(pca$sdev^2),3) # relative variance 
  cat("First component explains",Rel_variance_expl[1]*100, "% of variance in the data")
  
  # Transform the pheno data to a numeric matrix, including character columns. Should usually take a lot of attention, as some columns do 
  # not make sense to undergo this transformation. This is a lazy workaround. For this example its OK.
  pheno_numeric = apply(pheno_subset,2,function(x){as.numeric(as.factor(x))}) 
  pheno_numeric = pheno_numeric[,which(apply(pheno_numeric,2,function(x){sum(is.na(x))==0}))]
  
  # Correlate the PCs (max 10) to the phenodata, to see what effects can be explained by variance alone. Ideally, PC 1 is 
  # the intended experimental effect.
  correlation_frame = as.data.frame(na.omit(round(t(data.frame(cor(pca$x[,1:min(10,dim(mod)[2])],pheno_numeric))),2))) 
  correlation_frame[abs(correlation_frame)<0.3] = 0
  
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/", "QC11_PCA_", filesuffix, ".pdf", sep="")
  pdf(file=destination) 
  
  report_pc = 1:5 
  for(PC in report_pc){ # State what is correlating strongly to what 
    cat("\n",rownames(correlation_frame)[which.max(abs(correlation_frame[,PC]))],"(",which.max(abs(correlation_frame[,PC])),") is correlated 
        the most with PC",PC,"(",correlation_frame[which.max(abs(correlation_frame[,PC])),PC],")\n\n")
    
    # make colors based on factors
    lazycolors = as.numeric(as.factor(pheno_subset[,which.max(abs(correlation_frame[,PC]))]))
    
    # plot title
    main_title = paste0("PC ",PC," - ",rownames(correlation_frame)[which.max(abs(correlation_frame[,PC]))],"(",correlation_frame[which.max(abs(correlation_frame[,PC])),PC],")" )
    
    # plot pc score and color by correlate
    if(is.numeric(pheno_subset[,which.max(abs(correlation_frame[,PC]))])){
      plot(pca$x[,PC]~pheno_subset[,which.max(abs(correlation_frame[,PC]))],
           col=lazycolors,
           pch=19,
           xlab = rownames(correlation_frame)[which.max(abs(correlation_frame[,PC]))],
           ylab = paste0("PC ",PC," (",Rel_variance_expl[PC]*100,"%)"),
           main=main_title)
    }else{
      boxplot(pca$x[,PC]~pheno_subset[,which.max(abs(correlation_frame[,PC]))],
              col=lazycolors,
              pch=19,
              xlab = rownames(correlation_frame)[which.max(abs(correlation_frame[,PC]))],
              ylab = paste0("PC ",PC," (",Rel_variance_expl[PC]*100,"%)"),
              main=main_title)
    }
  }
  
  # plot PC1 and PC2, colours indicating cohort variable
  plot(pca$x[,2] ~ pca$x[,1], col=factor(pheno_subset$Cohort),
       xlab = paste0("PC", 1, " (", Rel_variance_expl[1]*100,"%)"),
       ylab = paste0("PC", 2, " (", Rel_variance_expl[2]*100,"%)"))
  
  pca_cohort <- cbind(as.data.frame(pca$x), "Cohort" = pheno_subset$Cohort)
  
  images <- as.list(1:2)
  
  images[[1]] <- ggplot(data = pca_cohort, mapping = aes(x = PC1, y = PC2, color = Cohort)) + 
    geom_point(shape = 10) +
    scale_color_manual(values=c("#505567","grey","blue"), 
                       name = "Cohort", 
                       labels = c("AD", "CTL", "MCI")) +
    xlab(paste0("PC", 1, " (", Rel_variance_expl[1]*100,"%)"))+
    ylab(paste0("PC", 2, " (", Rel_variance_expl[2]*100,"%)"))+
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15, face ="bold")
    ) 
  
  # plot PC1 and PC2, colours indicating gender variable
  pca_gender <- cbind(as.data.frame(pca$x), "Gender" = pheno_subset$Gender)
  
  images[[2]] <- ggplot(data = pca_gender, mapping = aes(x = PC1, y = PC2, color = Gender)) + 
    geom_point(shape = 10) +
    scale_color_manual(values=c("#505567","grey"), 
                       name = "Gender", 
                       labels = c("Female", "Male")) +
    xlab(paste0("PC", 1, " (", Rel_variance_expl[1]*100,"%)"))+
    ylab(paste0("PC", 2, " (", Rel_variance_expl[2]*100,"%)"))+
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15, face ="bold")
    ) 
  return(images)
  dev.off()
}

temp = list.files(file.path(s_OUT_dir, "Bedmethyl/"), pattern="_all.rds", 
                  full.names = TRUE)


pblapply(temp, function(bedfile) {
  
  # Read RDS file
  file = readRDS(bedfile)
  file_base_name <- sub("bedfiles_", "", sub("_all.rds*", "", basename(bedfile)))
  
  # Predicted sex
  pcaplots <- pca_correlations(file, file_base_name) # npcs = 16
  
  dev.off()
  
})



temp = list.files(file.path(s_OUT_dir, "Bedmethyl/"), pattern="_highconfidence.rds", 
                  full.names = TRUE)


pblapply(temp, function(bedfile) {
  
  # Read RDS file
  file = readRDS(bedfile)
  file_base_name <- sub("bedfiles_", "", basename(bedfile))
  
  # Predicted sex
  pcaplots <- pca_correlations(file, file_base_name) # npcs = 16
  
  dev.off()
})




#-----------------------------------------------------------------------------------------------------#
#					        	    19. M values                      ----
#-----------------------------------------------------------------------------------------------------#

# Modify modification values of zero and 100 

imputezerohundred <- function(mod) {
  # Determine the second smallest number to zero 
  min_value <- min(mod[mod > 0])
  
  # Determine the second largest number to 100
  max_value <- max(mod[mod < 100])

  # Find the smallest of the min_value and max_value by two. 
  # This number will be used to set the min_value and max_value at equal distances from zero and 100, respectively
  step_value <- min(min_value - 0, 100 - max_value)
  min_value <- 0 + step_value
  max_value <- 100 - step_value
  
  # Impute the min_value and max_value in the values with zero and 100, respectively 
  mod[mod==0] <- min_value
  mod[mod==100] <- max_value
  
  cat("min value =", min_value, " ") 
  cat("max_value =", max_value)
  
  return(mod)
}

# Make M values from modification values 
# (more normally distributed per CpG) -- (methylated value / unmethylated value)
mod2M <- function(mod){
  return(log2((mod) /(100 - mod)))
  } 

# Create list with 
temp = list.files(file.path(s_OUT_dir, "Bedmethyl/"), pattern="_all_filtered.rds", 
                  full.names = TRUE)


pblapply(temp, function(bedfile) {
  
  # Read RDS file
  file = readRDS(bedfile)
  file_base_name <- sub("_all_filtered.rds*", "", file)
  
  # M values
  Metas = mod2M(as.data.frame(imputezerohundred(file)))
  
  # Save processed data
  saveRDS(Metas, file = file.path(s_OUT_dir, "QC/", paste0(file_base_name, "_Metas.rds")))
  
  # Histogram of M values
  destination <- file.path(paste0(s_ROOT_dir, s_out_folder, "Plots/", "QC12_M_values", file_base_name,".pdf", sep=""))
  pdf(file=destination) 
  
  hist(as.matrix(Metas))
  dev.off()
  
  # which(is.infinite(hydroxymethyl_Metas),arr.ind = FALSE)
  
  # Free up memory
  rm(file, Metas)
  gc()
  
})

dev.off()

# # should result empty now
# which(is.infinite(hydroxymethyl_Metas),arr.ind = FALSE)
# which(is.infinite(methyl_Metas),arr.ind = TRUE)
# which(is.infinite(unmod_Metas))


#-----------------------------------------------------------------------------------------------------#
#					        	20. Save relevant data              ----
#-----------------------------------------------------------------------------------------------------#

# Only keep the useful variables in pheno. 
pheno_new <- QCmetrics[!SamplesFail,] 
# pheno_new <- subset(pheno_new, select = -c(Intensity, M.median, U.median, PredictedSex, PredictedStrain))

saveRDS(pheno_new, file = paste0(s_ROOT_dir, s_out_folder,"QC/Pheno_new.rds"))
write.csv(pheno_new, file = paste0(s_ROOT_dir, s_out_folder,"QC/Pheno_new.csv"), row.names = FALSE)

# Save work space
save.image(file=paste0(s_ROOT_dir,s_out_folder,"EPI-Alzheimer-ONT.RData"))

