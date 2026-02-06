#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#	QC.R
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
# 		    1. Main settings                                                        ----
#-----------------------------------------------------------------------------------------------------#

# Get all general settings 
s_ROOT_dir <<- "path/to/root/directory"
source(paste0(s_ROOT_dir,"Scripts/.Main/Settings_v2.R"))


# library(profvis)
# profvis({
#   data(diamonds, package = "ggplot2")
#   
#   plot(price ~ carat, data = diamonds)
#   m <- lm(price ~ carat, data = diamonds)
#   abline(m, col = "red")
# })

#-----------------------------------------------------------------------------------------------------#
# 				2. Load libraries needed                                                ----
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
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(pbapply)
library(RColorBrewer)
library(scales)
library(ggpubr)
# library(vcfR)
library(corrplot)
library(Hmisc)
library(matrixStats)
library(rstatix)
library(edgeR)
library(arsenal)
library(RColorBrewer)
library(svglite)
library(viridis)


# library(devtools)
# devtools::install_github("DevonDeRaad/SNPfiltR")
# library(SNPfiltR)

#-----------------------------------------------------------------------------------------------------#
# 				3. Load pheno data and determine CaseIDs                              ----
#-----------------------------------------------------------------------------------------------------#

# Load pheno data
load(file.path(paste0(s_ROOT_dir,s_out_folder,"Pheno/Pheno.Rdata"))) 

# Create temporary list with saved bedmethyl files
temp = list.files(file.path(s_project_folder, "Bedmethyl/"), pattern=paste0("\\.bam.bed.rds$"), 
                  full.names = TRUE)

# Extract the file basename (Sample Name)
file_base_name <- sub(".bam.bed.*", "", basename(temp))

# Match pheno information to samples included in data set
samples <- data.frame("Barcode" = sub(".*barcode","barcode", file_base_name),
                      "Flowcell" = sub("_barcode.*", "", file_base_name))


# Check for values that do not match 
cat("Values in pheno but not in samples: \n", setdiff(pheno$Flowcell, samples$Flowcell))  # Values in pheno not in samples
cat("Values in samples but not in pheno: \n", setdiff(samples$Flowcell, pheno$Flowcell))  # Values in samples not in pheno



# Set sample names
pheno_subset <- inner_join(samples, pheno, by = c("Flowcell", "Barcode")) 
sample_names <- paste0("sample_", pheno_subset$CaseID)
pheno_subset <- cbind(pheno_subset, sample_names)

# Free up memory
rm(pheno, sample_names, samples)
gc()


#-----------------------------------------------------------------------------------------------------#
#         4. Record of failure during QC                                        ----
#-----------------------------------------------------------------------------------------------------#

# SamplesFail will be the boolean record of which samples have already failed. 
# Create dataframe with quality metrics that can later be added to phenotype information.
QCmetrics <- pheno_subset 
SamplesFail <- as.logical(rep("FALSE", nrow(pheno_subset)))

# Entries will be changed to TRUE as samples fail 
Stepsummary <- as.data.frame(matrix(ncol=0, nrow=2)) 
rownames(Stepsummary) <- c("Failed This Step", "Total Failed")


#-----------------------------------------------------------------------------------------------------#
#        5. Removal of CpGs                                                       ----
#-----------------------------------------------------------------------------------------------------#
# Set cutoff for coverage
cutoff = 5

# Function to remove CpGs with overall low coverage
filter_cpgs_cpm <- function(scores_file, cutoff) {
  
  # Read file and impute zero in NA values (no reads)
  file <- readRDS(scores_file)
  
  # determine positions with min CUTOFF counts per million in at least 70% of samples
  keep <- filterByExpr(file, min.count = cutoff, group = NULL, min.prop=0.7)
  gc()
  positions <- rownames(file[keep,])
  cat(print(summary(positions)))
  return(positions)
}

positions_filtered <- filter_cpgs_cpm(scores_file = file.path(paste0(s_ROOT_dir, s_out_folder, "Bedmethyl_all/bedfiles_scores_all.rds")), 
                                      cutoff = cutoff)

# keep <- rowSums(scores_file >=5)
# min.prop = 0.7 * nrow(pheno_subset)
# keep <- keep > min.prop

# Free up memory
gc()

# cutoff 5 -- 1 826 331
# cutoff 6 -- 279 418
# cutoff 7 -- 16 459 
# cutoff 10 -- ~7500

# cutoff 10 -- 8222
# cutoff 5 -- 2 418 654 
# cutoff 6 -- 632 297 

# cutoff 6 --- 118 525 (13062025)


# Apply filter to all files (5mC, 5hmC, modified etc.)
temp = list.files(file.path(s_ROOT_dir, s_out_folder, "Bedmethyl_all/"), pattern="_all.rds", 
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
  saveRDS(file_filtered, file = file.path(s_ROOT_dir, s_out_folder, "Bedmethyl_all/", paste0("bedfiles_", file_base_name, "_all_", cutoff, "X.rds")))
  
  # Free up memory
  rm(file)
  gc()
  
}, cl = 4)

rm(positions_filtered)
gc()


#-----------------------------------------------------------------------------------------------------#
#			  	6a. Coverage distribution and number of calls                                       ----
#-----------------------------------------------------------------------------------------------------#

plot_coverages <- function(scores_file, cutoff) {
  
  # Load scores file
  file = readRDS(scores_file)

  # Log transform
  # file_trans <- round(log(file, base = 10), digits = 1)
  # summary(file_trans)
  # gc()
  
  # Subset data frame
  block_size <- 10 # creates subsets of 10
  index_list <- split(seq_len(ncol(file)), ceiling(seq_len(ncol(file)) / block_size))
  
  # Plot boxplots
  destination <- file.path(s_ROOT_dir, s_out_folder, paste0("Plots/QC6_Coverages_boxplot_", cutoff,"X.pdf"))
  pdf(destination, width = 11.7, height = 8.3)
  
  progress_bar <- txtProgressBar(min=0, max=length(index_list)) # Add progress bar
  
  for (i in 1:length(index_list)){
    index <- index_list[[i]]
    boxplot(file[,index],
            ylab = "Coverage",
            xlab = colnames(file)[index],
            ylim = c(0, 30),
            las = 2,
            col="steelblue", frame = TRUE)
    
    setTxtProgressBar(progress_bar, value = i)
  }
  dev.off()
  gc()
}

plot_coverages(file.path(s_OUT_dir, "Bedmethyl_all/", paste0("bedfiles_scores_all_", cutoff, "X.rds")), cutoff=cutoff)


#-----------------------------------------------------------------------------------------------------#
#			  	6b. Coverage distribution and number of calls                                       ----
#-----------------------------------------------------------------------------------------------------#

# Function to calculate the total number of (high confidence) calls per sample
calc_number_calls <- function(scores_file, cutoff) {
  
  # Calculate total number of calls
  file = readRDS(scores_file)
  
  number_calls <- t(data.frame("total_calls" = colSums(file != 0)/1000000, 
                               "hc_calls" = apply(file, 2, function(x) sum(x >= cutoff, na.rm = TRUE)/1000000)))
  
  gc
  
  number_calls_sorted <- number_calls[,order(number_calls[2,], decreasing = TRUE)] 
  
  # Plot number of calls 
  destination <- file.path(s_ROOT_dir, s_out_folder, paste0("Plots/QC6_NumberCalls.svg"))
  svglite(destination, width = 30, height = 10)
  
  graphics::barplot(number_calls_sorted, names.arg=colnames(number_calls_sorted),
          ylab = "Number of total calls (x 1,000,000)",
          xlab="",
          main=paste0("number of calls with >=", cutoff, "X"),
          horiz = FALSE, las=2, legend = TRUE, beside=TRUE, col = c("azure2", "slategray3"))
  
  dev.off()
  
  # Save number of calls
  write.csv(as.data.frame(number_calls), file.path(s_ROOT_dir, s_out_folder, paste0("QC/number_calls_",cutoff,"X.csv")))
  return(number_calls)
}



# Store number of calls
number_calls <- calc_number_calls(scores_file = file.path(s_OUT_dir, "Bedmethyl_all/bedfiles_scores_all.rds"), cutoff = cutoff)

# Add number of calls to QCmetrics
QCmetrics <- cbind(QCmetrics, "total_calls" =  number_calls[1,], "hc_calls" = number_calls[2,])

# free up memory 
rm(number_calls)
gc()


#-----------------------------------------------------------------------------------------------------#
#					 7. Basic stats and Batch effects             ----
#-----------------------------------------------------------------------------------------------------#

# Create list of bedmethyl files
temp = list.files(file.path(s_OUT_dir, "Bedmethyl_all/"), pattern=paste0(cutoff, "X.rds"), 
                  full.names = TRUE)

# Calculate basic stats
stats_all <- pblapply(temp, function(bedfile) {
  
  file = readRDS(bedfile)
  print(dim(file))
  
  mod_type <- sub("bedfiles_", "", sub(paste0("_all_", cutoff, "X.rds*"), "", basename(bedfile)))
  file_base_name <- sub("bedfiles_", "", sub(".rds*", "", basename(bedfile)))
  
  # Calculate median and standard deviation per sample
  median_mod <- colMedians(as.matrix(file), na.rm = TRUE)
  sd_mod <- colSds(as.matrix(file), na.rm = TRUE)
  
  stats <- data.frame(median_mod, sd_mod)
  colnames(stats) <- c(paste0("median_", mod_type), paste0("sd_", mod_type))
  
  # Calculate median per batch (flowcell)
  batch_median <- aggregate(median_mod, by = list(pheno_subset$Flowcell), FUN = median)
  
  # Plot each flowcell / Group as a boxplot
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC7_BoxplotFlowcell_", file_base_name, ".svg", sep="")
  svglite(file=destination)
  
  nCol<-length(unique(pheno_subset$Flowcell)) ## assumes there is a column called Flowcell in your phenotype file
  boxplot(median_mod ~ pheno_subset$Flowcell, ylab = paste0("Median ",file_base_name, " (%)"), 
          xlab = "", las = 2, horizontal = FALSE ,col = viridis(nCol)) 
  
  dev.off()
  
  
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC7_BoxplotGroup_", file_base_name, ".svg", sep="")
  svglite(file=destination)
  
  nCol <- length(unique(pheno_subset$Group))
  boxplot(median_mod ~ pheno_subset$Group, ylab = paste0("Median ",file_base_name, " (%)"), 
          xlab = "Group", las = 2, horizontal = FALSE, col = brewer.pal(n = nCol, name = "Blues")) 
  
  dev.off()
  
  # Free memory
  rm(file, median_mod, sd_mod, batch_median)
  gc()

  return(stats)
})

# Combine list of data frames into one data frame
stats_all <- bind_cols(stats_all)

rownames(stats_all) <- pheno_subset$sample_names

# Save stats_all
saveRDS(stats_all, file.path(s_ROOT_dir, s_out_folder, paste0("QC/stats_all_", cutoff, "X.rds")))
write.csv(as.data.frame(stats_all), file.path(s_ROOT_dir, s_out_folder, paste0("QC/stats_all_", cutoff, "X.csv")))

QCmetrics <- cbind(QCmetrics, stats_all)

# Alternatively color points in original scatterplot by flowcell 
# (relevant if more variation between samples)
destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC7_ScatterplotFlowcell_", cutoff, "X.svg", sep="")
svglite(file=destination)

nCol<-length(unique(pheno_subset$Flowcell))
plot(stats_all$median_methyl, stats_all$median_hydroxymethyl, pch = 16, xlab = "Median 5mC (%)", ylab = "Median 5hmC (%)", col = viridis(nCol)[factor(pheno_subset$Flowcell)])
legend("topright", levels(factor(pheno_subset$Flowcell)), col = viridis(nCol), pch = 16, cex = 0.4)

dev.off()

# Scatterplot - relation to scores

destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC7_ScatterplotScores_hydroxymethyl_", cutoff, "X.svg", sep="")
svglite(file=destination)

ggscatter(stats_all, x="median_scores", y="median_hydroxymethyl", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "#505050", fill = "darkgrey")) + 
  geom_point(colour="steelblue2", size= 3) + 
  theme_bw() + 
  stat_cor(method = "pearson", label.x = 15, label.y = 10)

dev.off()

destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC7_ScatterplotScores_methyl_", cutoff, "X.svg", sep="")
svglite(file=destination)

ggscatter(stats_all, x="median_scores", y="median_methyl", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "#505050", fill = "darkgrey")) + 
  geom_point(colour="steelblue2", size= 3) + 
  theme_bw() + 
  stat_cor(method = "pearson", label.x = 15, label.y = 10)

dev.off()

destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC7_ScatterplotScores_modified_", cutoff, "X.svg", sep="")
svglite(file=destination)

ggscatter(stats_all, x="median_scores", y="median_modified", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "#505050", fill = "darkgrey")) + 
  geom_point(colour="steelblue2", size= 3) + 
  theme_bw() + 
  stat_cor(method = "pearson", label.x = 15, label.y = 10)

dev.off()


#-----------------------------------------------------------------------------------------------------#
#        8. Filter sample outliers                                                 -----
#-----------------------------------------------------------------------------------------------------#

stats_all <- readRDS(file.path(s_ROOT_dir, s_out_folder, paste0("QC/stats_all_", cutoff, "X.rds")))

detect_outliers <- function(stats_df, cutoff) {
 
  # plot boxplots
  stats_df_median <- stats_df %>% dplyr::select(contains("median"))

  
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC8_OverviewOutliers_boxplot_", cutoff, "X.svg")
  svglite(file=destination)
  
  boxplot(stats_df_median, las = 2, main = "With outliers")
  
  dev.off()
  
  # Plot histograms
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC8_OverviewOutliers_hist_", cutoff, "X.svg")
  svglite(file=destination)

  par(mfrow=c(3,3))

  for (col in colnames(stats_df_median)) {
    hist(stats_df_median[,col], main = "", xlab=col)
  }

  dev.off()
  
  
  # Determine outliers in 5-hmC and 5-mC
  out_hm <- identify_outliers(stats_df, variable = "median_hydroxymethyl")
  out_m <- identify_outliers(stats_df, variable = "median_methyl")
  
  # Overlapping samples
  out_modification <- intersect(rownames(out_hm), rownames(out_m)) # in both 5mC and 5hmC 
  
  # # Determine outliers in coverage
  # out_coverage <- identify_outliers(stats_df, variable = "median_scores")

  # Determine low coverage samples (based on mean)
  # low_cov_samples <- rownames(stats_df)[stats_df$median_scores < 5]
  
  low_cov_samples <- pheno_subset$sample_names[pheno_subset$Coverage < 5]
  
  # no_cov_samples <- rownames(stats_df)[stats_df$median_scores < 1]

  # outliers <- unique(c(out_modification, rownames(out_coverage))) # in either 5mC or 5hmC
  outliers <- intersect(low_cov_samples, out_modification)
  # outliers <- unique(c(outliers, no_cov_samples))
  
  # plot post outlier
  
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC8_OverviewOutliers_boxplot_post_", cutoff, "X.svg")
  svglite(file=destination)
  
  boxplot(stats_df_median[!rownames(stats_df_median) %in% outliers,], las = 2, main = "Without outliers")
  
  dev.off()
  
  
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC8_OverviewOutliers_hist_post_", cutoff, "X.svg")
  svglite(file=destination)
  
  par(mfrow=c(3,3))
  
  for (col in colnames(stats_df_median)) {
    hist(stats_df_median[!rownames(stats_df_median) %in% outliers,col], main = "", xlab=col)
  }
  
  
  # hist(stats_df_median[!rownames(stats_df_median) %in% outliers,])
  
  dev.off()
  
  return(outliers)
}

outliers <- detect_outliers(stats_all, cutoff = cutoff)

outliers

# Mark samples that are outliers
outliers_filter <- rownames(stats_all) %in% outliers

# Samples that do not pass threshold
QCmetrics<-cbind(QCmetrics,outliers_filter)
SamplesFail[which(outliers_filter==TRUE)]<-TRUE
Step1<-c(length(which(outliers_filter==TRUE)),sum(SamplesFail))
Stepsummary<-cbind(Stepsummary,Step1)
print(Stepsummary)

# free up memory
rm(stats_all)
gc()


#-----------------------------------------------------------------------------------------------------#
#        9. Check sex                                                            ----
#-----------------------------------------------------------------------------------------------------#

# # Filter for sex chromosomes
# filter_sexchr <- function(scores_file) {
#   positions <- data.frame(rownames(readRDS(scores_file)))
#   colnames(positions) <- "position"
#   positions_sex <- positions %>% filter(str_detect(position, "chrX_")) # only chr X, works better
#   # positions_sex <- positions %>% filter(str_detect(position, "chrX_|chrY_")) # both chrX and chrY
#   return(positions_sex$position)
# }
#   
# positions_sexchr <- filter_sexchr(scores_file = file.path(paste0(s_OUT_dir, "Bedmethyl_all/bedfiles_scores_all.rds")))
# 
# 
# # Load function findGenderPC
# findGenderPC <- function(modifications, sex, npcs = 4, file_base_name){
#   
#   # Select no outlier cases 
#   modifications <- modifications[,!colnames(modifications) %in% outliers]
#   
#   # Select complete cases
#   modifications_all_com<-modifications[complete.cases(modifications),]
#   
#   # with 5X cutoff only leaves 1786 positions
#   
#   # Select positions on the sex chromosomes
#   modifications_all_com <- modifications_all_com[rownames(modifications_all_com) %in% positions_sexchr,]
#   pca<-prcomp(modifications_all_com)
#   
#   # Correlate prinicpal components with factor Gender
#   pca.cor <- rep(NA, npcs)
#   
#   for(i in 1:npcs){
#     pca.cor[i]<-cor(pca$rotation[,i], as.numeric(as.factor(sex)), use = "complete")  
#     }
#   
#   top<-order(abs(pca.cor), decreasing = TRUE)[1]
#   second<-order(abs(pca.cor), decreasing = TRUE)[2]
#   print(paste("Top correlated principal components with sex:", top, ",", second))
#   
#   # Plot the top principal components, color by sex 
#   destination <- file.path(s_ROOT_dir, s_out_folder, paste0("Plots/QC9_PredSex_", file_base_name, ".svg"))
#   svglite(destination)
#   
#   plot(pca$rotation[,top], pca$rotation[,second], pch = 16, col = c("green", "darkblue")[as.factor(sex)],
#     xlab = paste("PC", top), ylab = paste("PC", second))
#   legend("topright", levels(as.factor(sex)), pch = 16, col = c("green", "darkblue"))
#   
#   dev.off()
#   
#   # Predict sex
#   predSex <- rep(NA, length(sex))
#   options.sex<-levels(as.factor(sex))
#   
#   if(abs(pca.cor[top]) > 0.9){ print("Top PC has r > 0.9 with sex so good enough to confirm reported sexes") } else
#   {print(paste("Top PC has r =", round(abs(pca.cor[top]),2), "with sex so may not be good enough to confirm reported sexes")) }
#   
#   if(sign(pca.cor[top]) == 1){
#     predSex[which(pca$rotation[,top] < 0)] <- options.sex[1]
#     predSex[which(pca$rotation[,top] > 0)] <- options.sex[2]} 
#   
#   else {
#    predSex[which(pca$rotation[,top] < 0)] <-options.sex[2]
#    predSex[which(pca$rotation[,top] > 0)] <-options.sex[1]} 
#   
#   predSex <- data.frame(predSex) 
#   colnames (predSex) <- paste0("predSex_", file_base_name)
#   gc()
#   return(predSex)
# }
# 
# # Use findGenderPC function to predict sex 
# # ReportedSex <- pheno_subset$Gender # all samples
# ReportedSex <- pheno_subset[which(!pheno_subset$sample_names %in% outliers),"Gender"] # without outliers
# 
# temp = list.files(file.path(s_OUT_dir, "Bedmethyl_all/"), pattern=paste0("all_",cutoff, "X.rds"), 
#                   full.names = TRUE)
# 
# temp = str_subset(temp, pattern = "count", negate = TRUE)
# 
# PredictedSex <- pblapply(temp, function(bedfile) {
#   # Read RDS file
#   file = readRDS(bedfile)
#   file_base_name <- sub("bedfiles_", "", sub(".rds*", "", basename(bedfile)))
#   
#   # Predicted sex
#   predSex <- findGenderPC(file,ReportedSex, npc=4, file_base_name) # npcs = 16
#   return(predSex)
# })
# 
# # Add predicted sex to QCmetrics
# PredictedSex <- bind_cols(PredictedSex)
# # QCmetrics <- left_join(QCmetrics, PredictedSex)
# ############## 
# 
# # Update stepssummary
# # SamplesFail[which(PredictedSex!=ReportedSex)] <- TRUE 
# # Step3<-c(length(which(PredictedSex!=ReportedSex)),sum(SamplesFail)) 
# # Stepsummary<-cbind(Stepsummary,Step3)
# # print(Stepsummary)
# 
# # Free up memory
# rm(PredictedSex, positions_sexchr)
# gc()

# Error in svd(x, nu = 0, nv = k) : a dimension is zero
# No positions left in sex chromosomes - too many samples with NA values

#-----------------------------------------------------------------------------------------------------#
#                         10. Check SNPs                ----
#-----------------------------------------------------------------------------------------------------#
# 
# # Create list of VCF files
# temp = list.files(s_vcfs_folder, pattern="\\.vcf$", full.names = TRUE)
# 
# pblapply(temp, function(x) {
# 
#   # Load VCF file
#   vcf <- read.vcfR(x)
#   
#   # Extract DP coverage values
#   dp.matrix <- vcfR::extract.info(vcf, element = "DP", 
#                                   as.numeric = TRUE)
#   
#   # Extract values for which the depth coverage > 40
#   vcf@fix <- vcf@fix[dp.matrix > 45,]
#   vcf@gt <- vcf@gt[dp.matrix > 45,]
# 
#   # Extract alternative nucleotide
#   vcf_snps_num <- t(vcfR2loci(vcf, return.alleles=TRUE))
# 
#   vcf_snps_num <- data.frame("position" = rownames(vcf_snps_num), "alleles" = vcf_snps_num)
# 
#   # Extract the base name of the file
#   file_base_name <- basename(x)
# 
#   # Save processed data
#   saveRDS(vcf_snps_num, file = file.path(s_OUT_dir, "VCF/", paste0(file_base_name, "_alleles.rds")))
# 
#   # Free up memory
#   rm(vcf, vcf_snps_num)
#   gc()
# 
# }, cl = 4)
# 
# 
# # Combine the VCF files
# temp = list.files(file.path(s_OUT_dir, "VCF/"),
#                     pattern=paste0("\\.rds$"),
#                     full.names = TRUE)
# 
# vcf_all <-readRDS(temp[1]) # Start with first sample
# 
# progress_bar <- txtProgressBar(min=0, max=length(temp)) # Add progress bar
# 
# for (i in 2:length(temp)) {
#   # Read BED file
#   file <- data.frame(readRDS(temp[i]))
#   file_base_name <- sub(".bam_igv.vcf_scores*", "", basename(temp[i]))
# 
#   # Combine methyl values
#   vcf_all <- merge(vcf_all, file, by = "position", all=TRUE)
#   colnames(vcf_all)[1] <- c("position")
#   colnames(vcf_all)[i] <- sub(".bam_igv.vcf_alleles.rds*", "", basename(temp[i-1])) # existing data frame
#   colnames(vcf_all)[i+1] <- sub(".bam_igv.vcf_alleles.rds*", "", basename(temp[i])) # new data frame
# 
#   # Free memory
#   rm(file)
#   gc()
# 
#   setTxtProgressBar(progress_bar, value = i)
# }
# 
# # Move position values to rownames
# rownames(vcf_all) <- vcf_all$position
# vcf_all <- vcf_all %>% dplyr::select(!(position))
# 
# # Set sample names
# # names(vcf_all) <- paste0("sample_", pheno_subset$CaseID)
# # colnames(vcf_all_factor) <- paste0("sample_", seq.int(names(vcf_all)))
# 
# # select complete cases
# vcf_all_complete <- vcf_all[complete.cases(vcf_all),]
# head(vcf_all_complete)
# 
# # Factorize the rows based on genotypes (e.g. 0, 1,2)
# vcf_all_factor <- t(apply(vcf_all_complete, 1, function(row) as.numeric(factor(row))))
# 
# # save
# saveRDS(vcf_all_complete, file = file.path(s_OUT_dir, "VCF/", paste0("vcf_all_complete.rds")))
# saveRDS(vcf_all_factor, file = file.path(s_OUT_dir, "VCF/", paste0("vcf_all_factor.rds")))
# 

#-----------------------------------------------------------------------------------------------------#
#      11. Check genetically identical samples correlate across SNP probes        ----
#-----------------------------------------------------------------------------------------------------#

# # Load VCF dataframes
# vcf_all_complete <- readRDS(file = file.path(s_OUT_dir, "VCF/", paste0("vcf_all_complete.rds")))
# vcf_all_factor <- readRDS(file = file.path(s_OUT_dir, "VCF/", paste0("vcf_all_factor.rds")))
# 
# # ~ Genetic correlations -------------------------------------------
# snpCor<-cor(vcf_all_factor, use="pairwise.complete.obs")
# round(snpCor,3)
# 
# rcorr(vcf_all_factor, type = c("pearson","spearman")) # with p-values
# 
# #  ~  Create plots for correlation ---------------------------------
# 
# # Plot dendrogram
# hclust_res <- hclust(dist(t(vcf_all_factor)), method = "complete")
# plot(hclust_res)
# 
# # Visualize correlation matrix (with numbers)
# corrplot(snpCor, type = "upper", order = "hclust",
#          tl.col = "black", tl.srt = 45, method = 'number', addrect = 3,  col = COL2('RdBu', 20), col.lim=c(-1,1))
# 
# # Visualize correlation matrix (with circles)
# corrplot(snpCor, type = "upper", order = "hclust",
#          tl.col = "black", tl.srt = 45, addrect = 3,  col = COL2('RdBu', 20), col.lim=c(0,1))
# 
# # Visualize correlation matrix in heatmap
# col<- colorRampPalette(c("blue", "white", "red"))(20)
# heatmap(x = snpCor, col = col, symm = TRUE)
# 
# heatmap_vcfs <- function(factorized_alleles) {
#   # Get overview of rows with identical values across columns
#   rows_no_diff <- apply(vcf_all_factor, 1, function(row){
#     all(row == row[1])
#   })
# 
#   print(summary(rows_no_diff))
# 
#   # Get vcf_genotypes with differences
#   vcf_dif <- vcf_all_factor[rows_no_diff == FALSE,]
# 
#   # Make heatmap
#   heatmap(vcf_dif)
# }
# 
# heatmap_vcfs(vcf_all_factor)
# 
# 
# # ~ Calculate Max correlation -----------------------------
# # Ignore correlations between a sample and itself
# for(i in 1:ncol(vcf_all_factor)){
#   snpCor[i,i]<-NA
# }
# 
# corMax<-apply(snpCor, 1, max, na.rm = TRUE)
# hist(corMax, xlab = "Max. correlation with all other samples", main = "")
# corMax
# 
# # Add corMax to QCmetrics
# QCmetrics<-cbind(QCmetrics, corMax)
# SamplesFail[which(corMax>0.8)]<-TRUE
# Step5<-c(sum(corMax>0.8),sum(SamplesFail))
# Stepsummary<-cbind(Stepsummary,Step5)
# print(Stepsummary)


#-----------------------------------------------------------------------------------------------------#
#          12. Save QC Results of the failed and passed samples                 -----
#-----------------------------------------------------------------------------------------------------#

write.csv(QCmetrics[SamplesFail,], paste0(s_ROOT_dir, s_out_folder, "QC/SamplesFailedQC_",cutoff,"X.csv"), row.names = FALSE)
write.csv(QCmetrics[!SamplesFail,], paste0(s_ROOT_dir, s_out_folder, "QC/SamplesPassedQC_",cutoff,"X.csv"), row.names = FALSE)


#-----------------------------------------------------------------------------------------------------#
#          13. Filter passed data                                             -----
#-----------------------------------------------------------------------------------------------------#

# Passed samples pheno data
pheno_subset_pass <- pheno_subset[!SamplesFail,]

samples_pass <- pheno_subset_pass$sample_names

# Passed samples 5 mC and 5-hmc data
temp = list.files(file.path(s_OUT_dir, "Bedmethyl_all/"), pattern=paste0(cutoff, "X"), 
                  full.names = TRUE)

temp = str_subset(temp, "pass", negate = TRUE)


pblapply(temp, function(bedfile) {
  
  # Read file
  file <- readRDS(bedfile)
  cutoff = cutoff
  
  # Extract CpGs
  file_filtered <- file[,colnames(file) %in% samples_pass]
  
  # Extract the base name of the file
  file_base_name <- sub("bedfiles_", "", sub(".rds*", "", basename(bedfile)))
  
  # Save processed data
  saveRDS(file_filtered, file = file.path(s_OUT_dir, "Bedmethyl_all/", paste0("bedfiles_", file_base_name, "_pass.rds")))
  
  # Free up memory
  rm(file)
  gc()
  
}, cl = 4)


#-----------------------------------------------------------------------------------------------------#
#          14. Phenotype data                                  ----
#-----------------------------------------------------------------------------------------------------#

# Gender
gender_plot <- ggplot(pheno_subset_pass, aes(x = Gender, fill = Group)) +
  geom_bar(position = "dodge", width=0.8) + 
  scale_fill_brewer(palette="Blues", direction = -1) +  
  theme_bw()

# ApoE
ApoE_plot <- ggplot(pheno_subset_pass, aes(x = ApoE, fill = Group)) +
  geom_bar(position = "dodge", width=0.8) +
  scale_fill_brewer(palette="Blues", direction = -1) +  theme_bw()

# Braak score
braak_plot <- ggplot(pheno_subset_pass, aes(x = Braak_score, fill = Group)) +
  geom_bar(position = "dodge", width=0.8) +
  scale_fill_brewer(palette="Blues", direction = -1) +  theme_bw() +
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6"))


# MMSE scores
mmse_plot <- ggplot(pheno_subset_pass, aes(x = MMSE_score, fill = Group)) +
  geom_histogram(position = "identity", alpha = 0.9, bins = 10) +
  facet_wrap(~ Group) +
  scale_fill_brewer(palette="Blues", direction = -1) +
  theme_bw()

# Age
age_plot <- ggplot(pheno_subset_pass, aes(x = Age, fill = Group)) +
  geom_histogram(position = "identity", alpha = 1.0, bins = 20) +
  facet_wrap(~ Group) +
  scale_fill_brewer(palette="Blues", direction = -1) +
  theme_bw()


# DNA input
dna_plot <- ggplot(pheno_subset_pass, aes(x = DNA_input, fill = Group)) +
  geom_histogram(position = "identity", alpha = 1.0, bins = 10) +
  facet_wrap(~ Group) +
  scale_fill_brewer(palette="Blues", direction = -1) +
  theme_bw()


destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC14_phenotypeoverview_",cutoff,"X.svg", sep="")

svglite(destination, width = 15, height = 10)

ggarrange(gender_plot, ApoE_plot, braak_plot, mmse_plot, age_plot, dna_plot,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow=2)

dev.off()

rm(gender_plot, ApoE_plot, braak_plot, mmse_plot, age_plot, dna_plot)

gc()

#-----------------------------------------------------------------------------------------------------#
#          15. Create matrix of weights                                         ----
#-----------------------------------------------------------------------------------------------------#

# ~ weights based on all samples  --------------------------

weights_single_sample <- function(scores_file) {
  
  scores_filtered <- readRDS(scores_file)

  # Compute weights using vectorized row-wise division
  # scores_weights <- sweep(scores_filtered, 1, rowSums(scores_filtered, na.rm = TRUE), `/`)
  scores_weights <- sweep(scores_filtered, 1, rowSums(scores_filtered), `/`)
  
  print(head(scores_weights))
  saveRDS(scores_weights, file = paste0(s_ROOT_dir,s_out_folder,"QC/scores_global_weights_", cutoff,"X_pass.rds")) 
  
  # Free up memory
  rm(scores_filtered, scores_weights)
  gc()
}

weights_single_sample(scores_file = file.path(paste0(s_OUT_dir, "Bedmethyl_all/bedfiles_scores_all_", cutoff, "X_pass.rds")))


# ~ weights based on all samples within group -----------------------------

weights_global_sample <- function(scores_file, groups) {
  
  # Prepare list of groups
  group_weights_list <- list()
  unique_groups <- groups 
  
  # Read scores file
  scores_filtered <- readRDS(scores_file)
  
  for (group in unique_groups) {
    # samples that belong to a specific group
    group_samples <- pheno_subset_pass$sample_names[pheno_subset_pass$Group == group]
    
    # scores for these specific samples
    group_scores <- scores_filtered[,colnames(scores_filtered) %in% group_samples]
    
    # scores divided by total sum of scores
    group_weights <- round(sweep(group_scores, 1, rowSums(scores_filtered, na.rm=TRUE), `/`), digits = 3)
    
    # return weights
    group_weights_list[[group]] <- group_weights
  }
  
  # Combine group weights
  scores_weights <- do.call(cbind, group_weights_list)
  
  # remove prefix (AD, MCI, CTL)
  colnames(scores_weights) <- unlist(lapply(group_weights_list, function(x) colnames(x)))
  
  # Order weights based on Case ID
  scores_weights <- scores_weights[,colnames(scores_filtered)]
  print(head(scores_weights))
  
  # Save file
  saveRDS(scores_weights, file = paste0(s_ROOT_dir,s_out_folder,"QC/scores_group_weights_", cutoff,"X_pass.rds")) 
  
  # Free up memory
  rm(scores_weights,group_scores, group_weights, scores_filtered)
  gc()
}

weights_global_sample(scores_file = file.path(paste0(s_OUT_dir, "Bedmethyl_all/bedfiles_scores_all_", cutoff, "X_pass.rds")), 
                      groups = unique(pheno_subset_pass$Group))


#-----------------------------------------------------------------------------------------------------#
#					        16. Density plots per group                ----
#-----------------------------------------------------------------------------------------------------#

# ! Very heavy ! 
# temp = list.files(file.path(s_OUT_dir, "Bedmethyl_all/"), pattern=paste0(cutoff, "X"), 
#                   full.names = TRUE)
# 
# pblapply(temp, function(bedfile) {
#   
#   # Read RDS file
#   file = readRDS(bedfile)
# 
#   mod_type <- sub("bedfiles_", "", sub(paste0("_all_", cutoff, "X.rds*"), "", basename(bedfile)))
#   file_base_name <- sub("bedfiles_", "", sub(".rds*", "", basename(bedfile)))
#   
#   
#   # Data per group
#   dat_mod <- reshape2::melt(file)
#   colnames(dat_mod) <- c("sample_names", "value")
#   dat_mod <- dat_mod %>% left_join(pheno_subset[,c("sample_names","Group")], by = join_by(sample_names))
#   gc()
#   
#   # Basic density plots, separate graphs per group
#   images <- ggplot(dat_mod, aes(x=value, fill = Group)) + 
#     geom_density(alpha = 0.8) + scale_fill_manual(values=c("#505567", "#E5E6EB", "blue4"))+ 
#     facet_grid(Group ~ .) + xlab("Beta value") +
#     ylab("Density") +
#     theme_classic() +
#     theme(
#       axis.title.x = element_text(size = 15, face = "bold"),
#       axis.title.y = element_text(size = 15, face = "bold"),
#       axis.text.x = element_text(size = 15),
#       axis.text.y = element_text(size = 15),
#       legend.text = element_text(size = 15),
#       legend.title = element_text(size = 15, face ="bold"))
#   
#   ggexport(images,  width = 2000, height = 2500, res = 200,
#            filename = paste0(s_ROOT_dir, s_out_folder, "Plots/BetasDensities_separate_", file_base_name,".png"))
#   
#   # Basic density plot; groups overlapping
#   images <- ggplot(dat_mod, aes(x=value, fill = Group)) + 
#     geom_density(alpha = 0.8) + 
#     scale_fill_manual(values=c("#505567", "#E5E6EB", "blue4"), labels = c("CTL", "MCI", "AD")) +
#     xlab("Beta value") +
#     ylab("Density") + theme_classic() +
#     theme(
#       axis.title.x = element_text(size = 15, face = "bold"),
#       axis.title.y = element_text(size = 15, face = "bold"),
#       axis.text.x = element_text(size = 12),
#       axis.text.y = element_text(size = 12),
#       legend.text = element_text(size = 15),
#       legend.title = element_text(size = 15, face ="bold")
#     ) 
#   
#   ggexport(images, width = 3000, height = 1500, res = 200,
#            filename = paste0(s_ROOT_dir, s_out_folder, "Plots/BetasDensities", file_base_name, ".png"))
#   
#   # Free up memory
#   rm(images, dat_mod, file)
#   gc()
# })


#-----------------------------------------------------------------------------------------------------#
#					      17. Density plots per modification type           ------
#-----------------------------------------------------------------------------------------------------#

#### MAKE IT INTERACTIVE?

# List of samples 
temp = list.files(file.path(s_project_folder, "Bedmethyl/"), pattern="bam.bed.rds", 
                  full.names = TRUE)

file_base_name <- sub(".bam.bed.rds*", "", basename(temp))

# ~ Density plots per sample -------------------
destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC17_ModificationDensities_prefilter.pdf", sep="")
pdf(destination)

pblapply(file_base_name, function(file) {
  
  # Retrieve dataframes
  modified <- readRDS(file.path(s_project_folder, "Bedmethyl/", paste0(file, ".bam.bed_modified.rds")))
  methyl <- readRDS(file.path(s_project_folder, "Bedmethyl/", paste0(file, ".bam.bed_methyl.rds")))
  hydroxymethyl <- readRDS(file.path(s_project_folder, "Bedmethyl/", paste0(file, ".bam.bed_hydroxymethyl.rds")))
  
  # Plot the densities 
  plot(density(modified[,2], na.rm=TRUE), main = c("Density C modifications (%)", file, sep = ""), col = "darkgrey") 
  lines(density(methyl[,2], na.rm=TRUE), col = "darkblue")
  lines(density(hydroxymethyl[,2], na.rm=TRUE), col = "orange")
  legend("topright", legend=c("5mC+5hmC", "5mC", "5hmC"), lty=1, col=c("darkgrey","darkblue", "orange"))
  
  # Clear memory
  rm(modified, methyl, hydroxymethyl)
  gc()
  
})

dev.off()


# ~ Density plots on all data pre and post outlier filtering --------------------------------------

plot_densities <- function(mod_type, file_suffix) {
  
  # Get colors based on group
  colors <- as.character(factor(pheno_subset_pass$Group, labels = c("#A8CDECFF","#F6955EFF", "#682C37FF")))
  
  # Retrieve data
  file <- readRDS(file.path(paste0(s_OUT_dir, "Bedmethyl_all/", "bedfiles_", mod_type, "_all_", file_suffix, ".rds")))
  
  # Calculate max density 
  max_density <- max(apply(file, 2, function(x) {
    max(density(x, na.rm=TRUE)$y)
    }))
  
  # Plot densities
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/QC17_ModificationDensities_", mod_type, "_", file_suffix, ".svg")
  svglite(destination)
  
  plot(density(as.matrix(file)[,1], na.rm=TRUE), 
       main = c("Density C modifications (%)", mod_type, sep = ""), col = colors[1],
       ylim=c(0, max_density)) 
  
  for (i in 2:(length(file))) { lines(density(file[,i], na.rm=TRUE), col = colors[i])}
  
  legend("topright", legend=c("CTL", "MCI", "AD"), lty=1, col=c("#A8CDECFF","#F6955EFF", "#682C37FF"))
  
  dev.off()
  
  # Free up memory
  rm(file, max_density)
  gc()
}


# Prefilter outliers

plot_densities("methyl", file_suffix = paste0(cutoff, "X"))
plot_densities("hydroxymethyl", file_suffix = paste0(cutoff, "X"))
plot_densities("modified", file_suffix = paste0(cutoff, "X"))

# Post filtering outliers

plot_densities("methyl", file_suffix = paste0(cutoff, "X_pass"))
plot_densities("hydroxymethyl", file_suffix = paste0(cutoff, "X_pass"))
plot_densities("modified", file_suffix = paste0(cutoff, "X_pass"))


#-----------------------------------------------------------------------------------------------------#
#                       18. Check PCA correlations (first 10 comps)         ----
#-----------------------------------------------------------------------------------------------------#

pca_correlations <- function(mod, filesuffix) {
  
  ### Filter pheno file (outliers, complete cases)
  pheno_subset_out <- pheno_subset_pass[complete.cases(pheno_subset_pass),]
  pheno_subset_out <- pheno_subset_out[, !sapply(pheno_subset_out, function(x) all(duplicated(x)[-1L]))]
  
  # Remove Barcode factor 
  pheno_subset_out <- dplyr::select(pheno_subset_out, -c("Barcode"))
  
  # Select no outlier cases 
  mod <- mod[,colnames(mod) %in% pheno_subset_out$sample_names]
  
  # Convert to matrix 
  mod <- as.matrix(mod)
  
  # Recalculating PCA as the betas were updated 
  pca <- prcomp(t(mod[complete.cases(mod),]))
  
  # Plot some 
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/", "QC18_ScreePlot_", filesuffix, ".pdf", sep="")
  pdf(file=destination) 
  
  plot(pca) # scree plot, indicating which PC has most variance. 
  dev.off()
  
  Rel_variance_expl = round(pca$sdev^2 / sum(pca$sdev^2),3) # relative variance 
  cat("First component explains",Rel_variance_expl[1]*100, "% of variance in the data")
  
  # Transform the pheno data to a numeric matrix, including character columns. Should usually take a lot of attention, as some columns do 
  # not make sense to undergo this transformation. This is a lazy workaround. For this example its OK.
  pheno_numeric = apply(pheno_subset_out,2,function(x){as.numeric(as.factor(x))}) 
  pheno_numeric = pheno_numeric[,which(apply(pheno_numeric,2,function(x){sum(is.na(x))==0}))]
  
  # Correlate the PCs (max 10) to the phenodata, to see what effects can be explained by variance alone. Ideally, PC 1 is 
  # the intended experimental effect.
  correlation_frame = as.data.frame(na.omit(round(t(data.frame(cor(pca$x[,1:min(10,dim(mod)[2])],pheno_numeric))),2))) 
  correlation_frame[abs(correlation_frame)<0.3] = 0
  
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/", "QC18_PCA_", filesuffix, ".pdf", sep="")
  pdf(file=destination) 
  
  report_pc = 1:10 
  for(PC in report_pc){ # State what is correlating strongly to what 
    cat("\n",rownames(correlation_frame)[which.max(abs(correlation_frame[,PC]))],"(",which.max(abs(correlation_frame[,PC])),") is correlated 
        the most with PC",PC,"(",correlation_frame[which.max(abs(correlation_frame[,PC])),PC],")\n\n")
    
    # make colors based on factors
    lazycolors = as.numeric(as.factor(pheno_subset_out[,which.max(abs(correlation_frame[,PC]))]))
    
    # plot title
    main_title = paste0("PC ",PC," - ",rownames(correlation_frame)[which.max(abs(correlation_frame[,PC]))],"(",correlation_frame[which.max(abs(correlation_frame[,PC])),PC],")" )
    
    # plot pc score and color by correlate
    if(is.numeric(pheno_subset_out[,which.max(abs(correlation_frame[,PC]))])){
      plot(pca$x[,PC]~pheno_subset_out[,which.max(abs(correlation_frame[,PC]))],
           col=lazycolors,
           pch=19,
           xlab = rownames(correlation_frame)[which.max(abs(correlation_frame[,PC]))],
           ylab = paste0("PC ",PC," (",Rel_variance_expl[PC]*100,"%)"),
           main=main_title)
    }else{
      boxplot(pca$x[,PC]~pheno_subset_out[,which.max(abs(correlation_frame[,PC]))],
              col=lazycolors,
              pch=19,
              xlab = rownames(correlation_frame)[which.max(abs(correlation_frame[,PC]))],
              ylab = paste0("PC ",PC," (",Rel_variance_expl[PC]*100,"%)"),
              main=main_title)
    }
  }
  
  # plot PC1 and PC2, colours indicating group variable
  plot(pca$x[,2] ~ pca$x[,1], col=factor(pheno_subset_out$Group),
       xlab = paste0("PC", 1, " (", Rel_variance_expl[1]*100,"%)"),
       ylab = paste0("PC", 2, " (", Rel_variance_expl[2]*100,"%)"))
  
  pca_group <- cbind(as.data.frame(pca$x), "Group" = pheno_subset_out$Group)
  
  images <- as.list(1:2)
  
  images[[1]] <- ggplot(data = pca_group, mapping = aes(x = PC1, y = PC2, color = Group)) + 
    geom_point(shape = 10) +
    scale_color_manual(values=c("#A8CDECFF","#F6955EFF", "#682C37FF"), 
                       name = "Group", 
                       labels = c("CTL", "MCI", "AD")) +
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
  pca_gender <- cbind(as.data.frame(pca$x), "Gender" = pheno_subset_out$Gender)
  
  images[[2]] <- ggplot(data = pca_gender, mapping = aes(x = PC1, y = PC2, color = Gender)) + 
    geom_point(shape = 10) +
    scale_color_manual(values=c("#A8CDECFF","#F6955EFF"), 
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
  
  gc()
  dev.off()
}

temp = list.files(file.path(s_OUT_dir, "Bedmethyl_all/"), pattern=paste0(cutoff, "X_pass"), 
                  full.names = TRUE)

pblapply(temp[6:7], function(bedfile) {
  
  # Read RDS file
  file = readRDS(bedfile)
  file_base_name <- sub("bedfiles_", "", sub(".rds*", "", basename(bedfile)))
  
  # Predicted sex
  pcaplots <- pca_correlations(file, file_base_name) # npcs = 16
  
  dev.off()
  
  gc()
  
})




#-----------------------------------------------------------------------------------------------------#
#					        	    19. M values                      ----
#-----------------------------------------------------------------------------------------------------#

# Modify modification values of zero and 100 

imputezerohundred <- function(mod) {
  # Determine the second smallest number to zero 
  min_value <- min(mod[mod > 0], na.rm = TRUE)
  
  # Determine the second largest number to 100
  max_value <- max(mod[mod < 100], na.rm = TRUE)

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

# Create list with bedfiles with samples pass
temp = list.files(file.path(s_OUT_dir, "Bedmethyl_all/"), pattern=paste0(cutoff, "X_pass"), 
                  full.names = TRUE)
temp <- temp[!str_detect(temp, "counts|scores")] # omits data about counts or scores

pblapply(temp, function(bedfile) {
  
  # Read RDS file
  file = readRDS(bedfile)
  file_base_name <- sub("bedfiles_", "", sub(".rds*", "", basename(bedfile)))
  
  # M values
  Metas = mod2M(as.data.frame(imputezerohundred(file)))
  
  # Save processed data
  saveRDS(Metas, file = file.path(s_OUT_dir, "QC/", paste0(file_base_name, "_Metas.rds")))
  
  # Histogram of M values
  destination <- file.path(paste0(s_ROOT_dir, s_out_folder, "Plots/", "QC12_M_values", file_base_name,".pdf", sep=""))
  pdf(file=destination) 
  
  hist(as.matrix(Metas))
  dev.off()
  
  # Free up memory
  rm(file, Metas)
  gc()
  
})

# which(is.infinite(Metas),arr.ind = FALSE)  # does not work



#-----------------------------------------------------------------------------------------------------#
#					        	20. Save relevant data              ----
#-----------------------------------------------------------------------------------------------------#

# Only keep the useful variables in pheno. 
pheno_new <- QCmetrics[!SamplesFail,] 

# pheno_new <- subset(pheno_new, select = -c(Intensity, M.median, U.median, PredictedSex, PredictedStrain))

saveRDS(pheno_new, file = paste0(s_ROOT_dir, s_out_folder,"QC/Pheno_new_",cutoff, "X.rds"))
write.csv(pheno_new, file = paste0(s_ROOT_dir, s_out_folder,"QC/Pheno_new_",cutoff, ".csv"), row.names = FALSE)

# Save work space
save.image(file=paste0(s_ROOT_dir,s_out_folder,"EPI-Alzheimer-ONT.RData"))

