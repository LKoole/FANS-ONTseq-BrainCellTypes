#-----------------------------------------------------------------------------------------------------#
#							Settings                                         ----
#-----------------------------------------------------------------------------------------------------#

# Get all general settings 
s_ROOT_dir <<- "path/to/root/directory"

source(paste0(s_ROOT_dir,"Scripts/.Main/Settings_v2.R"))

#-----------------------------------------------------------------------------------------------------#
#							Libraries needed                                                  ----
#-----------------------------------------------------------------------------------------------------#
# Install packages if not present
# BiocManager::install(c("limma", "DMRcate", "qqman", "gplots", "enrichplot",
# "readr", "htmltools", "stringr", "clusterProfiler", "org.Hs.eg.db", "DOSE", "bacon",
# "annotatr", "TxDb.Hsapiens.UCSC.hg38.knownGene", "purrr", "readxl", "GenomicRanges",
# "EnhancedVolcano", "GenomeInfoDb", "GenomicDistributions", "GenomeInfoDbData", "GO.db"))

library(gplots) 
library(ggplot2)
library(ggpubr) # export gplots
library(svglite) # svg plots
library(purrr)
library(readxl) # read excel
library(pbapply)
library(dplyr) 
library(stringr)
library(data.table)
library(paletteer)
library(forcats)
library(plyranges)
library(smplot2)
library(RColorBrewer)
library(patchwork)
library(tidyr)
library(scales)
library(vtable)
library(methylKit)
library(RColorBrewer)
library(RColorBrewer)
library(bsseq)
library(DSS)
library(limma) # contrasts
library(bench)
library(profvis)
library(htmlwidgets)
library(ggprism)
library(emmeans)
library(ggridges)
library(car)

#-----------------------------------------------------------------------------------------------------#
#						Load pheno data                                                           ----
#-----------------------------------------------------------------------------------------------------#
# Set cutoff coverage
cutoff = 5

# Pheno file
pheno_subset_pass <- readRDS(file = paste0(s_ROOT_dir, s_out_folder,"QC/Pheno_new_",cutoff,"X.rds"))

#-----------------------------------------------------------------------------------------------------#
#					Functions                                                       ----
#-----------------------------------------------------------------------------------------------------#

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

mod2M <- function(mod){
  return(log2((mod) /(100 - mod)))
} 


# Limma without weights function
DM_analysis_single_limma <- function(Metas_file, mod_type, s_bootstrap_folder,
                                     s_CovFormulaImportantSV,temp_design) {
  
  Metas <- readRDS(Metas_file)
  Metas <- as.matrix(Metas) 
  
  # Make LMmodel without weights
  cat(paste0("Fitting linear model"),"\n\n")
  temp_fit = limma::lmFit(as.matrix(Metas),temp_design)
  
  # colSums(is.na(temp_fit$coefficients))
  # factor.Flowcell.RESEARCH_GSM0172RRMS_03042025  ISSUE
  
  # Show the coefs for first gene
  data.frame(temp_fit$coef[1,])
  
  cat(paste0("Rank genes using empirical Bayes method"),"\n\n")
  Limmaefit = limma::eBayes(temp_fit)
  
  # Save results limma
  # tryCatch(saveRDS(Limmaefit, file = file.path(paste0(s_ROOT_dir,s_out_folder,"DM/Limmaefit_Limma_", mod_type, ".rds"))), error = function(e) warning(e$message))
  
  loop_contrast = "GroupAD"
  s_coeff = loop_contrast 
  
  cat(paste0("Extract table of the top-ranked genes", "\n\n"))
  DM_Results = limma::topTable(Limmaefit, coef=s_coeff, num=Inf, sort.by = "P", adjust="BH") # Data of specified variables
  
  # Order based on p value
  DM_Results_sign=DM_Results[DM_Results$P.Value<0.01,]
  print(head(DM_Results,20))
  
  destination <- file.path(s_bootstrap_folder, paste0("DM_hist_Limma_",mod_type, "_",loop_contrast, ".pdf"))
  pdf(file=destination)
  
  # histogram of p values and adjusted p values
  hist(DM_Results$P.Value ,xlim=c(0,1),main=s_coeff)
  hist(DM_Results$adj.P.Val ,xlim=c(0,1),main=s_coeff)
  min(DM_Results$adj.P.Val)
  
  
  # Descriptives
  cat(paste0("Minimum adjusted p value of ", s_coeff, ": " ,round(min(DM_Results$adj.P.Val),digits = 3), "\n\n"))
  
  cat(paste0("Contrast: ",s_coeff,".\n  number of features below TH (0.01)\n  Pvalue:",
             sum(DM_Results$P.Value<0.01)," (",round((sum(DM_Results$P.Value<0.01)/dim(DM_Results)[1])*100,0),"%)",
             "\n  AdjPvalue:",sum(DM_Results$adj.P.Val<0.05)," (",round((sum(DM_Results$adj.P.Val<0.05)/dim(DM_Results)[1])*100,0),"%)","\n\n"))
  
  dev.off()
  
  cat(paste0("Save results"),"\n\n")
  
  tryCatch(saveRDS(DM_Results, file = file.path(s_bootstrap_folder,paste0("DM_Results_Limma_", mod_type, "_",loop_contrast, ".rds"))), error = function(e) warning(e$message))
  # tryCatch(fwrite(DM_Results,file = paste0(s_ROOT_dir,s_out_folder,"DM/DM_Results_",mod_type, "_", loop_contrast,".csv")), error = function(e) warning(e$message))
  
  rm(Limmaefit)
  gc()
  
  # return(DM_Results)
}



# Limma with weights function
DM_analysis_single_limmaWeights <- function(Metas_file, mod_type, scores_weights_file, pheno_file, 
                                            s_bootstrap_folder,s_CovFormulaImportantSV,temp_design) {
  
  # Load Metas
  Metas <- readRDS(Metas_file)
  Metas <- as.matrix(Metas)
  
  # Calculate weights
  scores_weights <- readRDS(scores_weights_file)
  scores_weights <- as.matrix(scores_weights)
  
  # Make LMmodel with weights
  cat(paste0("Fitting linear model"),"\n\n")
  temp_fit = limma::lmFit(as.matrix(Metas),temp_design, weights=scores_weights)
  
  # Show the coefs for first gene
  data.frame(temp_fit$coef[1,])
  
  cat(paste0("Rank genes using empirical Bayes method"),"\n\n")
  Limmaefit = limma::eBayes(temp_fit)
  
  # Save results limma
  # tryCatch(saveRDS(Limmaefit, file = file.path(paste0(s_ROOT_dir,s_out_folder,"DM/Limmaefit_LimmaW_", mod_type, ".rds"))), error = function(e) warning(e$message))
  
  loop_contrast = "GroupAD"
  s_coeff = loop_contrast 
  
  cat(paste0("Extract table of the top-ranked genes", "\n\n"))
  DM_Results = limma::topTable(Limmaefit, coef=s_coeff, num=Inf, sort.by = "P", adjust="BH") # Data of specified variables
  
  # Order based on p value
  # DM_Results_sign=DM_Results[DM_Results$P.Value<0.01,]
  print(head(DM_Results,20))
  
  destination <- file.path(s_bootstrap_folder, paste0("DM_hist_LimmaW_",mod_type, "_",loop_contrast, ".pdf"))
  pdf(file=destination)
  
  # histogram of p values and adjusted p values
  hist(DM_Results$P.Value ,xlim=c(0,1),main=s_coeff)
  hist(DM_Results$adj.P.Val ,xlim=c(0,1),main=s_coeff)
  min(DM_Results$adj.P.Val)
  
  
  # Descriptives
  cat(paste0("Minimum adjusted p value of ", s_coeff, ": " ,round(min(DM_Results$adj.P.Val),digits = 3), "\n\n"))
  
  cat(paste0("Contrast: ",s_coeff,".\n  number of features below TH (0.01)\n  Pvalue:",
             sum(DM_Results$P.Value<0.01)," (",round((sum(DM_Results$P.Value<0.01)/dim(DM_Results)[1])*100,0),"%)",
             "\n  AdjPvalue:",sum(DM_Results$adj.P.Val<0.05)," (",round((sum(DM_Results$adj.P.Val<0.05)/dim(DM_Results)[1])*100,0),"%)","\n\n"))
  
  dev.off()
  
  cat(paste0("Save results"),"\n\n")
  
  tryCatch(saveRDS(DM_Results, file = file.path(s_bootstrap_folder,paste0("DM_Results_LimmaW_", mod_type, "_",loop_contrast, ".rds"))), error = function(e) warning(e$message))
  # tryCatch(fwrite(DM_Results,file = paste0(s_ROOT_dir,s_out_folder,"DM/DM_Results_",mod_type, "_", loop_contrast,".csv")), error = function(e) warning(e$message))
  
  rm(Limmaefit)
  gc()
  
  # return(DM_Results)
}

# Creating weights function
weights_global_sample <- function(group_factor, coverage_file, pheno_file, output_file) {

  # Prepare list of groups
  group_weights_list <- list()
  unique_groups <- unique(group_factor)

  # Read scores file
  coverage <- readRDS(coverage_file)

  for (group in unique_groups) {
    # samples that belong to a specific group
    group_samples <- pheno_file$sample_names[pheno_file$Group == group]

    # scores for these specific samples
    group_scores <- coverage[,colnames(coverage) %in% group_samples]

    # scores divided by total sum of scores
    group_weights <- round(sweep(group_scores, 1, rowSums(coverage, na.rm=TRUE), `/`), digits = 4)

    # return weights
    group_weights_list[[group]] <- group_weights
  }

  # Combine group weights
  scores_weights <- do.call(cbind, group_weights_list)

  # remove prefix (AD, MCI, CTL)
  colnames(scores_weights) <- unlist(lapply(group_weights_list, function(x) colnames(x)))

  # Order weights based on Case ID
  scores_weights <- scores_weights[,colnames(coverage)]
  # print(head(scores_weights))

  # Save file
  saveRDS(scores_weights, file = output_file)

  return(scores_weights)

  # Free up memory
  rm(scores_weights,group_scores, group_weights)
  gc()

}


# DM_analysis_single_methylKit <- function(mod_type, counts_file, coverage, treatment_factor, s_bootstrap_folder){
#   
#   counts <- readRDS(counts_file)
#   
#   # Chr and position vectors
#   cat(paste0("Determining CpG positions", "\n\n"))
#   positions <- str_split_fixed(rownames(counts), "_", 2)
#   chr <- as.character(positions[,1])
#   pos <- as.numeric(positions[,2])
#   
#   idx_samples <- colnames(coverage)
#   
#   # Create methylRawlist
#   cat(paste0("Creating methylRaw data frames", "\n\n"))
#   methylRaw_list_subset <- pblapply(idx_samples, function(sample){
#     
#     df <- data.frame(
#       chr = chr,
#       start = pos,
#       end = pos,  # For CpG, start == end
#       strand = "+",
#       coverage = coverage[, sample],
#       numCs = counts[, sample],
#       numTs = coverage[, sample] - counts[, sample])
#     
#     methylRaw_obj <- new("methylRaw", df,
#                          sample.id = sample,
#                          assembly = "hg38",
#                          context = "CpG",
#                          resolution = "base")
#     
#     gc()
#     return(methylRaw_obj)
#   })
#   
#   rm(positions, chr, pos)
#   gc()
#   
#   # Create methyl Raw list object
#   cat(paste0("Creating methylRawList object", "\n\n"))
#   methylRaw_list_subset_sim <- methylRawList(methylRaw_list_subset, 
#                                              treatment = treatment_factor)
#   
#   # saveRDS(methylRaw_list_subset_sim, file = file.path(s_bootstrap_folder, paste0("methylRaw_list_",mod_type,".rds")))
#   
#   # Unite into one Base object
#   cat(paste0("Uniting methylRawList", "\n\n"))
#   methylBase_obj_sim <- methylKit::unite(methylRaw_list_subset_sim)
#   
#   # saveRDS(methylBase_obj_sim, file = file.path(s_bootstrap_folder, paste0("methylBase_obj_",mod_type,".rds")))
#   
#   rm(methylRaw_list_subset_sim)
#   gc()
#   
#   
#   # ---------------------------------------------------------------------------- # 
#   
#   # Modelling without correcting for overdispersion
#   cat(paste0("Running model: performing calculateDiffMeth", "\n\n"))
#   
#   DM_Results = calculateDiffMeth(methylBase_obj_sim, covariates = covariates, overdispersion = "none")
#   colnames(DM_Results) <- c("chr", "start", "end", "strand", "pvals", "fdrs", "meth.diff")
#   
#   s_coeff = "GroupAD"
#   
#   # Histogram of p values and adjusted p values
#   destination <- file.path(s_bootstrap_folder, paste0("DM_hist_Methylkit_", mod_type, ".pdf"))
#   pdf(file=destination)
#   
#   hist(DM_Results$pvals, xlim=c(0,1), main =s_coeff)
#   hist(DM_Results$fdrs, xlim=c(0,1), main =s_coeff)
#   
#   dev.off()
#   
#   # Descriptives
#   cat(paste0("Minimum adjusted p value of ", s_coeff, ": " ,round(min(DM_Results$fdrs),digits = 3), "\n\n"))
#   
#   cat(paste0("Contrast: ",s_coeff,".\n  number of features below TH (0.01)\n  Pvalue:",
#              sum(DM_Results$pvals<0.01, na.rm =TRUE)," (",round((sum(DM_Results$pvals<0.01, na.rm=TRUE)/dim(DM_Results)[1])*100,0),"%)",
#              "\n  AdjPvalue:",sum(DM_Results$fdrs<0.05, na.rm = TRUE)," (",round((sum(DM_Results$fdrs<0.05, na.rm = TRUE)/dim(DM_Results)[1])*100,0),"%)","\n\n"))
#   
#   # Saving DM Results
#   saveRDS(DM_Results, file = file.path(s_bootstrap_folder, paste0("DM_Results_methylkit_",mod_type,".rds")))
# }
# 
# 
# 
# # WITH CORRECTION
# DM_analysis_single_methylKit_MNcor <- function(mod_type, counts_file, coverage, treatment_factor, s_bootstrap_folder){
#   
#   counts <- readRDS(counts_file)
#   
#   # Chr and position vectors
#   cat(paste0("Determining CpG positions", "\n\n"))
#   positions <- str_split_fixed(rownames(counts), "_", 2)
#   chr <- as.character(positions[,1])
#   pos <- as.numeric(positions[,2])
#   
#   idx_samples <- colnames(coverage)
#   
#   # Create methylRawlist
#   cat(paste0("Creating methylRaw data frames", "\n\n"))
#   methylRaw_list_subset <- pblapply(idx_samples, function(sample){
#     
#     df <- data.frame(
#       chr = chr,
#       start = pos,
#       end = pos,  # For CpG, start == end
#       strand = "+",
#       coverage = coverage[, sample],
#       numCs = counts[, sample],
#       numTs = coverage[, sample] - counts[, sample])
#     
#     methylRaw_obj <- new("methylRaw", df,
#                          sample.id = sample,
#                          assembly = "hg38",
#                          context = "CpG",
#                          resolution = "base")
#     
#     gc()
#     return(methylRaw_obj)
#   })
#   
#   rm(positions, chr, pos)
#   gc()
#   
#   # Create methyl Raw list object
#   cat(paste0("Creating methylRawList object", "\n\n"))
#   methylRaw_list_subset_sim <- methylRawList(methylRaw_list_subset, 
#                                              treatment = treatment_factor)
#   
#   # saveRDS(methylRaw_list_subset_sim, file = file.path(s_bootstrap_folder, paste0("methylRaw_list_",mod_type,".rds")))
#   
#   # Unite into one Base object
#   cat(paste0("Uniting methylRawList", "\n\n"))
#   methylBase_obj_sim <- methylKit::unite(methylRaw_list_subset_sim)
#   
#   # saveRDS(methylBase_obj_sim, file = file.path(s_bootstrap_folder, paste0("methylBase_obj_",mod_type,".rds")))
#   
#   rm(methylRaw_list_subset_sim)
#   gc()
#   
#   
#   # ---------------------------------------------------------------------------- # 
#   
#   # Modelling with correction for overdispersion
#   cat(paste0("Running model: performing calculateDiffMeth with MN correction", "\n\n"))
#   
#   DM_Results_overd = calculateDiffMeth(methylBase_obj_sim, covariates = covariates, 
#                                        overdispersion = "MN", test = "Chisq")
#   
#   colnames(DM_Results_overd) <- c("chr", "start", "end", "strand", "pvals", "fdrs", "meth.diff")
#   
#   s_coeff = "GroupAD"
#   
#   # Histogram of p values and adjusted p values
#   destination <- file.path(s_bootstrap_folder, paste0("DM_hist_Methylkit_overd_", mod_type, ".pdf"))
#   pdf(file=destination)
#   
#   hist(DM_Results_overd$pvals, xlim=c(0,1), main =s_coeff)
#   hist(DM_Results_overd$fdrs, xlim=c(0,1), main =s_coeff)
#   
#   dev.off()
#   
#   # Descriptives
#   cat(paste0("Minimum adjusted p value of ", s_coeff, ": " ,round(min(DM_Results_overd$fdrs),digits = 3), "\n\n"))
#   
#   cat(paste0("Contrast: ",s_coeff,".\n  number of features below TH (0.01)\n  Pvalue:",
#              sum(DM_Results_overd$pvals<0.01, na.rm =TRUE)," (",round((sum(DM_Results_overd$pvals<0.01, na.rm=TRUE)/dim(DM_Results_overd)[1])*100,0),"%)",
#              "\n  AdjPvalue:",sum(DM_Results_overd$fdrs<0.05, na.rm = TRUE)," (",round((sum(DM_Results_overd$fdrs<0.05, na.rm = TRUE)/dim(DM_Results_overd)[1])*100,0),"%)","\n\n"))
#   
#   saveRDS(DM_Results_overd, file = file.path(s_bootstrap_folder, paste0("DM_Results_methylkit_MNcorrection_",mod_type,".rds")))
#   
# }



# WITH CORRECTION
DM_analysis_single_methylKit_combined <- function(mod_type, counts_file, coverage, treatment_factor, s_bootstrap_folder){
  
  counts <- readRDS(counts_file)
  
  # Chr and position vectors
  cat(paste0("Determining CpG positions", "\n\n"))
  positions <- str_split_fixed(rownames(counts), "_", 2)
  chr <- as.character(positions[,1])
  pos <- as.numeric(positions[,2])
  
  idx_samples <- colnames(coverage)
  
  # Create methylRawlist
  cat(paste0("Creating methylRaw data frames", "\n\n"))
  methylRaw_list_subset <- pblapply(idx_samples, function(sample){
    
    df <- data.frame(
      chr = chr,
      start = pos,
      end = pos,  # For CpG, start == end
      strand = "+",
      coverage = coverage[, sample],
      numCs = counts[, sample],
      numTs = coverage[, sample] - counts[, sample])
    
    methylRaw_obj <- new("methylRaw", df,
                         sample.id = sample,
                         assembly = "hg38",
                         context = "CpG",
                         resolution = "base")
    
    gc()
    return(methylRaw_obj)
  })
  
  rm(positions, chr, pos)
  gc()
  
  # Create methyl Raw list object
  cat(paste0("Creating methylRawList object", "\n\n"))
  methylRaw_list_subset_sim <- methylRawList(methylRaw_list_subset, 
                                             treatment = treatment_factor)
  
  # saveRDS(methylRaw_list_subset_sim, file = file.path(s_bootstrap_folder, paste0("methylRaw_list_",mod_type,".rds")))
  
  # Unite into one Base object
  cat(paste0("Uniting methylRawList", "\n\n"))
  methylBase_obj_sim <- methylKit::unite(methylRaw_list_subset_sim)
  
  # saveRDS(methylBase_obj_sim, file = file.path(s_bootstrap_folder, paste0("methylBase_obj_",mod_type,".rds")))
  
  rm(methylRaw_list_subset_sim)
  gc()
  
  
  
  # ---------------------------------------------------------------------------- # 
  
  # Modelling without correction for overdispersion
  cat(paste0("Running model: performing calculateDiffMeth", "\n\n"))
  
  DM_Results = calculateDiffMeth(methylBase_obj_sim, covariates = covariates, overdispersion = "none")
  colnames(DM_Results) <- c("chr", "start", "end", "strand", "pvals", "fdrs", "meth.diff")
  
  s_coeff = "GroupAD"
  
  # Descriptives
  cat(paste0("Minimum adjusted p value of ", s_coeff, ": " ,round(min(DM_Results$fdrs),digits = 3), "\n\n"))
  
  cat(paste0("Contrast: ",s_coeff,".\n  number of features below TH (0.01)\n  Pvalue:",
             sum(DM_Results$pvals<0.01, na.rm =TRUE)," (",round((sum(DM_Results$pvals<0.01, na.rm=TRUE)/dim(DM_Results)[1])*100,0),"%)",
             "\n  AdjPvalue:",sum(DM_Results$fdrs<0.05, na.rm = TRUE)," (",round((sum(DM_Results$fdrs<0.05, na.rm = TRUE)/dim(DM_Results)[1])*100,0),"%)","\n\n"))
  
  # Saving DM Results
  saveRDS(DM_Results, file = file.path(s_bootstrap_folder, paste0("DM_Results_methylkit_",mod_type,".rds")))
  
  rm(DM_Results)
  
  # ---------------------------------------------------------------------------- # 
  
  # Modelling with correction for overdispersion
  cat(paste0("Running model: performing calculateDiffMeth with MN correction", "\n\n"))
  
  DM_Results_overd = calculateDiffMeth(methylBase_obj_sim, covariates = covariates, 
                                       overdispersion = "MN", test = "Chisq")
  
  colnames(DM_Results_overd) <- c("chr", "start", "end", "strand", "pvals", "fdrs", "meth.diff")
  
  s_coeff = "GroupAD"
  
  # Descriptives
  cat(paste0("Minimum adjusted p value of ", s_coeff, ": " ,round(min(DM_Results_overd$fdrs),digits = 3), "\n\n"))
  
  cat(paste0("Contrast: ",s_coeff,".\n  number of features below TH (0.01)\n  Pvalue:",
             sum(DM_Results_overd$pvals<0.01, na.rm =TRUE)," (",round((sum(DM_Results_overd$pvals<0.01, na.rm=TRUE)/dim(DM_Results_overd)[1])*100,0),"%)",
             "\n  AdjPvalue:",sum(DM_Results_overd$fdrs<0.05, na.rm = TRUE)," (",round((sum(DM_Results_overd$fdrs<0.05, na.rm = TRUE)/dim(DM_Results_overd)[1])*100,0),"%)","\n\n"))
  
  saveRDS(DM_Results_overd, file = file.path(s_bootstrap_folder, paste0("DM_Results_methylkit_MNcorrection_",mod_type,".rds")))
  
}


# Differential 5-mC and 5-hmC locus analysis
DM_analysis_single_DSS <- function(mod_type, counts_file, coverage, s_bootstrap_folder,
                                   s_CovFormulaImportantSV,temp_design) {
  
  cat(paste0("Loading in data...", "\n\n"))
  # Load counts data (5-mC or 5-hmC)
  counts <- readRDS(counts_file)
  counts <- as.matrix(counts)
  
  coverage <- as.matrix(coverage)
  
  # Determine positions
  position <- str_split_fixed(rownames(counts), "_", n = 2)
  
  # Create BSobj
  cat(paste0("Creating BSobj...", "\n\n"))
  BSobj <- BSseq(chr = as.character(position[,1]), 
                 pos = as.integer(position[,2]), 
                 M = counts, 
                 Cov = coverage)
  
  # tryCatch(saveRDS(BSobj, file = file.path(s_bootstrap_folder, paste0("BSobj_",mod_type, ".rds"))), error = function(e) warning(e$message))
  
  # Free up memory
  rm(counts, position)
  gc()
  
  
  # Run DML
  cat(paste0("Running DMLfit...", "\n\n"))
  
  DMLfit = DMLfit.multiFactor(BSobj, design=pheno_subset_benchmark, formula=s_CovFormulaImportantSV, smoothing = FALSE)
  # tryCatch(saveRDS(DMLfit, file = file.path(s_bootstrap_folder,"DMLfit_",mod_type,".rds")), error = function(e) warning(e$message))
  
  # Contrast analysis
  cat(paste0("Contrast analysis on DMLfit"),"\n\n")
  
  loop_contrast <- "GroupAD"
  s_coeff = loop_contrast
  
  # Obtain statistics for contrasts
  
  DM_Results <- DMLtest.multiFactor(DMLfit, coef = s_coeff)
  print(head(DM_Results[order(DM_Results$pvals, decreasing = FALSE),]))
  
  # Histogram of p values and adjusted p values
  destination <- file.path(s_bootstrap_folder, paste0("DM_hist_DSS_",loop_contrast, mod_type, ".pdf"))
  pdf(file=destination)
  
  hist(DM_Results$pvals, xlim=c(0,1), main =s_coeff)
  hist(DM_Results$fdrs, xlim=c(0,1), main =s_coeff)
  
  dev.off()
  
  # Descriptives
  cat(paste0("Minimum adjusted p value of ", s_coeff, ": " ,round(min(DM_Results$fdrs),digits = 3), "\n\n"))
  
  cat(paste0("Contrast: ",s_coeff,".\n  number of features below TH (0.01)\n  Pvalue:",
             sum(DM_Results$pvals<0.01, na.rm =TRUE)," (",round((sum(DM_Results$pvals<0.01, na.rm=TRUE)/dim(DM_Results)[1])*100,0),"%)",
             "\n  AdjPvalue:",sum(DM_Results$fdrs<0.05, na.rm = TRUE)," (",round((sum(DM_Results$fdrs<0.05, na.rm = TRUE)/dim(DM_Results)[1])*100,0),"%)","\n\n"))
  
  tryCatch(saveRDS(DM_Results, file = file.path(s_bootstrap_folder, paste0("DM_Results_DSS_", mod_type, "_",loop_contrast, ".rds"))), error = function(e) warning(e$message))
  # tryCatch(fwrite(DM_Results,file = paste0(s_ROOT_dir,s_out_folder,"DM/DM_Results_",mod_type, "_", loop_contrast,".csv")), error = function(e) warning(e$message))
  
  
  cat(paste0("Save results"),"\n\n")
  
  rm(DMLfit)
  gc()
}


gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
  lambda
}


# Function to get cpg_ids
get_cpg_ids <- function(DM_Results) {
  
  cpg_ids <- if ("chr" %in% names(DM_Results) & "pos" %in% names(DM_Results))  {
    paste0(DM_Results$chr, "_", DM_Results$pos) # DSS
  } 
  else if ("chr" %in% names(DM_Results) & "start" %in% names(DM_Results))  {
    paste0(DM_Results$chr, "_", DM_Results$start) # MethylKit
  } 
  else if (identical(names(DM_Results), c("logFC", "AveExpr", "stat", "pvals", "fdrs", "B" ))) {
    rownames(DM_Results) # Limma
  } else {
    cat("Column names not having chr and pos ")
  }
  
  return(cpg_ids)
}


# Statistics
get_stats_txt <- function(data, depend_var, indep_var1, indep_var2) {
  
  if (is.null(indep_var2)) {
    
    file_output = file.path(s_OUT_dir, "QC", paste0("QC_stats_", depend_var, "_", indep_var1,".txt"))
    
    stats_df <- group_by(data, .data[[indep_var1]]) %>%
      summarise(
        N_set = length(.data[[depend_var]]),
        mean = round(mean(.data[[depend_var]], na.rm = TRUE), digits = 3),
        sd = round(sd(.data[[depend_var]], na.rm = TRUE), digits = 3),
        se = round(sd(.data[[depend_var]], na.rm = TRUE)/sqrt(length((.data[[depend_var]]))), digits= 3),
        median = round(median(.data[[depend_var]], na.rm = TRUE),digits =3 ),
        q25 = round(quantile(.data[[depend_var]], probs = 0.25, na.rm = TRUE), digits = 3),
        q75 = round(quantile(.data[[depend_var]], probs = 0.75, na.rm= TRUE), digits = 3))
    
    write.table(stats_df, file = file_output, quote = FALSE, sep = "\t")
    
    # One-way ANOVA with interaction
    # save model
    formula <- as.formula(paste(depend_var, "~", indep_var1))
    mod <- aov(formula, data = data)
    
    # Unbalanced designs
    # Type II ANOVA - no significant interaction
    # Type III ANOVA - with significant interaction
    
    # print results
    summary_df <- utils::capture.output(summary(mod))
    # write.table(summary_df, file = file_output, append = TRUE, quote = FALSE, sep = "\t")
    
    
    # Pairwise comparisons
    cat(paste0("Pairwise comparisons (Sidak): ~ ", indep_var1, "\n"))
    formula <- as.formula(paste("~", indep_var1))
    em <- emmeans(mod, formula)
    em_pairs <- pairs(em, adjust = "sidak")
    
    # write.table(em_pairs, file = file_output, append = TRUE, quote = FALSE, sep = "\t")
    
    em_pairs <- as.data.frame(em_pairs)

    
  } else {
    
    file_output = file.path(s_OUT_dir, "QC", paste0("QC_stats_", depend_var, "_", indep_var1, "_",indep_var2, ".txt"))
    
    # Interaction effect
    stats_df <- group_by(data, .data[[indep_var1]], .data[[indep_var2]]) %>%
      summarise(
        N_set = length(.data[[depend_var]]),
        mean = round(mean(.data[[depend_var]], na.rm = TRUE), digits = 3),
        sd = round(sd(.data[[depend_var]], na.rm = TRUE), digits = 3),
        se = round(sd(.data[[depend_var]], na.rm = TRUE)/sqrt(length((.data[[depend_var]]))), digits= 3),
        median = round(median(.data[[depend_var]], na.rm = TRUE),digits =3 ),
        q25 = round(quantile(.data[[depend_var]], probs = 0.25, na.rm = TRUE), digits = 3),
        q75 = round(quantile(.data[[depend_var]], probs = 0.75, na.rm= TRUE), digits = 3))
    
    write.table(stats_df, file = file_output, append = FALSE, quote = FALSE, sep = "\t")
    
    
    # Main effect independent variable 1
    stats_df <- group_by(data, .data[[indep_var1]]) %>%
      summarise(
        N_set = length(.data[[depend_var]]),
        mean = round(mean(.data[[depend_var]], na.rm = TRUE), digits = 3),
        sd = round(sd(.data[[depend_var]], na.rm = TRUE), digits = 3),
        se = round(sd(.data[[depend_var]], na.rm = TRUE)/sqrt(length((.data[[depend_var]]))), digits= 3),
        median = round(median(.data[[depend_var]], na.rm = TRUE),digits =3 ),
        q25 = round(quantile(.data[[depend_var]], probs = 0.25, na.rm = TRUE), digits = 3),
        q75 = round(quantile(.data[[depend_var]], probs = 0.75, na.rm= TRUE), digits = 3))
    
    write.table(stats_df, file = file_output, append = TRUE, quote = FALSE, sep = "\t")
    
    
    # Main effect independent variable 2
    stats_df <- group_by(data, .data[[indep_var2]]) %>%
      summarise(
        N_set = length(.data[[depend_var]]),
        mean = round(mean(.data[[depend_var]], na.rm = TRUE), digits = 3),
        sd = round(sd(.data[[depend_var]], na.rm = TRUE), digits = 3),
        se = round(sd(.data[[depend_var]], na.rm = TRUE)/sqrt(length((.data[[depend_var]]))), digits= 3),
        median = round(median(.data[[depend_var]], na.rm = TRUE),digits =3 ),
        q25 = round(quantile(.data[[depend_var]], probs = 0.25, na.rm = TRUE), digits = 3),
        q75 = round(quantile(.data[[depend_var]], probs = 0.75, na.rm= TRUE), digits = 3))
    
    write.table(stats_df, file = file_output, append = TRUE, quote = FALSE, sep = "\t")
    
    # Two-way ANOVA with interaction
    # save model
    
    formula <- as.formula(paste(depend_var, "~", indep_var1, "*", indep_var2))
    mod <- aov(formula, data = data)
    
    # Unbalanced designs
    # Type II ANOVA - no significant interaction
    # Type III ANOVA - with significant interaction
    
    Anova_df <- utils::capture.output(Anova(mod, type = "II"))
    write.table(Anova_df, file = file_output, append = TRUE, quote = FALSE, sep = "\t")
    
    # print results
    summary_df <- utils::capture.output(summary(mod))
    write.table(summary_df, file = file_output, append = TRUE, quote = FALSE, sep = "\t")
    
    
    # Pairwise comparisons
    cat(paste0("Pairwise comparisons (Sidak): ~ ", indep_var2, " * ", indep_var1, "\n"))
    formula <- as.formula(paste("~", indep_var2, "|", indep_var1))
    
    em <- emmeans(mod, formula)  # "group" means within each "set"
    em_pairs <- utils::capture.output(pairs(em, adjust = "sidak"))
    
    write.table(em_pairs, file = file_output, append = TRUE, quote = FALSE, sep = "\t")    
    
    
    # Pairwise comparisons
    cat(paste0("Pairwise comparisons (Sidak): ~ ", indep_var1, " * ", indep_var2, "\n"))
    formula <- as.formula(paste("~", indep_var1, "|", indep_var2))
    
    em <- emmeans(mod, formula)  # "set" means within each "group"
    em_pairs <- utils::capture.output(pairs(em, adjust = "sidak"))
    
    write.table(em_pairs, file = file_output, append = TRUE, quote = FALSE, sep = "\t")    
    
    
    # Main variable
    cat(paste0("Pairwise comparisons (Sidak): ~ ", indep_var1, "\n"))
    formula <- as.formula(paste("~", indep_var1))
    
    em <- emmeans(mod, formula)
    em_pairs <- utils::capture.output(pairs(em, adjust = "sidak"))
    
    write.table(em_pairs, file = file_output, append = TRUE, quote = FALSE, sep = "\t")    
    
    # Data for plot
    # Pairwise comparisons
    em_pairs <- pairs(em, adjust = "sidak")
    em_pairs <- as.data.frame(em_pairs)
  }
  return(em_pairs)
    
}

# convert p values to asterix and NS
convert_p_to_asterix <- function(values) {
  
  p_signifs <- lapply(values, function(value) {
    
    p_signif <- if (value <= 0.001) 
    {"***"} else if (value <= 0.01) 
    {"**"} else if (value <= 0.05) 
    {"*"} else if (value > 0.05) {"NS"}
    return(p_signif)
  })
  
  p_signifs <- do.call(rbind, p_signifs)
  return(p_signifs)
  
}



plot_boxplots <- function(data, x_var, y_var, group_var, pheno_file) {
  
  
  # Long format sequencing depths
  data_long <- data %>% pivot_longer(cols = everything(), names_to = x_var,
                                                       values_to = y_var)
  
  data_long <- left_join(data_long, pheno_file[,c(x_var, group_var)])
  data_long <- data_long[complete.cases(data_long),]
  
  
  # Plot sequencing depths distribution per sample
  my_colors <- c("#084594","#2171B5","#6BAED6", "#718FAB", "#343434", "#C6DBEF")
  my_shapes <- c(16, 15, 17, 18)
  
  destination <- file.path(s_OUT_dir, paste0("Plots/Bootstrapping_", x_var,"_",y_var ,".svg"))
  svglite(file=destination, width = 18, height = 10)
  
  image <- ggplot(data, aes(x =  reorder(.data[[x_var]], .data[[y_var]], FUN=median), y = .data[[y_var]],
                                          fill = .data[[group_var]], shape = .data[[group_var]])) +
    geom_boxplot(color = "black", outlier.shape = NA, coef = 0) +
    coord_cartesian(ylim =  c(0, 100)) +
    theme(legend.position = "top",
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank()) +
    
    xlab("Samples") +
    ylab("DNA methylation (%)") +
    scale_fill_manual(name = group_var, values = c("#084594","#6BAED6"), labels = label_text) +
    scale_shape_manual(name =group_var, values = my_shapes, labels = label_text)
  
  print(image)
  
  dev.off()
  
}


#-----------------------------------------------------------------------------------------------------#
#					1. Select samples                                                       ----
#-----------------------------------------------------------------------------------------------------#

# Get groups and sample names
groups <- pheno_subset_pass$Group
sample_names <- pheno_subset_pass$sample_names

# samples that belong to a specific cohort
group_samples <- list()

# Get sample names per group
for (group in unique(groups)) {
  group_samples[[group]] <- sample_names[groups == group]
}

# Subset to Control and AD samples only
group_samples <- group_samples[c("CTL", "AD")]
idx_samples <- unlist(group_samples)

idx_samples_CTL <- group_samples$CTL
idx_samples_AD <- group_samples$AD

rm(groups, group_samples, sample_names)


#-----------------------------------------------------------------------------------------------------#
#					2. Load subset data to compare                                                           ----
#-----------------------------------------------------------------------------------------------------#


# Pheno subset for benchmarking
pheno_subset_benchmark <- pheno_subset_pass[idx_samples,]
pheno_subset_benchmark$Group <- factor(pheno_subset_benchmark$Group, levels = c("CTL", "AD"))
pheno_subset_benchmark$Gender <- factor(pheno_subset_benchmark$Gender)
pheno_subset_benchmark$Flowcell <- factor(pheno_subset_benchmark$Flowcell)

# Treatment factor coded 0 to 1
treatment_factor <- factor(pheno_subset_benchmark$Group, labels = c(0,1), levels = c("CTL", "AD"))

# Covariates
covariates <- pheno_subset_benchmark[,c("Age", "Gender", "PMI", "Flowcell")]
covariates$Flowcell<- factor(covariates$Flowcell)
covariates$Gender <- factor(covariates$Gender)

mod_type = "methyl"


rm(pheno_subset_pass)
gc()

#-----------------------------------------------------------------------------------------------------#
#						3. Select 200000 CpG sites to be tested                                                     ----
#-----------------------------------------------------------------------------------------------------#

# Load scores
scores <- as.matrix(readRDS(file = paste0(s_ROOT_dir, s_out_folder,"Bedmethyl_all/bedfiles_scores_all_5X_pass_nosexchr.rds")))
betas <-  as.matrix(readRDS(file = paste0(s_ROOT_dir, s_out_folder,"Bedmethyl_all/bedfiles_methyl_all_5X_pass_nosexchr.rds")))

# Randomly select cpg sites
n_cpg_sites = 210000

# Randomly select cpg sites

set.seed(72)
idx_cpgs <- sample(nrow(scores), n_cpg_sites)


# Subset data
scores_subset <- as.data.frame(scores[idx_cpgs, c(idx_samples_CTL, idx_samples_AD)])
betas_subset <- as.data.frame(betas[idx_cpgs, c(idx_samples_CTL, idx_samples_AD)])

saveRDS(scores_subset, file = paste0(s_ROOT_dir, s_out_folder,"Bedmethyl_all/bedfiles_scores_all_5X_pass_nosexchr_subset_210.rds"))
saveRDS(betas_subset, file = paste0(s_ROOT_dir, s_out_folder,"Bedmethyl_all/bedfiles_methyl_all_5X_pass_nosexchr_subset_210.rds"))


# Long format sequencing depths
scores_subset_long <- scores_subset %>% pivot_longer(cols = everything(), names_to = "sample_names",
                                                   values_to = "Coverage")

scores_subset_long <- left_join(scores_subset_long, pheno_subset_benchmark[,c("sample_names", "Group")])
scores_subset_long <- scores_subset_long[complete.cases(scores_subset_long),]


# Plot sequencing depths distribution per sample
destination <- file.path(s_OUT_dir, paste0("Plots/Bootstrapping_Coverage_Sample_names_210.svg"))
svglite(file=destination, width = 18, height = 10)

image <- ggplot(scores_subset_long, aes(x =  reorder(sample_names, Coverage, FUN=median), y = Coverage,
                          fill = Group, shape = Group)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  coord_cartesian(ylim =  c(0, 45)) +
  theme(legend.position = "top",
         axis.text.x = element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("Samples") +
  ylab("Sequencing depth (X) per CpG site") +
  scale_fill_manual(name = "Group", values = c("#084594","#6BAED6"), labels = c("CTL", "AD"))

print(image)

dev.off()



# Long format sequencing depths
betas_subset_long <- betas_subset %>% pivot_longer(cols = everything(), names_to = "sample_names",
                                                     values_to = "Methylation")

betas_subset_long <- left_join(betas_subset_long, pheno_subset_benchmark[,c("sample_names", "Group")])
betas_subset_long <- betas_subset_long[complete.cases(betas_subset_long),]


# Plot sequencing depths distribution per sample
destination <- file.path(s_OUT_dir, paste0("Plots/Bootstrapping_Methylation_Sample_names_210.svg"))
svglite(file=destination, width = 18, height = 10)

image <- ggplot(betas_subset_long, aes(x =  reorder(sample_names, Methylation, FUN=median), y = Methylation,
                                        fill = Group, shape = Group)) +
  geom_boxplot(color = "black", outlier.shape = NA, coef = 0) +
  coord_cartesian(ylim =  c(0, 100)) +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank()) +
  
  xlab("Samples") +
  ylab("DNA methylation (%)") +
  scale_fill_manual(name = "Group", values = c("#084594","#6BAED6"), labels = c("CTL", "AD"))

print(image)

dev.off()

# Free up memory
rm(scores, betas, scores_subset_long, betas_subset_long)
gc()



#-----------------------------------------------------------------------------------------------------#
#					4. Zero change                                                        ----
#-----------------------------------------------------------------------------------------------------#
# Scores weights
weights_global_sample(pheno_subset_benchmark$Group, 
                      file.path(s_OUT_dir, "Bedmethyl_all", "bedfiles_scores_all_5X_pass_nosexchr_subset_210.rds"),
                      pheno_subset_benchmark,
                      output_file = file.path(s_OUT_dir,"QC", "scores_group_weights_benchmarking_pre.rds"))



# Dataset (Base)
betas_subset <- readRDS(file.path(s_OUT_dir, "Bedmethyl_all", "bedfiles_methyl_all_5X_pass_nosexchr_subset_210.rds"))
coverage = readRDS(file.path(s_OUT_dir, "Bedmethyl_all", "bedfiles_scores_all_5X_pass_nosexchr_subset_210.rds"))


# settings for the bootstrapping
n_samples <- ncol(betas_subset)
n_true_cpgs <- 500  # Number of true DMPs to spike in
effect_size <- 0  # e.g., 20% methylation difference

pblapply(4, function(repetition) {
  
  # New directory 
  s_bootstrap_folder <- paste0(s_OUT_dir, "DM/","bootstrap_pre_ES", effect_size, "_", repetition)
  dir.create(s_bootstrap_folder, showWarnings = FALSE)
  
  true_diff_cpgs <- sample(1:length(idx_cpgs), n_true_cpgs)
  saveRDS(true_diff_cpgs, file = file.path(s_bootstrap_folder, paste0("true_diff_cpgs",repetition,".rds")))
  
  # Inject effect
  betas_sim <- betas_subset
  
  for (i in true_diff_cpgs) {
    if (runif(1) > 0.5) {
      betas_sim[i, idx_samples_AD] <- pmin(betas_subset[i, idx_samples_AD] + effect_size, 100)
    } else {
      betas_sim[i, idx_samples_AD] <- pmax(betas_subset[i, idx_samples_AD] - effect_size, 0)
    }
  }
  
  # Convert to M values
  Metas_sim <- mod2M(as.data.frame(imputezerohundred(betas_sim)))
  
  # Convert to count values
  counts_sim <- round(betas_sim * coverage / 100)
  counts_sim[is.na(counts_sim)] <- 0
  
  # Save objects
  saveRDS(betas_sim, file = file.path(s_bootstrap_folder, paste0("betas_subset_bootstrap",repetition,".rds")))
  saveRDS(Metas_sim, file = file.path(s_bootstrap_folder, paste0("Metas_subset_bootstrap",repetition,".rds")))
  saveRDS(counts_sim, file = file.path(s_bootstrap_folder, paste0("counts_subset_bootstrap",repetition,".rds")))
  
  rm(betas_sim, Metas_sim, counts_sim)
  gc()
  
  
  # ~ Methylkit combined ----
  
  DM_analysis_single_methylKit_combined(mod_type = mod_type,
                                        s_bootstrap_folder = s_bootstrap_folder,
                                        counts_file = file.path(s_bootstrap_folder, paste0("counts_subset_bootstrap",repetition,".rds")),
                                        coverage = coverage,
                                        treatment_factor = treatment_factor)
  
  
  #-----------------------------------------------------------------------------------------------------#
  # ~ DSS ----
  
  # Make design matrix formula; correcting for Subject, determining the effect of Tissue
  s_CovFormulaImportantSV = ~ Group + Age + factor(Gender) + PMI + factor(Flowcell)
  
  # make design matrix
  temp_design = model.matrix(as.formula(s_CovFormulaImportantSV),data=pheno_subset_benchmark)
  colnames(temp_design) = make.names(colnames(temp_design))
  head(temp_design)
  
  # Run the DM analysis
  DM_analysis_single_DSS(mod_type = mod_type,
                         counts_file = file.path(s_bootstrap_folder, paste0("counts_subset_bootstrap",repetition,".rds")),
                         coverage = coverage,
                         s_bootstrap_folder = s_bootstrap_folder,
                         s_CovFormulaImportantSV = s_CovFormulaImportantSV,
                         temp_design = temp_design)
  
  #-----------------------------------------------------------------------------------------------------#
  # ~ Limma regular ----
  
  # Make design matrix formula; correcting for Subject, determining the effect of Tissue
  s_CovFormulaImportantSV = ~ Group + Age + factor(Gender) + PMI + factor(Flowcell)
  
  # make design matrix
  temp_design = model.matrix(as.formula(s_CovFormulaImportantSV),data=pheno_subset_benchmark)
  colnames(temp_design) = make.names(colnames(temp_design))
  head(temp_design)
  
  # Running analysis
  DM_analysis_single_limma(Metas_file = file.path(s_bootstrap_folder, paste0("Metas_subset_bootstrap",repetition,".rds")),
                           mod_type = mod_type,
                           s_bootstrap_folder = s_bootstrap_folder,
                           s_CovFormulaImportantSV = s_CovFormulaImportantSV,
                           temp_design = temp_design)
  
  gc()
  
  #-----------------------------------------------------------------------------------------------------# 
  # ~ Limma with weights ----
  
  # Make design matrix formula; correcting for Subject, determining the effect of Tissue
  s_CovFormulaImportantSV = ~ Group + Age + factor(Gender) + PMI + factor(Flowcell)
  
  # make design matrix
  temp_design = model.matrix(as.formula(s_CovFormulaImportantSV),data=pheno_subset_benchmark)
  colnames(temp_design) = make.names(colnames(temp_design))
  head(temp_design)
  
  DM_analysis_single_limmaWeights(file.path(s_bootstrap_folder, paste0("Metas_subset_bootstrap",repetition,".rds")),
                                  mod_type = mod_type, 
                                  # coverage_file = file.path(s_OUT_dir, "Bedmethyl_all", "bedfiles_scores_all_5X_pass_nosexchr_subset.rds"), 
                                  # group_factor = pheno_subset_benchmark$Group, 
                                  scores_weights_file = file.path(s_OUT_dir,"QC", "scores_group_weights_benchmarking_pre.rds"),
                                  pheno_file = pheno_subset_benchmark,
                                  s_bootstrap_folder = s_bootstrap_folder,
                                  s_CovFormulaImportantSV = s_CovFormulaImportantSV,
                                  temp_design = temp_design)
  
  gc()
  
  
  #-----------------------------------------------------------------------------------------------------#
  #  ~ List of DM Results ----
  temp = list.files(file.path(s_bootstrap_folder), pattern="DM_Results_",
                    full.names = TRUE)
  
  file_base_name <- sub(".rds", "", sub("DM_Results_", "", basename(temp)))
  
  DM_Results_list = lapply(temp, function(file){
    DM_Results <- readRDS(file)
    
    if(identical(names(DM_Results), c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B" ))) {
      colnames(DM_Results) <- c("logFC", "AveExpr", "stat", "pvals", "fdrs", "B" )
    } 
    
    return(DM_Results)
  })
  
  
  names(DM_Results_list) <- file_base_name
  
  saveRDS(DM_Results_list, file.path(s_OUT_dir, "BootstrapSum", paste0("DM_Results_list_bootstrap",repetition,".rds")))
  
  #-----------------------------------------------------------------------------------------------------#
  #  ~ TPR and FDR ----
  
  # True positives
  true_diff_cpgs_names <- rownames(coverage)[true_diff_cpgs]
  
  # True negatives
  true_neg_cpg_names <- rownames(coverage)[-c(true_diff_cpgs)]
  
  temp = names(DM_Results_list)
  
  fdr_stats <- lapply(temp, function(DM_Results_method){
    
    # Load DM Results
    DM_Results <- DM_Results_list[[DM_Results_method]]
    
    na_values <- sum(is.na(DM_Results$fdrs))
    
    # Get Significant results
    DM_Results_sign <- DM_Results[DM_Results$fdrs < 0.05, ]
    cpg_ids_sign <- get_cpg_ids(DM_Results_sign)
    
    
    # Get not significant Cpg sites
    DM_Results_notsign <- DM_Results[DM_Results$fdrs >= 0.05, ]
    cpg_ids_notsign <- get_cpg_ids(DM_Results_notsign)
    
    
    # True positives
    true_pos <- length(intersect(cpg_ids_sign, true_diff_cpgs_names))
    
    # False positives
    false_pos <- length(setdiff(cpg_ids_sign, true_diff_cpgs_names))
    
    # True negatives
    true_neg <- length(intersect(cpg_ids_notsign, true_neg_cpg_names))
    
    # False negatives
    false_neg <- length(setdiff(cpg_ids_notsign, true_neg_cpg_names))
    
    # True Positive Rate (TPR)
    TPR <- true_pos / (true_pos + false_neg)
    
    # False discovery rate (FDR)
    FDR <- false_pos / (true_pos+ false_pos)
    
    # Lambda value qqplot
    lambda <-median(qchisq(1-DM_Results$pvals,df= 1))/qchisq(0.5,df = 1)
    
    # Combine into lambda
    stats <- data.frame(true_pos, false_pos, true_neg, false_neg, TPR, FDR, na_values, lambda)
    return(stats)
  })
  
  names(fdr_stats) <- names(DM_Results_list)
  fdr_stats <- do.call("rbind", fdr_stats)
  
  saveRDS(fdr_stats, file.path(s_OUT_dir, "BootstrapSum", paste0("Bootstrap_results_", repetition, ".rds")))
  write.table(fdr_stats, file.path(s_OUT_dir, "BootstrapSum", paste0("Bootstrap_results_", repetition, ".txt")), 
              sep = "\t", quote = FALSE)
  
  rm(DM_Results_list)
  gc()
  
})


#-----------------------------------------------------------------------------------------------------#
#						5. Extract non-significant CpG sites                                                 ----
#-----------------------------------------------------------------------------------------------------#


# DM Results list
DM_Results_list <- readRDS(file = file.path(s_ROOT_dir, s_out_folder, "BootstrapSum", "0_perc_210", "DM_Results_list_bootstrap4.rds"))

temp = names(DM_Results_list) 

significant_cpgs <- pblapply(temp, function(method){
  
  DM_Results <- DM_Results_list[[method]]
  
  # Significant results
  DM_Results_sign <- DM_Results[DM_Results$fdrs < 0.05,]
  
  print(DM_Results_sign)
  
  if (nrow(DM_Results_sign) > 0) {
    top_probes <- if ("start" %in% colnames(DM_Results_sign)) {
      
      paste0(DM_Results_sign$chr, "_", DM_Results_sign$start)
      
    } else if ("pos" %in% colnames(DM_Results_sign)) {
      paste0(DM_Results_sign$chr, "_", DM_Results_sign$pos)
      
    } else {
      rownames(DM_Results_sign)
    }

  } else {NULL}
  
  })

significant_cpgs <- unique(unlist(significant_cpgs))

non_significant_cpgs <- rownames(coverage)[!rownames(coverage) %in% significant_cpgs]
saveRDS(non_significant_cpgs, file = paste0(s_ROOT_dir, s_out_folder,"QC", "non_significant_cpgs_210.rds"))

rm(DM_Results, DM_Results_list, DM_Results_sign)


#-----------------------------------------------------------------------------------------------------#
#						6. Select 200000 CpG sites                                                      ----
#-----------------------------------------------------------------------------------------------------#

# Load scores
scores <- as.matrix(readRDS(file = paste0(s_ROOT_dir, s_out_folder,"Bedmethyl_all/bedfiles_scores_all_5X_pass_nosexchr.rds")))
betas <- as.matrix(readRDS(file = paste0(s_ROOT_dir, s_out_folder,"Bedmethyl_all/bedfiles_methyl_all_5X_pass_nosexchr.rds")))


# Non significant CpGs, previously tested
non_significant_cpgs <- readRDS(file = paste0(s_ROOT_dir, s_out_folder,"QC", "non_significant_cpgs_210.rds"))

# Randomly select cpg sites
n_cpg_sites = 200000

set.seed(72)
idx_cpgs <- sample(non_significant_cpgs, n_cpg_sites)


# Subset data
scores_subset <- as.data.frame(scores[idx_cpgs, c(idx_samples_CTL, idx_samples_AD)])
betas_subset <- as.data.frame(betas[idx_cpgs, c(idx_samples_CTL, idx_samples_AD)])

saveRDS(scores_subset, file = paste0(s_ROOT_dir, s_out_folder,"Bedmethyl_all/bedfiles_scores_all_5X_pass_nosexchr_subset.rds"))
saveRDS(betas_subset, file = paste0(s_ROOT_dir, s_out_folder,"Bedmethyl_all/bedfiles_methyl_all_5X_pass_nosexchr_subset.rds"))

rm(betas, scores)
gc()

# Long format sequencing depths
scores_subset_long <- scores_subset %>% pivot_longer(cols = everything(), names_to = "sample_names",
                                                     values_to = "Coverage")

scores_subset_long <- left_join(scores_subset_long, pheno_subset_benchmark[,c("sample_names", "Group")])
scores_subset_long <- scores_subset_long[complete.cases(scores_subset_long),]


# Plot sequencing depths distribution per sample
destination <- file.path(s_OUT_dir, paste0("Plots/Bootstrapping_Coverage_Sample_names_200_new.svg"))
svglite(file=destination, width = 18, height = 10)

image <- ggplot(scores_subset_long, aes(x =  reorder(sample_names, Coverage, FUN=median), y = Coverage,
                                        fill = Group, shape = Group)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  coord_cartesian(ylim =  c(0, 45)) +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("Samples") +
  ylab("Sequencing depth (X) per CpG site") +
  scale_fill_manual(name = "Group", values = c("#084594","#6BAED6"), labels = c("CTL", "AD"))

print(image)

dev.off()



# Long format sequencing depths
betas_subset_long <- betas_subset %>% pivot_longer(cols = everything(), names_to = "sample_names",
                                                   values_to = "Methylation")

betas_subset_long <- left_join(betas_subset_long, pheno_subset_benchmark[,c("sample_names", "Group")])
betas_subset_long <- betas_subset_long[complete.cases(betas_subset_long),]


# Plot sequencing depths distribution per sample
destination <- file.path(s_OUT_dir, paste0("Plots/Bootstrapping_Methylation_Sample_names_200_new.svg"))
svglite(file=destination, width = 18, height = 10)

image <- ggplot(betas_subset_long, aes(x =  reorder(sample_names, Methylation, FUN=median), y = Methylation,
                                       fill = Group, shape = Group)) +
  geom_boxplot(color = "black", outlier.shape = NA, coef = 0) +
  coord_cartesian(ylim =  c(0, 100)) +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank()) +
  
  xlab("Samples") +
  ylab("DNA methylation (%)") +
  scale_fill_manual(name = "Group", values = c("#084594","#6BAED6"), labels = c("CTL", "AD"))

print(image)

dev.off()


# Free up memory
rm(scores_subset, betas_subset, scores_subset_long, betas_subset_long, non_significant_cpgs , image)
gc()


#-----------------------------------------------------------------------------------------------------#
#					7. Impute differences and DM analysis (first check zero)                                                       ----
#-----------------------------------------------------------------------------------------------------#

# temp = list.files(file.path(s_OUT_dir, "Bedmethyl_all/"), pattern=paste0(cutoff, "X_pass_nosexchr_subset"),
#                   full.names = TRUE)
# temp = str_subset(temp, "count|scores", negate=TRUE)
# betas_subset <- readRDS(temp[2])


# Dataset (Base)
betas_subset <- readRDS(file.path(s_OUT_dir, "Bedmethyl_all", "bedfiles_methyl_all_5X_pass_nosexchr_subset.rds"))
coverage = readRDS(file.path(s_OUT_dir, "Bedmethyl_all", "bedfiles_scores_all_5X_pass_nosexchr_subset.rds"))

# Scores weights
weights_global_sample(pheno_subset_benchmark$Group, 
                      file.path(s_OUT_dir, "Bedmethyl_all", "bedfiles_scores_all_5X_pass_nosexchr_subset.rds"),
                      pheno_subset_benchmark,
                      output_file = file.path(s_OUT_dir,"QC", "scores_group_weights_benchmarking.rds"))


# settings for the bootstrapping
n_samples <- ncol(betas_subset)
n_true_cpgs <- 500  # Number of true DMPs to spike in
effect_size <- 5  # e.g., 20% methylation difference

# How many times
n_bootstrap <- 50

gc()

pblapply(1:n_bootstrap, function(repetition) {
  
  # New directory 
  s_bootstrap_folder <- paste0(s_OUT_dir, "DM/","bootstrap_ES", effect_size, "_", repetition)
  dir.create(s_bootstrap_folder, showWarnings = FALSE)

  true_diff_cpgs <- sample(1:length(idx_cpgs), n_true_cpgs)
  saveRDS(true_diff_cpgs, file = file.path(s_bootstrap_folder, paste0("true_diff_cpgs",repetition,".rds")))
  
  # Inject effect
  betas_sim <- betas_subset
  
  for (i in true_diff_cpgs) {
    if (runif(1) > 0.5) {
      betas_sim[i, idx_samples_AD] <- pmin(betas_subset[i, idx_samples_AD] + effect_size, 100)
    } else {
      betas_sim[i, idx_samples_AD] <- pmax(betas_subset[i, idx_samples_AD] - effect_size, 0)
    }
  }
  
  # Convert to M values
  Metas_sim <- mod2M(as.data.frame(imputezerohundred(betas_sim)))

  # Convert to count values
  counts_sim <- round(betas_sim * coverage / 100)
  counts_sim[is.na(counts_sim)] <- 0

  # Save objects
  saveRDS(betas_sim, file = file.path(s_bootstrap_folder, paste0("betas_subset_bootstrap",repetition,".rds")))
  saveRDS(Metas_sim, file = file.path(s_bootstrap_folder, paste0("Metas_subset_bootstrap",repetition,".rds")))
  saveRDS(counts_sim, file = file.path(s_bootstrap_folder, paste0("counts_subset_bootstrap",repetition,".rds")))
  
  rm(betas_sim, Metas_sim, counts_sim)
  gc()
  
  # 
  # #-----------------------------------------------------------------------------------------------------#
  # # ~ Methylkit MN cor ----
  # 
  # DM_analysis_single_methylKit_MNcor(mod_type = mod_type,
  #                                      s_bootstrap_folder = s_bootstrap_folder,
  #                                      counts_file = file.path(s_bootstrap_folder, paste0("counts_subset_bootstrap",repetition,".rds")),
  #                                      coverage = coverage,
  #                                      treatment_factor = treatment_factor)
  # 
  # #-----------------------------------------------------------------------------------------------------#
  # # ~ Methylkit  ----
  # 
  # DM_analysis_single_methylKit(mod_type = mod_type,
  #                                s_bootstrap_folder = s_bootstrap_folder,
  #                                counts_file = file.path(s_bootstrap_folder, paste0("counts_subset_bootstrap",repetition,".rds")),
  #                                coverage = coverage,
  #                                treatment_factor = treatment_factor)
  # 
  
  # ~ Methylkit combined ----

  DM_analysis_single_methylKit_combined(mod_type = mod_type,
                                       s_bootstrap_folder = s_bootstrap_folder,
                                       counts_file = file.path(s_bootstrap_folder, paste0("counts_subset_bootstrap",repetition,".rds")),
                                       coverage = coverage,
                                       treatment_factor = treatment_factor)


  #-----------------------------------------------------------------------------------------------------#
  # ~ DSS ----

  # Make design matrix formula; correcting for Subject, determining the effect of Tissue
  s_CovFormulaImportantSV = ~ Group + Age + factor(Gender) + PMI + factor(Flowcell)

  # make design matrix
  temp_design = model.matrix(as.formula(s_CovFormulaImportantSV),data=pheno_subset_benchmark)
  colnames(temp_design) = make.names(colnames(temp_design))
  head(temp_design)

  # Run the DM analysis
  # coverage <- readRDS(file.path(s_OUT_dir, paste0("Bedmethyl_all/", "bedfiles_scores_all_", cutoff, "X_pass_nosexchr_subset.rds")))
  # coverage <- as.matrix(coverage)

  DM_analysis_single_DSS(mod_type = mod_type,
                           counts_file = file.path(s_bootstrap_folder, paste0("counts_subset_bootstrap",repetition,".rds")),
                           coverage = coverage,
                           s_bootstrap_folder = s_bootstrap_folder,
                           s_CovFormulaImportantSV = s_CovFormulaImportantSV,
                           temp_design = temp_design)

  #-----------------------------------------------------------------------------------------------------#
  # ~ Limma regular ----

  # Make design matrix formula; correcting for Subject, determining the effect of Tissue
  s_CovFormulaImportantSV = ~ Group + Age + factor(Gender) + PMI + factor(Flowcell)

  # make design matrix
  temp_design = model.matrix(as.formula(s_CovFormulaImportantSV),data=pheno_subset_benchmark)
  colnames(temp_design) = make.names(colnames(temp_design))
  head(temp_design)

  # Running analysis
  DM_analysis_single_limma(Metas_file = file.path(s_bootstrap_folder, paste0("Metas_subset_bootstrap",repetition,".rds")),
                             mod_type = mod_type,
                             s_bootstrap_folder = s_bootstrap_folder,
                             s_CovFormulaImportantSV = s_CovFormulaImportantSV,
                             temp_design = temp_design)

  gc()
  
  #-----------------------------------------------------------------------------------------------------# 
  # ~ Limma with weights ----

  # Make design matrix formula; correcting for Subject, determining the effect of Tissue
  s_CovFormulaImportantSV = ~ Group + Age + factor(Gender) + PMI + factor(Flowcell)

  # make design matrix
  temp_design = model.matrix(as.formula(s_CovFormulaImportantSV),data=pheno_subset_benchmark)
  colnames(temp_design) = make.names(colnames(temp_design))
  head(temp_design)

  DM_analysis_single_limmaWeights(file.path(s_bootstrap_folder, paste0("Metas_subset_bootstrap",repetition,".rds")),
                                    mod_type = mod_type, 
                                    # coverage_file = file.path(s_OUT_dir, "Bedmethyl_all", "bedfiles_scores_all_5X_pass_nosexchr_subset.rds"), 
                                    # group_factor = pheno_subset_benchmark$Group, 
                                    scores_weights_file = file.path(s_OUT_dir,"QC", "scores_group_weights_benchmarking.rds"),
                                    pheno_file = pheno_subset_benchmark,
                                    s_bootstrap_folder = s_bootstrap_folder,
                                    s_CovFormulaImportantSV = s_CovFormulaImportantSV,
                                    temp_design = temp_design)

  gc()
  

  #-----------------------------------------------------------------------------------------------------#
  #  ~ List of DM Results ----
  temp = list.files(file.path(s_bootstrap_folder), pattern="DM_Results_",
                    full.names = TRUE)
  
  file_base_name <- sub(".rds", "", sub("DM_Results_", "", basename(temp)))
  
  DM_Results_list = lapply(temp, function(file){
    DM_Results <- readRDS(file)
    
    if(identical(names(DM_Results), c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B" ))) {
      colnames(DM_Results) <- c("logFC", "AveExpr", "stat", "pvals", "fdrs", "B" )
    } 
    
    return(DM_Results)
  })
  
  
  names(DM_Results_list) <- file_base_name
  
  saveRDS(DM_Results_list, file.path(s_OUT_dir, "BootstrapSum", paste0("DM_Results_list_bootstrap",repetition,".rds")))
  
  #-----------------------------------------------------------------------------------------------------#
  #  ~ TPR and FDR ----
  
  # True positives
  true_diff_cpgs_names <- rownames(coverage)[true_diff_cpgs]
  
  # True negatives
  true_neg_cpg_names <- rownames(coverage)[-c(true_diff_cpgs)]
  
  temp = names(DM_Results_list)
  
  fdr_stats <- lapply(temp, function(DM_Results_method){
    
    # Load DM Results
    DM_Results <- DM_Results_list[[DM_Results_method]]
    
    na_values <- sum(is.na(DM_Results$fdrs))
    
    # Get Significant results
    DM_Results_sign <- DM_Results[DM_Results$fdrs < 0.05, ]
    cpg_ids_sign <- get_cpg_ids(DM_Results_sign)
    
    
    # Get not significant Cpg sites
    DM_Results_notsign <- DM_Results[DM_Results$fdrs >= 0.05, ]
    cpg_ids_notsign <- get_cpg_ids(DM_Results_notsign)
    
    
    # True positives
    true_pos <- length(intersect(cpg_ids_sign, true_diff_cpgs_names))
    
    # False positives
    false_pos <- length(setdiff(cpg_ids_sign, true_diff_cpgs_names))
    
    # True negatives
    true_neg <- length(intersect(cpg_ids_notsign, true_neg_cpg_names))
    
    # False negatives
    false_neg <- length(setdiff(cpg_ids_notsign, true_neg_cpg_names))
    
    # True Positive Rate (TPR)
    TPR <- true_pos / (true_pos + false_neg)
    
    # False discovery rate (FDR)
    FDR <- false_pos / (true_pos+ false_pos)
    
    # Lambda value qqplot
    lambda <-median(qchisq(1-DM_Results$pvals,df= 1))/qchisq(0.5,df = 1)
    
    # Combine into lambda
    stats <- data.frame(true_pos, false_pos, true_neg, false_neg, TPR, FDR, na_values, lambda)
    return(stats)
  })
  
  names(fdr_stats) <- names(DM_Results_list)
  fdr_stats <- do.call("rbind", fdr_stats)
  
  saveRDS(fdr_stats, file.path(s_OUT_dir, "BootstrapSum", paste0("Bootstrap_results_", repetition, ".rds")))
  write.table(fdr_stats, file.path(s_OUT_dir, "BootstrapSum", paste0("Bootstrap_results_", repetition, ".txt")), 
              sep = "\t", quote = FALSE)
  
  rm(DM_Results_list)
  gc()
  
})

rm(betas_subset, covariates, coverage)
gc()


#-----------------------------------------------------------------------------------------------------#
#					8. Compare DM analyses                                                      ----
#-----------------------------------------------------------------------------------------------------#
# 
# # -----------------------------------------------------------------------------------------------------#
# # ~ Combine TPR and FDR results ----
# 
# temp = list.files(file.path(s_OUT_dir, "BootstrapSum/"), pattern="Bootstrap_results",
#                   full.names = TRUE)
# 
# temp = str_subset(temp, ".rds", negate= FALSE)
# 
# bootstrap_results_list <- pblapply(temp, function(file){
#   bootstrap_result <- readRDS(file)
#   method_category = sub("_methyl", "", rownames(bootstrap_result))
#   bootstrap_result <- cbind(bootstrap_result, method_category)
#   return(bootstrap_result)
# })
# 
# bootstrap_results <- do.call(rbind, bootstrap_results_list)
# 
# # Relevel and rename the method_category variable
# # bootstrap_results$method_category <- factor(bootstrap_results$method_category,
# #                                             levels = c("methylkit", "methylkit_MNcorrection", "Limma_GroupAD", "LimmaW_GroupAD", "DSS_GroupAD"),
# #                                             labels = c("methylKit", "methylKit_MNcor" , "DSS","limma", "limma_weights"))
# 
# bootstrap_results$method_category <- factor(bootstrap_results$method_category,
#                                             levels = c("Limma_GroupAD", "LimmaW_GroupAD", "DSS_GroupAD"),
#                                             labels = c("limma", "limma_weights", "DSS"))
# 

# -----------------------------------------------------------------------------------------------------#
# ~ Combine TPR and FDR results (different effect sizes) ----

effect_sizes <- c("5_perc", "10_perc", "20_perc")

bootstrap_results_list <- pblapply(effect_sizes, function(effect_size){
  
  
  temp = list.files(file.path(s_OUT_dir, "BootstrapSum", effect_size), pattern="Bootstrap_results",
                    full.names = TRUE)
  
  temp = str_subset(temp, ".rds", negate= FALSE)
  
  bootstrap_results_list <- pblapply(temp, function(file){
    bootstrap_result <- readRDS(file)
    method_category = sub("_methyl", "", rownames(bootstrap_result))
    bootstrap_result <- cbind(bootstrap_result, method_category)
    return(bootstrap_result)
  })
  
  # return(bootstrap_results_list)
  bootstrap_results <- do.call(rbind, bootstrap_results_list)
  
  # Relevel and rename the method_category variable
  # bootstrap_results$method_category <- factor(bootstrap_results$method_category, 
  #                                             levels = c("methylkit", "methylkit_MNcorrection", "Limma_GroupAD", "LimmaW_GroupAD", "DSS_GroupAD"),
  #                                             labels = c("methylKit", "methylKit_MNcor" , "DSS","limma", "limma_weights"))
  
  bootstrap_results$method_category <- factor(bootstrap_results$method_category, 
                                              levels = c("methylkit", "methylkit_MNcorrection","DSS_GroupAD", "Limma_GroupAD", "LimmaW_GroupAD"),
                                              labels = c("methylKit", "methylKit (MN)" ,  "DSS", "limma", "limma (WLS)"))
  
  return(bootstrap_results)
})


names(bootstrap_results_list) <- c("5perc", "10perc", "20perc")

# Create data frame with gene and contrast
bootstrap_results_list_comb <- bind_rows(
  lapply(names(bootstrap_results_list), function(effect_size) {
    tibble(bootstrap_results_list[[effect_size]], ES = effect_size)
  }))



# -----------------------------------------------------------------------------------------------------#
# ~ Plot TPR results ----
my_colors <- c("#084594","#2171B5","#6BAED6", "#C6DBEF")
my_colors <- c("#084594","#2171B5","#6BAED6", "#718FAB", "#343434", "#C6DBEF")

my_shapes <- c(16, 15, 17, 18, 25)

# Get statistics
get_stats_txt(bootstrap_results_list_comb, "TPR", "method_category", "ES")
get_stats_txt(bootstrap_results_list_comb, "FDR", "method_category", "ES")

my_sum <- bootstrap_results_list_comb %>% 
  group_by(method_category, ES) %>% 
  summarise(TPR_mean = mean(TPR),
            FDR_mean = mean(FDR),
            TPR_median = median(TPR),
            FDR_median = median(FDR),
            TPR_sd = sd(TPR))


barplot_pheno_nosig <- function(data, x_var, y_var, group_var, label_text, xlab_text = "", ylab_text) {
  
  
  destination <- file.path(s_OUT_dir, paste0("Plots/Pheno_barplot_",x_var,"_",y_var,".svg"))
  svglite(file=destination)
  
      # Calculates mean, sd, se and IC
      my_sum <- data %>%
        mutate(ES = fct_relevel(ES, "5perc", "10perc", "20perc")) %>%
        group_by(.data[[x_var]], .data[[group_var]]) %>%
        summarise( 
          N_set = length(.data[[y_var]]),
          mean = round(mean(.data[[y_var]], na.rm = TRUE), digits = 3),
          sd = round(sd(.data[[y_var]], na.rm = TRUE), digits = 3),
          se = round(sd(.data[[y_var]], na.rm = TRUE)/sqrt(length((.data[[y_var]]))), digits= 3),
          .groups = "drop"
        )
      
      print(my_sum)
      
    image <- data %>%
      mutate(ES = fct_relevel(ES, "5perc", "10perc", "20perc")) %>%
      ggplot() +
      geom_bar(data = my_sum, aes(x = .data[[x_var]], y = mean, color = .data[[group_var]]), fill = "white",
               stat= "identity",position=position_dodge(width = 0.8), width = 0.6, linewidth = 0.5) + 
      # geom_jitter(width=0.2, size = 2)+
      geom_point(data = data, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[group_var]], shape = .data[[group_var]]),
                 position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +
      theme_gray()+
      theme(legend.position = "right") +
      xlab(xlab_text) + 
      ylab(ylab_text) +  
      geom_hline(yintercept=1, linetype = "dashed", color = alpha("black", 0.5)) +
      # geom_hline(yintercept=0.05, linetype = "dashed", color = alpha("red", 0.5)) +
      
      geom_errorbar(data = my_sum, aes(x = .data[[x_var]], ymin=mean-se, ymax=mean+se, group = .data[[group_var]], color = .data[[group_var]]), 
                    width=0.2, linewidth = 1.0, alpha=1.0, position = position_dodge((width = 0.8))) +
      # add_pvalue(df_p_val, y.position = "y.position", colour = "black",
      #            tip.length= 0, bracket.nudge.y = 2,
      #            xmin = "xmin",
      #            xmax = "xmax") +
      
      # add_pvalue(df_p_val2, y.position = "y.position", colour = "black",
                 # xmin = "xmin", xmax = "xmax", tip.length = 0, bracket.nudge.y = 2) +
      scale_fill_manual(name = "Methods", values = my_colors, labels = label_text) +
      scale_shape_manual(name = "Methods", values = my_shapes, labels = label_text) +
      scale_color_manual(name = "Methods", values = my_colors, labels = label_text)
    
    
  
  
  print(image)
  dev.off()
  
  cat("Saved plot to: ", destination, "\n")
  
  print(image)
  return(image)
}



violinplot_pheno_nosig <- function(data, x_var, y_var, group_var, label_text, xlab_text = "", ylab_text) {
  
  
  destination <- file.path(s_OUT_dir, paste0("Plots/Pheno_violinplot_",x_var,"_",y_var,".svg"))
  svglite(file=destination)
  
  # Calculates mean, sd, se and IC
  my_sum <- data %>%
    mutate(ES = fct_relevel(ES, "5perc", "10perc", "20perc")) %>%
    group_by(.data[[x_var]], .data[[group_var]]) %>%
    summarise( 
      N_set = length(.data[[y_var]]),
      mean = round(mean(.data[[y_var]], na.rm = TRUE), digits = 3),
      sd = round(sd(.data[[y_var]], na.rm = TRUE), digits = 3),
      se = round(sd(.data[[y_var]], na.rm = TRUE)/sqrt(length((.data[[y_var]]))), digits= 3),
      .groups = "drop"
    )
  
  print(my_sum)
  
  image <- data %>%
    mutate(ES = fct_relevel(ES, "5perc", "10perc", "20perc")) %>%
    ggplot() +
    geom_violin(data = data,  aes(x = .data[[x_var]],  y = .data[[y_var]], color = .data[[group_var]]), fill = "white",
               position=position_dodge(width = 0.8), width = 1.5, linewidth = 0.5)+
  
     # geom_jitter(width=0.2, size = 2)+
    geom_point(data = data, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[group_var]], shape = .data[[group_var]]),
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +
    theme_gray()+
    theme(legend.position = "right") +
    xlab(xlab_text) + 
    ylab(ylab_text) +  
    geom_hline(yintercept=1, linetype = "dashed", color = alpha("black", 0.5)) +
    # geom_hline(yintercept=0.05, linetype = "dashed", color = alpha("red", 0.5)) +
    
    geom_errorbar(data = my_sum, aes(x = .data[[x_var]], ymin=mean-se, ymax=mean+se, group = .data[[group_var]], color = .data[[group_var]]), 
                  width=0.2, linewidth = 1.0, alpha=1.0, position = position_dodge((width = 0.8))) +
    # add_pvalue(df_p_val, y.position = "y.position", colour = "black",
    #            tip.length= 0, bracket.nudge.y = 2,
    #            xmin = "xmin",
    #            xmax = "xmax") +
    
    # add_pvalue(df_p_val2, y.position = "y.position", colour = "black",
    # xmin = "xmin", xmax = "xmax", tip.length = 0, bracket.nudge.y = 2) +
    scale_fill_manual(name = "Methods", values = my_colors, labels = label_text) +
    scale_shape_manual(name = "Methods", values = my_shapes, labels = label_text) +
    scale_color_manual(name = "Methods", values = my_colors, labels = label_text)
  
  
  
  
  print(image)
  dev.off()
  
  cat("Saved plot to: ", destination, "\n")
  
  print(image)
  return(image)
}


# Ridge plot
bootstrap_results_list_comb %>%
  ggplot( aes(y=method_category, x=TPR,  fill=method_category)) +
  # geom_density_ridges(alpha=0.6, stat="binline", bins=20) +
  geom_density_ridges(alpha=0.6) +
  
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("") +
  ylab("Assigned Probability (%)") + facet_wrap("~ES", nrow=3) 



# -----------------------------------------------------------------------------------------------------#
# ~ Plot TPR results ----

# df_p_val2 <- bootstrap_results_list_comb %>%
#   mutate(ES = fct_relevel(ES, "5perc", "10perc", "20perc")) %>%
#   rstatix::group_by(ES) %>%
#   rstatix::t_test(as.formula(paste0("TPR", " ~ ", "method_category"))) %>%
#   rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
#   rstatix::add_significance(p.col = "p.adj") %>%
#   rstatix::add_xy_position(x = "ES", dodge = 0.8)
# 
# print(df_p_val2)
# 
# names(df_p_val2)[3] <- "method_category"

methods <- levels(bootstrap_results_list_comb$method_category)


# Plot statistics
p_TPR_method <- bootstrap_results_list_comb %>%
  mutate(ES = fct_relevel(ES, "5perc", "10perc", "20perc")) %>%
  ggplot(aes(x = ES, y = TPR, colour = method_category, shape = method_category)) +
  geom_boxplot(fill = "white") + 
  # geom_jitter(width=0.2, size = 2)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +
  theme_gray()+
  geom_hline(yintercept=1, linetype = "dashed", color = alpha("black", 0.5)) +
  # theme(axis.line = element_line(color='black'),
  #       plot.background = element_blank(),
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank()) +
  theme(legend.position = "right") +
  xlab("Method") + 
  ylab("True Positive Rate (TPR)") +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  # add_pvalue(df_p_val2, y.position = "y.position", colour = "black",
  #            xmin = "xmin", xmax = "xmax", tip.length = 0) +
  scale_color_manual(name = "Methods", values = my_colors, labels = methods) +
  scale_shape_manual(name = "Methods", values = my_shapes, labels = methods) 
  


p_TPR_method



p_TPR_method_nobar <- barplot_pheno_nosig(bootstrap_results_list_comb, 
                                  x_var = "ES",
                                  y_var="TPR", 
                                  group_var = "method_category", 
                                  label_text = levels(bootstrap_results_list_comb$method_category), 
                                  xlab_text = "",
                                  ylab_text = "True Positive Rate (TPR)")

p_TPR_method_violin <- violinplot_pheno_nosig(bootstrap_results_list_comb, 
                                          x_var = "ES",
                                          y_var="TPR", 
                                          group_var = "method_category", 
                                          label_text = levels(bootstrap_results_list_comb$method_category), 
                                          xlab_text = "",
                                          ylab_text = "True Positive Rate (TPR)")


# -----------------------------------------------------------------------------------------------------#
# ~ Plot FDR results ----

df_p_val2 <- bootstrap_results_list_comb %>%
  mutate(ES = fct_relevel(ES, "5perc", "10perc", "20perc")) %>%
  rstatix::group_by(ES) %>%
  rstatix::t_test(as.formula(paste0("FDR", " ~ ", "method_category"))) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
  rstatix::add_significance(p.col = "p.adj") %>%
  rstatix::add_xy_position(x = "ES", dodge = 0.8)

print(df_p_val2)

names(df_p_val2)[3] <- "method_category"

# Plot statistics
p_FDR_method <- bootstrap_results_list_comb %>%
  mutate(ES = fct_relevel(ES, "5perc", "10perc", "20perc")) %>%
  ggplot(aes(x = ES, y = FDR, colour = method_category, shape = method_category)) +
  geom_boxplot(fill = "white") + 
  # geom_jitter(width=0.2, size = 2)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +
  theme_gray()+
  geom_hline(yintercept=1, linetype = "dashed", color = alpha("black", 0.5)) +
  geom_hline(yintercept=0.05, linetype = "dashed", color = alpha("red", 0.5)) +
  # theme(axis.line = element_line(color='black'),
  #       plot.background = element_blank(),
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank()) +
  theme(legend.position = "right") +
  xlab("Method") + 
  ylab("False Discovery Rate (FDR)") +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  # add_pvalue(df_p_val2, y.position = "y.position", colour = "black",
  #            xmin = "xmin", xmax = "xmax", tip.length = 0) +
  scale_color_manual(name = "Methods", values = my_colors, labels = methods) +
  scale_shape_manual(name = "Methods", values = my_shapes, labels = methods) 


p_FDR_method

p_FDR_method_bar <- barplot_pheno_nosig(bootstrap_results_list_comb, 
                                        x_var = "ES",
                                        y_var="FDR", 
                                        group_var = "method_category", 
                                        label_text = levels(bootstrap_results_list_comb$method_category), 
                                        xlab_text = "",
                                        ylab_text = "False Discovery Rate (FDR)")



# Combine the plots 
destination <- file.path(s_OUT_dir, paste0("Plots/FDR_TPR_Method6.svg"))
svglite(file=destination, width = 15, height = 6)

ggarrange(p_TPR_method_bar, p_FDR_method_bar,
          labels = c("H", "I"),
          ncol = 2, nrow=1)
dev.off()





# -----------------------------------------------------------------------------------------------------#
#  ~ QQ plots ----

effect_sizes <- c("5_perc", "10_perc", "20_perc")

# effect_sizes <- "0_perc"
p_qqplots <- pblapply(effect_sizes, function(effect_size){
  
  temp = list.files(file.path(s_OUT_dir, "BootstrapSum", effect_size), pattern="DM_Results_list_bootstrap",
                    full.names = TRUE)
  

  # QQ plots
  images <- lapply(temp[1], function(file) {
    
    DM_Results_list <- readRDS(file)
    
    bootstrap_number <- sub(".rds", "", sub("DM_Results_list_", "", basename(file)))
    
    names(DM_Results_list) <- c("DSS","limma", "limma (WLS)", "methylKit", "methylKit (MN)")
    
    methods <- names(DM_Results_list)[c(4,5,1,2,3)]
    
    qqplots <- lapply(methods, function(DM_Results_method){
      
      # Load DM Results
      DM_Results <- DM_Results_list[[DM_Results_method]]
      
      # DM_Results_method_name <- sub("_methyl", "", sub("_GroupAD", "", DM_Results_method))
      DM_Results_method_name <- DM_Results_method
      
      
      pvals <- DM_Results$pvals
      
      # Calculate lambda value
      # alpha <- calc_lambda(pvals)
      
      # Plot QQ-plot
      image <- gg_qqplot(pvals) +
        theme_bw() +
        ggtitle(DM_Results_method_name) +
        annotate(
          geom = "text",
          x = -Inf,
          y = Inf,
          hjust = -0.15,
          vjust = 1 + 0.15 * 3,
          label = sprintf("λ = %.2f", inflation(pvals)) # ,size = 8
        ) +
        theme(
          # axis.ticks = element_line(size = 0.5),
          panel.grid = element_blank()
          # panel.grid = element_line(size = 0.5, color = "grey80")
        )
      
      # destination <- file.path(s_OUT_dir, "BootstrapSum", "qqplots", paste0("QQplots_",bootstrap_number, "_", DM_Results_method_name,"_", effect_size,".png", sep=""))
      # png(destination)
      # 
      # print(image)
      # dev.off()
      
      return(image)
    })
    
    return(qqplots)
  
  })
  return(images)
})


# svglite(destination)

combined <- ggarrange(plotlist = unlist(p_qqplots), nrow = 3, ncol = 5)

# destination <- file.path(s_OUT_dir, "BootstrapSum", "qqplots", paste0("QQplots_",bootstrap_number, "_combined.svg"))
# ggsave(destination, combined, device = "svg", width = 15, height = 3.5)

destination <- file.path(s_OUT_dir, "BootstrapSum", "qqplots", paste0("QQplots_combined_allperc_200.png"))
ggsave(destination, combined, device = "png", width = 15, height = 10)

destination <- file.path(s_OUT_dir, "BootstrapSum", "qqplots", paste0("QQplots_combined_allperc_200.svg"))
svglite(destination, width = 15, height = 10)

print(combined)
dev.off()



# temp = list.files(file.path(s_OUT_dir, "BootstrapSum/10_perc"), pattern="DM_Results_list_bootstrap",
#                   full.names = TRUE)


# QQ plots
# pblapply(temp[10], function(file) {
#   
#   DM_Results_list <- readRDS(file)
#   
#   bootstrap_number <- sub(".rds", "", sub("DM_Results_list_", "", basename(file)))
#   
#   methods <- names(DM_Results_list)
#   
#   methods <- methods[c(4,5,1,2,3)]
#   
#   qqplots <- pblapply(methods, function(DM_Results_method){
#     
#     # Load DM Results
#     DM_Results <- DM_Results_list[[DM_Results_method]]
#     
#     DM_Results_method_name <- sub("_methyl", "", sub("_GroupAD", "", DM_Results_method))
#     
#     pvals <- DM_Results$pvals
#     
#     # Calculate lambda value
#     # alpha <- calc_lambda(pvals)
#     
#     # Plot QQ-plot
#     image <- gg_qqplot(pvals) +
#       theme_bw() +
#       ggtitle(DM_Results_method_name) +
#       annotate(
#         geom = "text",
#         x = -Inf,
#         y = Inf,
#         hjust = -0.15,
#         vjust = 1 + 0.15 * 3,
#         label = sprintf("λ = %.2f", inflation(pvals)) # ,size = 8
#       ) +
#       theme(
#         # axis.ticks = element_line(size = 0.5),
#         panel.grid = element_blank()
#         # panel.grid = element_line(size = 0.5, color = "grey80")
#       )
#     
#     # destination <- file.path(s_OUT_dir, "BootstrapSum", "qqplots", paste0("QQplots_",bootstrap_number, "_", DM_Results_method_name,".png", sep=""))
#     # png(destination)
#     # 
#     # print(image)
#     # dev.off()
#     
#     return(image)
#   })
#   
#   
#   # svglite(destination)
#   
#   combined <- ggarrange(plotlist = qqplots, labels = c("B", "C", "D", "E", "F"), nrow = 1)
#   
#   destination <- file.path(s_OUT_dir, "BootstrapSum", "qqplots", paste0("QQplots_",bootstrap_number, "_combined.svg"))
#   ggsave(destination, combined, device = "svg", width = 15, height = 3.5)
#   
#   destination <- file.path(s_OUT_dir, "BootstrapSum", "qqplots", paste0("QQplots_",bootstrap_number, "_combined.png"))
#   ggsave(destination, combined, device = "png", width = 15, height = 3.5)
# })


# -----------------------------------------------------------------------------------------------------#
#  ~ Calculate Time and memory allocation ----

s_OUT_dir_time <- "G:/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/Results/RESEARCH_GSM0172RRMS_ALL/13062025_analysis_Benchmarking/"

# Get list of profiles
temp = list.files(file.path(s_OUT_dir_time, "BootstrapSum/first"), pattern="profile",
                  full.names = TRUE)

# N = 10

temp = str_subset(temp, "profile.html|profile_files", negate = TRUE)
profile_names <- sub(".rds", "", basename(temp))


# Interval is every 0.01 sec = 10 msec

# Calculate time and memory allocation based on each profile
profile_list <- pblapply(temp, function(file) {
  
  # read profile
  profile <- readRDS(file)

  # Extract data
  profdat <- profile$x$message$prof
  
  # Data associated with the DM analysis
  idx_times <- str_detect(profdat$label, "DM_analysis_", negate = FALSE)
  
  # Calculate duration in msecs
  time_data   <- profdat$time[idx_times]
  duration <- length(time_data) * 10 / 1000 # because of interval - in sec
  
  # calculate incremental memory allocation
  mem_data <- profdat$meminc[profdat$time %in% time_data & profdat$meminc >0 ]
  meminc <- round(sum(mem_data), 4) / 1000 # in gigabytes
  
  # Return results
  profile_results <- cbind(duration, meminc)
  return(profile_results)
})

# combine results
profile_stats <- do.call("rbind", profile_list)

method_category <- sub("_profile.*", "", basename(temp))
profile_stats <- data.frame(profile_stats, method_category)

profile_stats$method_category <- factor(profile_stats$method_category, 
                                            levels = c("methylKit", "methylKit_MNcor","DSS", "limma", "limma_weights"),
                                            labels = c("methylKit", "methylKit (MN)" ,  "DSS", "limma", "limma (WLS)"))



# Plot TIME
# my_colors <- c("#00223b", "#084594","#2171B5","#6BAED6", "#C6DBEF")
my_colors <- c("#084594","#2171B5","#6BAED6", "#718FAB", "#343434", "#C6DBEF")

my_shapes <- c(16, 15, 17, 18, 25)
# my_shapes <- c(16, 15, 17, 18)



# -----------------------------------------------------------------------------------------------------#
#  ~ Plot time ----

# Calculate statistics
stats <- get_stats_txt(profile_stats, "duration", "method_category", NULL)
group1 <- str_split_fixed(stats$contrast, " - ",2 )[,1]
group2 <- str_split_fixed(stats$contrast, " - ",2 )[,2]

df_p_val <- data.frame(group1, group2, "p.adj" = convert_p_to_asterix(stats$p.value))
df_p_val <- df_p_val[df_p_val$p.adj != "NS",]
df_p_val <- cbind(df_p_val, "y.position" = seq(600, by = 50, length.out = nrow(df_p_val)))


# Plot statistics (n =)
p_time_method <- profile_stats %>%
  ggplot(aes(x = method_category, y = duration, colour = method_category)) +
  geom_boxplot(fill = "white") + 
  # geom_jitter(width=0.2, size = 2)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.8) +
  theme_gray()+
  # theme(axis.line = element_line(color='black'),
  #       plot.background = element_blank(),
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank()) +
  theme(legend.position = "right") +
  scale_color_manual(name = "Methods", values = my_colors) +
  scale_shape_manual(name = "Methods", values = my_shapes) +
  # add_pvalue(df_p_val, y.position = "y.position", colour = "black", tip.length= 0) 
  xlab("Method") + 
  ylab("Duration (sec)")

p_time_method


# -----------------------------------------------------------------------------------------------------#
#  ~ Plot memory allocation ----

# Calculate statistics
stats <- get_stats_txt(profile_stats, "meminc", "method_category", NULL)
group1 <- str_split_fixed(stats$contrast, " - ",2 )[,1]
group2 <- str_split_fixed(stats$contrast, " - ",2 )[,2]

df_p_val <- data.frame(group1, group2, "p.adj" = convert_p_to_asterix(stats$p.value))
df_p_val <- df_p_val[df_p_val$p.adj != "NS",]
df_p_val <- cbind(df_p_val, "y.position" = seq(250, by = 10, length.out = nrow(df_p_val)))

p_meminc_method <- profile_stats %>%
  ggplot(aes(x = method_category, y = meminc, colour = method_category)) +
  geom_boxplot(fill = "white") + 
  # geom_jitter(width=0.2, size = 2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.8) +
  theme_gray()+
  # theme(axis.line = element_line(color='black'),
  #       plot.background = element_blank(),
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank()) +
  theme(legend.position = "right") +
  scale_color_manual(name = "Methods", values = my_colors) +
  scale_shape_manual(name = "Methods", values = my_shapes) +
  # add_pvalue(df_p_val, y.position = "y.position", colour = "black", tip.length= 0) +
  xlab("Method") + 
  ylab("Incremental memory allocation (Gb)")

p_meminc_method

destination <- file.path(s_OUT_dir, paste0("Plots/Memory_Time_Method3.svg"))
svglite(file=destination, width = 15, height = 6)

ggarrange(p_time_method, p_meminc_method,
          labels = c("F", "G"),
          ncol = 2, nrow=1)
dev.off()


# -----------------------------------------------------------------------------------------------------#
#  ~ Combine all plots ----

destination <- file.path(s_OUT_dir, paste0("Plots/Memory_Time_FDR_TPR_Method4.svg"))
svglite(file=destination, width = 19, height = 5)

ggarrange(p_time_method + theme(legend.position="none"), 
          p_meminc_method + theme(legend.position="none"), 
          p_TPR_method_bar + theme(legend.position="none"), p_FDR_method_bar + theme(legend.position="none"),
          labels = c("F", "G", "H", "I"),
          ncol = 4, nrow=1)

dev.off()



destination <- file.path(s_OUT_dir, paste0("Plots/Memory_Time_FDR_TPR_Method_vertical2.svg"))
svglite(file=destination, width = 7, height = 13)
# 
# ggarrange(p_time_method + theme(legend.position="none"), 
#           p_meminc_method + theme(legend.position="none"), 
#           p_TPR_method + theme(legend.position="none"), p_FDR_method + theme(legend.position="none"),
#           labels = c("C", "D", "E", "F"),
#           ncol = 1, nrow=4)


ggarrange(p_time_method + theme(legend.position="right"), 
          p_meminc_method + theme(legend.position="right"), 
          p_TPR_method_bar + theme(legend.position="right"), p_FDR_method_bar + theme(legend.position="right"),
          labels = c("C", "D", "E", "F"),
          ncol = 1, nrow=4)
dev.off()



# 
# # -----------------------------------------------------------------------------------------------------#
# # ~ Plot TPR results ----
# my_colors <- c("#084594","#2171B5","#6BAED6", "#C6DBEF")
# my_shapes <- c(16, 15, 17, 18)
# 
# # Get statistics
# get_stats_txt(bootstrap_results, "TPR", "method_category", NULL)
# get_stats_txt(bootstrap_results, "FDR", "method_category", NULL)
# 
# bootstrap_results %>% 
#   group_by(method_category) %>% 
#   summarise(TPR_mean = mean(TPR),
#             FDR_mean = mean(FDR),
#             TPR_median = median(TPR),
#             FDR_median = median(FDR),
#             TPR_sd = sd(TPR))
# 
# 
# 
# # Calculate statistics for p values
# stats <- get_stats_txt(bootstrap_results, "TPR", "method_category", NULL)
# group1 <- str_split_fixed(stats$contrast, " - ",2 )[,1]
# group2 <- str_split_fixed(stats$contrast, " - ",2 )[,2]
# 
# df_p_val <- data.frame(group1, group2, "p.adj" = convert_p_to_asterix(stats$p.value))
# df_p_val <- df_p_val[df_p_val$p.adj != "NS",]
# df_p_val <- cbind(df_p_val, "y.position" = seq(1.1, by = 0.1, length.out = nrow(df_p_val)))
# 
# # Plot statistics
# p_TPR_method <- bootstrap_results %>%
#   ggplot(aes(x = method_category, y = TPR, colour = method_category)) +
#   geom_boxplot(fill = "white") + 
#   # geom_jitter(width=0.2, size = 2)+
#   geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.8) +
#   theme_gray()+
#   geom_hline(yintercept=1, linetype = "dashed") +
#   # theme(axis.line = element_line(color='black'),
#   #       plot.background = element_blank(),
#   #       panel.grid.major = element_blank(),
#   #       panel.grid.minor = element_blank()) +
#   theme(legend.position = "right") +
#   scale_color_manual(name = bootstrap_results$method_category, values = my_colors) +
#   scale_shape_manual(name = bootstrap_results$method_category, values = my_shapes) +
#   xlab("Method") + 
#   ylab("True Positive Rate (TPR)") +
#   add_pvalue(df_p_val, y.position = "y.position", colour = "black", tip.length= 0)
# 
# # -----------------------------------------------------------------------------------------------------#
# # ~ Plot FDR results ----
# 
# # Calculate statistics FDR
# stats <- get_stats_txt(bootstrap_results, "FDR", "method_category", NULL)
# group1 <- str_split_fixed(stats$contrast, " - ",2 )[,1]
# group2 <- str_split_fixed(stats$contrast, " - ",2 )[,2]
# 
# df_p_val <- data.frame(group1, group2, "p.adj" = convert_p_to_asterix(stats$p.value))
# df_p_val <- df_p_val[df_p_val$p.adj != "NS",]
# df_p_val <- cbind(df_p_val, "y.position" = seq(0.2, by = 0.1, length.out = nrow(df_p_val)))
# 
# # Plot statistics
# p_FDR_method <- bootstrap_results %>%
#   ggplot(aes(x = method_category, y = FDR, colour = method_category)) +
#   geom_boxplot(fill = "white") + 
#   # geom_jitter(width=0.2, size = 2)+
#   geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.8) +
#   theme_gray()+
#   # geom_hline(yintercept=1, linetype = "dashed") +
#   # theme(axis.line = element_line(color='black'),
#   #       plot.background = element_blank(),
#   #       panel.grid.major = element_blank(),
#   #       panel.grid.minor = element_blank()) +
#   theme(legend.position = "right") +
#   scale_color_manual(name = bootstrap_results$method_category, values = my_colors) +
#   scale_shape_manual(name = bootstrap_results$method_category, values = my_shapes) +
#   xlab("Method") + 
#   ylab("False Discovery Rate (TPR)") +
#   add_pvalue(df_p_val, y.position = "y.position", colour = "black", tip.length= 0)
# 
# 
# destination <- file.path(s_OUT_dir, paste0("Plots/FDR_TPR_Method2.svg"))
# svglite(file=destination, width = 15, height = 6)
# 
# ggarrange(p_TPR_method, p_FDR_method,
#           labels = c("H", "I"),
#           ncol = 2, nrow=1)
# dev.off()


