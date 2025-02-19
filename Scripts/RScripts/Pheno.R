#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#		Pheno.R
#
#	Purpose 
#		This code is for Pheno prep of DNA methylation and hydroxymethylation data
#
# Author comment:
#	lisa.koole@maastrichtuniversity.nl
#
#	Personal comment:
#
#
#-----------------------------------------------------------------------------------------------------#
#							Main settings or load
#-----------------------------------------------------------------------------------------------------#
# Get all general settings 
s_ROOT_dir <<- "G:/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/"
source(paste0(s_ROOT_dir,"Scripts/.Main/Settings.R"))


#-----------------------------------------------------------------------------------------------------#
#							Libraries
#-----------------------------------------------------------------------------------------------------#
# load library
library("readxl")
library("data.table")
library("dplyr")


#-----------------------------------------------------------------------------------------------------#
#							Load metadata about cases
#-----------------------------------------------------------------------------------------------------#

annotation_Cases = read_excel(paste0(s_ROOT_dir, "Anno/BReIN lijst finaal.xlsx"), sheet=1)

# Change CN (control) to CTL
annotation_Cases$DX[annotation_Cases$DX == "CN"] <- "CTL"

# Select Cohort, Sample number, and CaseID 
annotation_Cases = annotation_Cases[,1:3]

# Remove rows lacking sample information
annotation_Cases$Number <- as.numeric(annotation_Cases$Number)
annotation_Cases <- annotation_Cases[is.na(annotation_Cases$Number) == FALSE,]


#-----------------------------------------------------------------------------------------------------#
#							Load metadata about details per case (raw annotation)
#-----------------------------------------------------------------------------------------------------#

annotation_raw = read_excel(paste0(s_ROOT_dir, "Anno/BReIN lijst finaal.xlsx"), sheet="details per case")

# get the columns of interest
# 1 is CaseID
# 2 is race
# 3 is gender
# 4 is expired age
# 6 is PMI (post mortem interval)
# 7 is last mmse test score
# 9 is Control (yes or no)
# 10 is AD (yes or no)
# 27 is MCI (yes or no)
# 31 is ApoE
# 32 is Plaque density
# 43 is Braak score

anno = data.frame(
  "CaseID" = annotation_raw[,1],
  "Gender" = annotation_raw[,3],
  "Age" = suppressWarnings(as.numeric(unlist(annotation_raw[,4]))),
  "PMI" = suppressWarnings(as.numeric(unlist(annotation_raw[,6]))),
  "MMSE_score" = suppressWarnings(as.numeric(unlist(annotation_raw[,7]))),
  # "Control" = Annotation_raw[,9],
  # "AD" = Annotation_raw[,10],
  #  "MCI" = Annotation_raw[,27],
  "ApoE" = annotation_raw[,31],
  "Plaque_density" = annotation_raw[,32],
  "Braak_score" = suppressWarnings(as.numeric(unlist(annotation_raw[,34])))
)  

colnames(anno) <- c("CaseID","Gender","Age", "PMI", "MMSE_score", "ApoE",
                    "Plaque_density", "Braak_score")

# Remove rows lacking information about cases (NA values or headers)
anno <- anno[-(which(is.na(anno$CaseID) | anno$CaseID == "CaseID")),]

# Add annotation_cases information to the Anno dataframe  
annotation_Cases <- annotation_Cases[match(anno$CaseID, annotation_Cases$case),]

anno <- cbind(anno, "Cohort"= annotation_Cases$DX, "Subject_number" = annotation_Cases$Number) 


#-----------------------------------------------------------------------------------------------------#
#						Change Gender variable (female - male)
#-----------------------------------------------------------------------------------------------------#
# Check coding of Gender in Anno file
table(anno$Gender)

# Convert 1 into female and 2 into male 
if (isTRUE(anno$Gender == 1 | anno$Gender == 2) == FALSE) {
anno$Gender <- ifelse(anno$Gender == 1, "female", "male")} else {
   print("Change the coding for variable Gender to 1 and 2") }


#-----------------------------------------------------------------------------------------------------#
#							Select samples
#-----------------------------------------------------------------------------------------------------#
# e.g. AD and MCI samples
# Samples_of_interest_index = c(anno$Cohort == "AD" | anno$Cohort ==  "MCI")

# e.g. all samples
Samples_of_interest_index = c(anno$Cohort == "AD" | anno$Cohort ==  "MCI" | anno$Cohort == "CTL")

# select samples to make it less heavy on RAM.
# Samples_of_interest_index = Samples_of_interest_index & (anno$CaseID %in% 110:114)

# make Pheno object
pheno = anno[Samples_of_interest_index,]

pheno$Age <- as.numeric(pheno$Age)
pheno$Cohort <- factor(pheno$Cohort, levels = c("CTL", "MCI", "AD"))

# So there is a special case with subject 114. Has different length in CpG probes, which is weird. excluding 
pheno = pheno[pheno $CaseID != "16-51",]

table(pheno$Cohort)
# AD CTL MCI 
# 75  74  50 

# Overview of samples included
table(pheno[!is.na(pheno$Flowcell), "Cohort"])
# CTL MCI  AD 
# 23  21  22 

#-----------------------------------------------------------------------------------------------------#
#							Annotation ONT sequencing run
#-----------------------------------------------------------------------------------------------------#

annotation_ONT <- read_excel(paste0(s_ROOT_dir, "Anno/GSM1072RRMS_overview_pheno.xlsx"))

annotation_ONT$CaseID <- gsub(annotation_ONT$CaseID, pattern = "_", replace ="-")

pheno[match(annotation_ONT$CaseID, pheno$CaseID), ]

pheno <- left_join(pheno, annotation_ONT, by = "CaseID")

#-----------------------------------------------------------------------------------------------------#
#							Save anno
#-----------------------------------------------------------------------------------------------------# 
# Save the Pheno object in the results folder.
save(pheno,file = paste0(s_ROOT_dir,s_out_folder,"Pheno/Pheno.Rdata")) 

