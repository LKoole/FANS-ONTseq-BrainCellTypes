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
s_ROOT_dir <<- "path/to/root/directory"
source(paste0(s_ROOT_dir,"Scripts/.Main/Settings_v2.R"))


#-----------------------------------------------------------------------------------------------------#
#							Libraries
#-----------------------------------------------------------------------------------------------------#
# load library
library(readxl)
library(data.table)
library(dplyr)
library(ggplot2)
library(svglite)
library(ggpubr)
library(vtable)

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

anno <- cbind(anno, "Group"= annotation_Cases$DX, "Subject_number" = annotation_Cases$Number) 


#-----------------------------------------------------------------------------------------------------#
#						Change Gender variable (female - male)
#-----------------------------------------------------------------------------------------------------#
# Check coding of Gender in Anno file
table(anno$Gender)

# Convert 1 into female and 2 into male 
if (isTRUE(anno$Gender == 1 | anno$Gender == 2) == FALSE) {
anno$Gender <- ifelse(anno$Gender == 1, "male", "female")} else {
   print("Change the coding for variable Gender to 1 and 2") }


#-----------------------------------------------------------------------------------------------------#
#							Select samples
#-----------------------------------------------------------------------------------------------------#
# e.g. AD and MCI samples
# Samples_of_interest_index = c(anno$Group == "AD" | anno$Group ==  "MCI")

# e.g. all samples
Samples_of_interest_index = c(anno$Group == "AD" | anno$Group ==  "MCI" | anno$Group == "CTL")

# select samples to make it less heavy on RAM.
# Samples_of_interest_index = Samples_of_interest_index & (anno$CaseID %in% 110:114)

# make Pheno object
pheno = anno[Samples_of_interest_index,]

pheno$Age <- as.numeric(pheno$Age)
pheno$Group <- factor(pheno$Group, levels = c("CTL", "MCI", "AD"))

# So there is a special case with subject 16-51. Has different length in CpG probes, which is weird. excluding 
pheno = pheno[pheno $CaseID != "16-51",]

table(pheno$Group)
# CTL MCI  AD 
# 74  50  75 


#-----------------------------------------------------------------------------------------------------#
#							Annotation ONT sequencing run
#-----------------------------------------------------------------------------------------------------#

annotation_ONT <- read_excel(paste0(s_ROOT_dir, "Anno/GSM0172RRMS_overview_pheno.xlsx"), sheet =1)

annotation_ONT$CaseID <- sub("_NEW|_OLD", "", annotation_ONT$CaseID)

annotation_ONT$CaseID <- gsub(annotation_ONT$CaseID, pattern = "_", replace ="-")

# Select flowcell (all)
annotation_ONT <- annotation_ONT[!is.na(annotation_ONT$Coverage),]


# Select cell type 
annotation_ONT <- annotation_ONT[annotation_ONT$Cell_type=="NeuN",]


pheno[match(annotation_ONT$CaseID, pheno$CaseID), ]

pheno <- left_join(pheno, annotation_ONT, by = "CaseID")


# Overview of samples included (all cell types)
table(pheno[, "Group"])


# Select NeuN
table(pheno[pheno$Cell_type == "NeuN" , "Group"])

# CTL MCI  AD 
# 37  34  41

# CTL MCI  AD 
# 54  41  54 

# CTL MCI  AD 
# 69  45  66 


#-----------------------------------------------------------------------------------------------------#
#							Plot metadata
#-----------------------------------------------------------------------------------------------------#

# Gender
gender_plot <- ggplot(pheno, aes(x = Gender, fill = Group)) +
  geom_bar(position = "dodge", width=0.8) + 
  scale_fill_brewer(palette="Blues", direction = -1) +  
  theme_bw()

# ApoE
ApoE_plot <- ggplot(pheno, aes(x = ApoE, fill = Group)) +
  geom_bar(position = "dodge", width=0.8) +
  scale_fill_brewer(palette="Blues", direction = -1) +  theme_bw()

# Braak score
braak_plot <- ggplot(pheno, aes(x = Braak_score, fill = Group)) +
  geom_bar(position = "dodge", width=0.8) +
  scale_fill_brewer(palette="Blues", direction = -1) +  theme_bw() +
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6"))


# MMSE scores
mmse_plot <- ggplot(pheno, aes(x = MMSE_score, fill = Group)) +
  geom_histogram(position = "identity", alpha = 0.9, bins = 10) +
  facet_wrap(~ Group) +
  scale_fill_brewer(palette="Blues", direction = -1) +
  theme_bw()

# Age
age_plot <- ggplot(pheno, aes(x = Age, fill = Group)) +
  geom_histogram(position = "identity", alpha = 1.0, bins = 20) +
  facet_wrap(~ Group) +
  scale_fill_brewer(palette="Blues", direction = -1) +
  theme_bw()

# DNA input
dna_plot <- ggplot(pheno, aes(x = DNA_input, fill = Group)) +
  geom_histogram(position = "identity", alpha = 1.0, bins = 10) +
  facet_wrap(~ Group) +
  scale_fill_brewer(palette="Blues", direction = -1) +
  theme_bw()

destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/Pheno_phenotypeoverview.svg", sep="")

svglite(destination, width = 15, height = 10)

ggarrange(gender_plot, ApoE_plot, braak_plot, mmse_plot, age_plot,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow=2)

dev.off()

rm(gender_plot, ApoE_plot, braak_plot, mmse_plot, age_plot, dna_plot)

gc()


st(pheno, group = "Group", add.median = TRUE, group.test=TRUE)

st(pheno, group = "Group", add.median = TRUE, group.test=list(format = paste0("{name}={stat}", " p = ", "{pval}{stars}")))

st(pheno, group = "Group", add.median = TRUE, group.test=list('test'= vtable:::groupt.it))


chisq.test(pheno$Group, pheno$Gender)

chisq.test(pheno$Group, pheno$Age)

get_stats_txt(pheno,"Age",  "Group", NULL)

pheno_no_na <- pheno[!is.na(pheno$Barcode),]

st(pheno_no_na, group = "Group", add.median = TRUE, group.test=TRUE)


#-----------------------------------------------------------------------------------------------------#
#							Save anno
#-----------------------------------------------------------------------------------------------------# 
# Save the Pheno object in the results folder.
save(pheno,file = file.path(paste0(s_ROOT_dir,s_out_folder,"Pheno/Pheno.Rdata"))) 

