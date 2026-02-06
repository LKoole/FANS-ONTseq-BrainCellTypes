#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#		MatchAtlas_v2.R
#
#	Purpose 
#		This code was made for comparing the human methylation atlas (Loyfer et al., 2023) and NeuN+ fraction
#
# Author comment:
#
#
#-----------------------------------------------------------------------------------------------------#
#							Main settings or load                                             ---- 
#-----------------------------------------------------------------------------------------------------#
# Set root folder
s_ROOT_dir <<- "path/to/root/directory"

# Get all general settings 
source(paste0(s_ROOT_dir,"Scripts/.Main/Settings_v2.R"))

#-----------------------------------------------------------------------------------------------------#
#							Libraries needed                                                  ----
#-----------------------------------------------------------------------------------------------------#


library(readxl) # read excel
library(data.table)
library(dplyr)
library(stringr)
library(annotatr)
library(GenomicRanges) # create GenomicRanges
library(GenomicDistributionsData) # TSS
library(ggplot2)
library(tidyr)
library(dplyr)
library(svglite)


#-----------------------------------------------------------------------------------------------------#
#							Select control samples                                                 ----
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


sample_names_CTL <- pheno_subset$sample_names[pheno_subset$Group == "CTL"]


#-----------------------------------------------------------------------------------------------------#
#						Get modified (5mC + 5hmC) values                                                  ----
#-----------------------------------------------------------------------------------------------------#

# Load modified data (5mC + 5hmC levels combined)
modified <- readRDS("G:/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/Results/RESEARCH_GSM0172RRMS_ALL/13062025_analysis_Limma/Bedmethyl_all/bedfiles_modified_all_5X.rds")

modified_CTL <- modified[,colnames(modified) %in% sample_names_CTL]




#-----------------------------------------------------------------------------------------------------#
#						Get CpG positions                                                 ----
#-----------------------------------------------------------------------------------------------------#

position <- str_split_fixed(rownames(modified_CTL), "_", n = 2)
cpg_ranges <- GRanges(seqnames = position[,1],
                      IRanges(start = as.numeric(position[,2]),
                              end = as.numeric(position[,2])+1))

cpg_ids <- as.data.table(rownames(modified_CTL))

# Free up memory
rm(position)
gc()


#-----------------------------------------------------------------------------------------------------#
#						Read the human methylation atlast file and create GRanges object         ----
#-----------------------------------------------------------------------------------------------------#

atlas <- read_excel(paste0(s_ROOT_dir, "Anno/celltype_atlas.xlsx"), sheet = 4, skip = 2)

position <- str_split_fixed(atlas$position_hg38, ":", 2)

position <- str_split_fixed(position[,2], "-", 2)


atlas_ranges <- GRanges(seqnames = atlas$chr,
                        IRanges(start = as.numeric(position[,1]),
                                end = as.numeric(position[,2])))

mcols(atlas_ranges) <- data.frame(atlas)


#-----------------------------------------------------------------------------------------------------#
#					Overlap atlas regions and nanopore data     ----
#-----------------------------------------------------------------------------------------------------#
# Find overlapping regions (more hits per cpg site)
overlaps <- findOverlaps(cpg_ranges, atlas_ranges, ignore.strand=TRUE, select = "all")

# Extract annotation data
ann_data <- as.data.table(atlas_ranges@elementMetadata@listData)
ann_data[, ann_ids := .I]  # Keep original row index


valid_cpg_ids <- cpg_ids[queryHits(overlaps)]

# Build overlaps table
overlap_dt <- data.table(
  cpg_id = valid_cpg_ids,
  ann_ids = overlaps@to
)

names(overlap_dt) <- c("cpg_id", "ann_ids")


#-----------------------------------------------------------------------------------------------------#
#				Beta values of atlas     ----
#-----------------------------------------------------------------------------------------------------#
ann_data_merged <- merge(overlap_dt, ann_data, by = "ann_ids", all.x = TRUE)

ann_data_merged <- ann_data_merged[!duplicated(ann_data_merged$cpg_id),]

target_betas <- ann_data_merged$Target.meth. * 100
ann_data_merged <- cbind(ann_data_merged, target_betas)


#-----------------------------------------------------------------------------------------------------#
#				Beta values of nanopore data     ----
#-----------------------------------------------------------------------------------------------------#

modified_CTL_filt <- modified_CTL[rownames(modified_CTL) %in% ann_data_merged$cpg_id,]

mean_betas <- data.frame(apply(modified_CTL_filt, 1, function(x) {mean(x, na.rm = TRUE)} ))
mean_betas <- cbind(mean_betas, rownames(mean_betas))
colnames(mean_betas) = c("mean_betas", "cpg_id")

median_betas <- data.frame(apply(modified_CTL_filt, 1, function(x) {median(x, na.rm = TRUE)} ))
median_betas <- cbind(median_betas, rownames(median_betas))
colnames(median_betas) = c("median_betas", "cpg_id")


#-----------------------------------------------------------------------------------------------------#
#			Merge beta values of nanopore data and atlas data    ----
#-----------------------------------------------------------------------------------------------------#

ann_data_merged <- left_join(ann_data_merged, mean_betas)
ann_data_merged <- left_join(ann_data_merged, median_betas)

# Aggregate beta values belonging to same region (mean value)
my_sum <- ann_data_merged %>%
  group_by(ann_ids) %>%
  summarise( 
    mean_betas_agreg = round(mean(mean_betas, na.rm = TRUE), digits = 2)
  ) 

ann_data_merged <- left_join(ann_data_merged, my_sum)


#-----------------------------------------------------------------------------------------------------#
#			Boxplots beta values per cell type    ----
#-----------------------------------------------------------------------------------------------------#

boxplot(ann_data_merged$mean_betas ~ ann_data_merged$Type, las = 2) # Nanopore data
boxplot(ann_data_merged$target_betas ~ ann_data_merged$Type, las = 2) # Atlas data



#-----------------------------------------------------------------------------------------------------#
#			Heatmap beta values per cell type    ----
#-----------------------------------------------------------------------------------------------------#

data <- reshape2::melt(ann_data_merged, id.vars = c("cpg_id", "Type"), measure.vars = c("target_betas", "mean_betas"))

# Mean beta values; cpg_ids all
ggplot(data, aes(x = variable, y = cpg_id, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") + facet_wrap("~Type") + theme(axis.text.y = element_blank() )

# Median beta value; cpg_ids all
data <- reshape2::melt(ann_data_merged, id.vars = c("cpg_id", "Type"), measure.vars = c("target_betas", "median_betas"))

ggplot(data, aes(x = variable, y = cpg_id, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu")


# Regions specifically for neuron and oligodend
data <- reshape2::melt(ann_data_merged_NeuN, id.vars = c("ann_ids", "Type"), measure.vars = c("target_betas", "mean_betas"))

data$ann_ids <- as.character(data$ann_ids)

ggplot(data, aes(x = variable, y = ann_ids, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") + facet_wrap("~Type") + theme(axis.text.y = element_blank())



#-----------------------------------------------------------------------------------------------------#
#			Heatmap beta values per cell type (final plot)  ----
#-----------------------------------------------------------------------------------------------------#
data <- reshape2::melt(ann_data_merged, id.vars = c("ann_ids", "Type", "Gene"), measure.vars = c("target_betas", "mean_betas_agreg"))
data <- data[!duplicated(data),]

# Relabel beta value sets to Reference (Ref.) vs NeuN+ data
data$variable <- factor(data$variable, levels = c("target_betas", "mean_betas_agreg"), labels = c("Ref.", "NeuN+"))


# Shorten cell type names
data$Type[data$Type == "Skeletal-Musc:Smooth-Musc"] <- "Skeletal/Smooth-Musc"
data$Type[data$Type == "Pancreas-Alpha:Pancreas-Beta:Pancreas-Delta"] <- "Pancreas-Alpha/Beta/Delta"
data$Type[data$Type == "Lung-Ep-Alveo:Lung-Ep-Bron"] <- "Lung-Ep-Alveo:Bron"
data$Type[data$Type == "Heart-Cardio:Skeletal-Musc:Smooth-Musc"] <- "Heart-Car:Skel/Smooth-Musc"
data$Type[data$Type == "Colon-Ep:Gastric-Ep:Small-Int-Ep"] <- "Colon:Gastric:Small-Int-Ep"
data$Type[data$Type == "Breast-Basal-Ep:Breast-Luminal-Ep"] <- "Breast-Basal/Luminal-Ep"





# All cell types including genes overlapping multiple cell types
data %>% 
  arrange(Type) %>%
  mutate(ann_ids = factor(ann_ids, levels = unique(ann_ids))) %>%
  ggplot(aes(x = variable, y = ann_ids, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") + facet_wrap("~Type", axes = "all_x", scales = "free_y") + theme(axis.text.y = element_blank()) +
  theme_bw()



destination = file.path(s_OUT_dir, "Plots", "Atlas_heatmap.svg")

svglite(destination, width = 19, height = 11)

data %>% 
  arrange(Type) %>%
  mutate(ann_ids = factor(ann_ids, levels = unique(ann_ids))) %>%
  ggplot(aes(x = variable, y = Gene, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") + facet_wrap("~Type", axes = "all_x", scales = "free_y") + 
  theme_bw() +
  xlab("") +
  labs(fill = "Mean 5mC+5hmC") +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "top"
  )

dev.off()


# All cell types excluding the combined cell types e.g. Neuron:Oligodend

data_subset <- data[str_detect(data$Type, ":", negate = TRUE),] 

data_subset$Type <- sub("-Ep", "-Epithelium", data_subset$Type)
data_subset$Type <- sub("-Fibro", "-Fibroblast", data_subset$Type)
data_subset$Type <- sub("-Musc", "-Muscle", data_subset$Type)

table(data_subset$Type)


# Extend names

data_subset$Type[data_subset$Type == "Blood-Granul"] <- "Blood-Granulocytes"
data_subset$Type[data_subset$Type == "Blood-Mono+Macro"] <- "Blood-Monocytes/Macrophages"

data_subset$Type[data_subset$Type == "Bone-Osteob"] <- "Bone-Osteoblast"

data_subset$Type[data_subset$Type == "Endothel"] <- "Endothelium"
data_subset$Type[data_subset$Type == "Epid-Kerat"] <- "Epidermal-Keratinocytes"
data_subset$Type[data_subset$Type == "Eryth-prog"] <- "Erythrocyte-Progenitors"
data_subset$Type[data_subset$Type == "Heart-Cardio"] <- "Heart-Cardiocytes"
data_subset$Type[data_subset$Type == "Liver-Hep"] <- "Liver-Hepatocytes"
data_subset$Type[data_subset$Type == "Lung-Epithelium-Bron"] <- "Lung-Bronchial-Epithelium"
data_subset$Type[data_subset$Type == "Lung-Epithelium-Alveo"] <- "Lung-Alveolar-Epithelium"
data_subset$Type[data_subset$Type == "Neuron"] <- "Neurons"

data_subset$Type[data_subset$Type == "Oligodend"] <- "Oligodendrocytes"
data_subset$Type[data_subset$Type == "Small-Int-Epithelium"] <- "Small-Intenstine-Epithelium"



# All cell types including genes overlapping multiple cell types
data_subset %>% 
  arrange(Type) %>%
  mutate(ann_ids = factor(ann_ids, levels = unique(ann_ids))) %>%
  ggplot(aes(x = variable, y = Gene, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") + facet_wrap("~Type", axes = "all_x", scales = "free_y") + theme(axis.text.y = element_blank()) +
  theme_bw()



destination = file.path(s_OUT_dir, "Plots", "Atlas_heatmap2.svg")

svglite(destination, width = 19, height = 11)

data_subset %>% 
  arrange(Type) %>%
  mutate(ann_ids = factor(ann_ids, levels = unique(ann_ids))) %>%
  ggplot(aes(x = variable, y = Gene, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") + facet_wrap("~Type", axes = "all_x", scales = "free_y") + 
  theme_bw() +
  xlab("") +
  labs(fill = "Mean 5mC+5hmC") +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "top"
  )

dev.off()



