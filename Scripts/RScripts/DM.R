#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#	EWAS.R
#
#	Purpose 
#		This code was made for performing EWAS
#
# Author comment:
#	Lisa Koole
#	lisa.koole@maastrichtuniversity.nl
#
#
#-----------------------------------------------------------------------------------------------------#
#							Main settings or load
#-----------------------------------------------------------------------------------------------------#
# Get all general settings 

s_ROOT_dir <<- "G:/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/"
source(paste0(s_ROOT_dir,"Scripts/.Main/Settings.R"))


#-----------------------------------------------------------------------------------------------------#
#							        Libraries needed
#-----------------------------------------------------------------------------------------------------#
library(limma) # DM
library(dplyr) 
library(DMRcate) # DM		  
library(missMethyl) # DM and GO
library(qqman)
library(gplots) 
library(readr)
library(htmltools)
library(stringr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(purrr)
library(readxl)
library(GenomicRanges)
library(EnhancedVolcano)


BiocManager::install("GenomicDistributionsData")
library(GenomeInfoDb)
library(GenomicDistributions)
library(GenomicDistributionsData)




#-----------------------------------------------------------------------------------------------------#
#					            	Load data
#-----------------------------------------------------------------------------------------------------#

load(paste0(s_ROOT_dir,s_out_folder,"QC/methyl_Metas.Rdata")) 
load(paste0(s_ROOT_dir,s_out_folder,"QC/hydroxymethyl_Metas.Rdata"))
load(paste0(s_ROOT_dir, s_out_folder, "QC/unmod_Metas.Rdata"))
load(paste0(s_ROOT_dir, s_out_folder, "QC/scores_filtered.Rdata"))
load(paste0(s_ROOT_dir, s_out_folder, "QC/positions_filtered.Rdata"))
load(paste0(s_ROOT_dir, s_out_folder, "QC/scores_filtered.Rdata"))
load(paste0(s_ROOT_dir, s_out_folder, "Pheno/Pheno.Rdata"))


# Create matrices
methyl_Metas <- as.matrix(methyl_Metas)
hydroxymethyl_Metas <- as.matrix(hydroxymethyl_Metas)
unmod_Metas <- as.matrix(unmod_Metas)


# Match samples ... 
sample_names = sub("sample_", "", colnames(methyl_Metas))
pheno_subset <- pheno[match(sample_names, pheno$CaseID),]


#-----------------------------------------------------------------------------------------------------#
#							        Design Model
#-----------------------------------------------------------------------------------------------------#
# Make design matrix formula; correcting for Subject, determining the effect of Tissue
s_CovFormulaImportantSV = "~ Cohort + Age + Gender"

# make design matrix
temp_design = model.matrix(as.formula(s_CovFormulaImportantSV),data=pheno_subset)
colnames(temp_design) = make.names(colnames(temp_design))


#-----------------------------------------------------------------------------------------------------#
#				1. Differential methylation and hydroxymethylation analysis using limma
#-----------------------------------------------------------------------------------------------------#

DM_analysis <- function(Metas, mod_type) {
  
  # Make LMmodel 
  temp_fit = limma::lmFit(as.matrix(Metas),temp_design)
  
  # Make LMmodel with weights
  #temp_fit = limma::lmFit(as.matrix(methyl_Metas),temp_design, weights=scores_weights)
  
  # Show the coefs for first gene
  data.frame(temp_fit$coef[1,])
  
  temp_contrast_matrix <- limma::makeContrasts(
    MCI_vs_HC = CohortMCI,
    AD_vs_HC = CohortAD,
    AD_vs_MCI = CohortAD-CohortMCI,
    levels = colnames(temp_design))
  
  # Re-orientate the fitted model object from the coefficients of the original design matrix to
  # any set of contrasts of the original coefficients. 
  temp_fit2 <- limma::contrasts.fit(temp_fit, temp_contrast_matrix)
  
  Limmaefit = limma::eBayes(temp_fit2)
  
  for(loop_contrast in colnames(temp_contrast_matrix)){
    
    s_coeff = loop_contrast 
    DM_Results = limma::topTable(Limmaefit, coef=s_coeff, num=Inf, sort.by = "P", adjust="BH") # Data of specified variables
    print(head(DM_Results,20))
    DM_Results_sign=DM_Results[DM_Results$P.Value<0.01,]
    
    destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DM_hist_",loop_contrast, mod_type, ".pdf")
    pdf(file=destination)
    
    # histogram of p values and adjusted p values
    hist(DM_Results$P.Value ,xlim=c(0,1),main=s_coeff)
    hist(DM_Results$adj.P.Val ,xlim=c(0,1),main=s_coeff)
    min(DM_Results$adj.P.Val)
    cat(paste0("Contrast: ",s_coeff,".\n  number of features below TH (0.01)\n  Pvalue:",
               sum(DM_Results$P.Value<0.01)," (",round((sum(DM_Results$P.Value<0.01)/dim(DM_Results)[1])*100,0),"%)",
               "\n  AdjPvalue:",sum(DM_Results$adj.P.Val<0.05)," (",round((sum(DM_Results$adj.P.Val<0.05)/dim(DM_Results)[1])*100,0),"%)","\n\n"))
    dev.off()
    
    save(DM_Results,file = paste0(s_ROOT_dir,s_out_folder,"DM/DM_results_",loop_contrast, "_", mod_type,".Rdata")) 
    write.csv(DM_Results,file = paste0(s_ROOT_dir,s_out_folder,"DM/DM_results_",loop_contrast, "_", mod_type,".csv")) 
    save(DM_Results_sign,file = paste0(s_ROOT_dir,s_out_folder,"DM/DM_results_",loop_contrast,"_", mod_type,"_sign.Rdata")) 
  }
}

DM_analysis(methyl_Metas, "5mC")
gc()

DM_analysis(hydroxymethyl_Metas, "5hmC")
gc()

DM_analysis(unmod_Metas, "unmodC")
gc()


# Obtain all 5hmC Results
temp = list.files(path=paste0(s_ROOT_dir,s_out_folder,"DM"), pattern = "5hmC.csv")
DM_Results_5hmC <- as.list(temp)


DM_Results_5hmC <- lapply(temp, function(DM_result) {
  read.csv(paste0(s_ROOT_dir,s_out_folder,"DM/",DM_result), row.names=1)
})

# TO BE CONTINUED

#-----------------------------------------------------------------------------------------------------#
#							2. Enrichment at probe level of all significant
#-----------------------------------------------------------------------------------------------------#
# get the CpGs; significant
best_CpGs = rownames(DM_Results[DM_Results$adj.P.Val<0.05,])
all_CpGs = rownames(DM_Results)

head(all_CpGs)


#-----------------------------------------------------------------------------------------------------#
#							QQ plot of P values
#-----------------------------------------------------------------------------------------------------#
destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DM5_qqPlot.pdf", sep="")
pdf(file=destination)

library(qqman)
#qq(DM_Results$P.Value)
dev.off()


#-----------------------------------------------------------------------------------------------------#
#					Generate column with chr number and MAPINFO
#-----------------------------------------------------------------------------------------------------#

# Obtain genomic coordinates for filtered positions 
position <- data.frame(unlist(str_split_fixed(rownames(DM_Results), pattern = "_", n = 2)))
position <- cbind(position, "end_pos0" = position$start_pos0 + 1)

colnames(position) <- c("chrom", "start_pos0")
position$start_pos0 <- as.numeric(position$start_pos0)
chrom <- sub("_.*", "", sub("chr", "", rownames(DM_Results)))
DM_Results <- cbind(DM_Results, "MAPINFO" = position$start_pos0, "CHR" = chrom, "Name" = rownames(DM_Results))


#-----------------------------------------------------------------------------------------------------#
#                      Manhattan plot 
#-----------------------------------------------------------------------------------------------------#

manhattanplot <- function(DM_Result) {
DM_Results <- DM_Results[which(DM_Results$CHR != "Y"),] ## to also exclude the X chromosome repeat and edit this line

DM_Results$CHR <- as.character(DM_Results$CHR)
DM_Results$CHR[which(DM_Results$CHR == "X")] <- 23 ## to also recode Y chromosome repeat and edit this line
DM_Results$CHR <- as.numeric(DM_Results$CHR)
DM_Results <- DM_Results[which(DM_Results$CHR != ""),]

bonfP<-0.05/nrow(DM_Results)

destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DM Plots/DMP_4_Manhattan.png", sep="")
png(file=destination)

manhattan(DM_Results, p = "P.Value", snp ="Name", bp = "MAPINFO", chr = "CHR", 
          genomewide = -log10(bonfP), suggestiveline = -log10(5e-5), logp=TRUE, col=c("grey","#505567"), ylim=c(0,8), annotatePval = 0.05)
}

manhattanplot(DM_Results)

dev.off()


#-----------------------------------------------------------------------------------------------------#
#					     6. Boxplots signficant probes based on FDR adjusted (Cohort)
#-----------------------------------------------------------------------------------------------------#
# Order results table by p-value
DM_Results<-DM_Results[order(DM_Results$adj.P.Val),]

boxplot_FDR <- function(Metas) {
  #destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DMP_6_BoxplotGenomeWideSig_FDR.pdf", sep="")
  #pdf(file=destination, width = 10, height = 8)

  par(mfrow = c(2,3)) ## will plot 6 figures per page in 2 rows by 3 columns
  
  for(i in which(DM_Results$adj.P.Val < 0.05)[1:10]){
  boxplot(Metas[rownames(DM_Results)[i],] ~ pheno_subset$Cohort, ylab = "DNA Methylation (M value)", 
          main = rownames(DM_Results)[i],
          xlab = NULL) 
  mtext(paste("P = ", signif(DM_Results$adj.P.Val[i],3)), side = 3, adj = 1) ## add P value to plot rounded to 3 sf
  }

}

boxplot_FDR(methyl_Metas)
dev.off()


#-----------------------------------------------------------------------------------------------------#
#					     6. Boxplots signficant probes based on FDR adjusted (Braak stage)
#-----------------------------------------------------------------------------------------------------#

boxplot_FDR_braak <- function(Metas) {
  #destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DMP_6_BoxplotGenomeWideSig_FDR.pdf", sep="")
  #pdf(file=destination, width = 10, height = 8)
  
  par(mfrow = c(2,3)) ## will plot 6 figures per page in 2 rows by 3 columns
  
  for(i in which(DM_Results$adj.P.Val < 0.05)[1:10]){
    boxplot(Metas[rownames(DM_Results)[i],] ~ pheno_subset$Braak_score, ylab = "DNA Methylation (M value)", 
            main = rownames(DM_Results)[i],
            xlab = NULL) 
    mtext(paste("P = ", signif(DM_Results$adj.P.Val[i],3)), side = 3, adj = 1) ## add P value to plot rounded to 3 sf
  }
  
}

boxplot_FDR_braak(methyl_Metas)
dev.off()

#-----------------------------------------------------------------------------------------------------#
#					    7. Scatterplot signficant probes based on adj P value (Age)
#-----------------------------------------------------------------------------------------------------#

scatterplot_FDR_age <- function(Metas) {
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DMP_7_ScatterplotGenomeWideSig.pdf", sep="")
  pdf(file=destination, width = 10, height = 8)
  
  par(mfrow = c(2,3)) ## will plot 6 figures per page in 2 rows by 3 columns
  
  for(i in which(DM_Results$adj.P.Val < 0.05)[1:10]){
    plot(pheno_subset$Age, Metas[rownames(DM_Results)[i],] , xlab = "Age", ylab = "DNA methylation (M value)", pch = 16, main = rownames(DM_Results)[i]) 
    mtext(paste("P = ", signif(DM_Results$adj.P.Val[i],3)), side = 3, adj = 1) 
    
    # add P value to plot rounded to 3 sf
  }
  dev.off()
}

scatterplot_FDR_age(methyl_Metas)
dev.off()


#-----------------------------------------------------------------------------------------------------#
#					    7. Scatterplot signficant probes based on adj P value (DNA input)
#-----------------------------------------------------------------------------------------------------#

scatterplot_FDR_DNAinput <- function(Metas) {
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DMP_7_ScatterplotGenomeWideSig_DNA_Input.pdf", sep="")
  pdf(file=destination, width = 10, height = 8)
  
  par(mfrow = c(2,3)) ## will plot 6 figures per page in 2 rows by 3 columns
  
  for(i in which (DM_Results$adj.P.Val < 0.05)[1:10]){
    plot(pheno_subset$DNA_input, Metas[rownames(DM_Results)[i],] , xlab = "Age", ylab = "DNA methylation (M value)", pch = 16, main = rownames(DM_Results)[i]) 
    mtext(paste("P = ", signif(DM_Results$adj.P.Val[i],3)), side = 3, adj = 1) 
    
    # add P value to plot rounded to 3 sf
  }
  dev.off()
}

scatterplot_FDR_DNAinput(methyl_Metas)
dev.off()



#-----------------------------------------------------------------------------------------------------#
#					           8. Hierachical clustering
#-----------------------------------------------------------------------------------------------------#

# this code is clustering the 5000 probes with the largest SD i.e the most variable probes in this data set

hierarchicalcluster <- function(Metas) {
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DMP_8_HierachicalClustering.pdf", sep="")
  pdf(file=destination, width = 10, height = 8)

  sigma <- apply(Metas, 1, sd) # this calculates the standard deviation of all probes

  plot(hclust(dist(t(Metas[order(sigma, decreasing = TRUE)[1:5000],]))), 
     labels = paste(pheno_subset$Cohort, sep = "_"), cex = 0.68, main = "") 

  dev.off()
}

hierarchicalcluster(methyl_Metas)
dev.off()

#-----------------------------------------------------------------------------------------------------#
#					           9. Plot heatmaps
#-----------------------------------------------------------------------------------------------------#

heatmaps <- function(Metas) {
  sigma <- apply(Metas, 1, sd) # this calculates the standard deviation of all probes
  
  # Heatmap type 1
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DMP_9_Heatmap1.pdf", sep="")
  pdf(file=destination, width = 10, height = 8)
  
  heatmap.2(Metas[order(sigma, decreasing = TRUE)[1:5000],], trace = "none", 
            labCol = pheno_subset$CaseID, dendrogram = "column", labRow = "", mar = c(5,1), 
            density.info = "none", scale = "none")
  dev.off()
  
  #This code can be useful if have samples some with missing data
  #betas.com<-betas[complete.cases(betas),]
  
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DMP_9_Heatmap2.pdf", sep="")
  pdf(file=destination, width = 10, height = 8)
  
  
  # Heatmap type 2
  vec<-rep(NA, length(as.character(pheno_subset$Cohort)))
  vec[which(as.character(pheno_subset$Cohort)=="CTL")]<-"lightgrey"
  vec[which(as.character(pheno_subset$Cohort)=="MCI")]<-"black"
  vec[which(as.character(pheno_subset$Cohort)=="AD")]<-"#505567"
  
  my_palette11 <- colorRampPalette(c("darkblue", "white", "darkred"))(n = 1000)
  
  # lab_rows <- DM_Results[match(rownames(DM_Results), row.names(Metas[order(sigma, decreasing = TRUE)[1:5000],])),"genesUniq"][1:5000]
  lab_rows <- DM_Results[match(rownames(DM_Results), row.names(Metas[order(sigma, decreasing = TRUE)[1:5000],])), "Name"][1:5000]
  
  heatmap.2(Metas[order(sigma, decreasing = TRUE)[1:5000],], col=my_palette11, main="probes",keysize=2,key.par=list(cex=0.7, cex.main=1), trace= "none",dendrogram="column",tracecol=NA, cexCol=1.3,cexRow=0.7,ColSideColors = vec,margins=c(5,13), labCol=pheno_subset$CaseID)
  legend("topright", col = c("lightgrey", "black", "#505567"), pch = 16, c("CTL", "MCI", "AD"), horiz = FALSE, cex=0.8)
  
  dev.off()

}

heatmaps(methyl_Metas)
dev.off()
#-----------------------------------------------------------------------------------------------------#
#					Generate GRanges for included CpG positions
#-----------------------------------------------------------------------------------------------------#

# Obtain genomic coordinates for filtered positions 

position_cpg_gr <- GRanges(seqnames = position$chrom, ranges=IRanges(start = position$start_pos0, end = position$end_pos0))


#-----------------------------------------------------------------------------------------------------#
#                       annotation gene symbols
#-----------------------------------------------------------------------------------------------------#

annotations = build_annotations(genome = "hg38", annotations=c("hg38_cpgs","hg38_basicgenes"))


overlaps <- findOverlaps(position_cpg_gr, genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
overlapping_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)[subjectHits(overlaps)]

gene_symbols <- mapIds(org.Hs.eg.db, keys = as.character(overlapping_genes$gene_id), column = "SYMBOL",  keytype = "ENTREZID", multiVals = "first")

# Add gene symbols to `position_cpg_gr` element metadata based on overlap index
position_cpg_gr@elementMetadata$symbol <- NA  # Initialize with NA
position_cpg_gr@elementMetadata$symbol[queryHits(overlaps)] <- gene_symbols


position_cpg_gr

#-----------------------------------------------------------------------------------------------------#
#                       annotation regions
#-----------------------------------------------------------------------------------------------------#
bedfiles_gr_ann <- annotate_regions(bedfiles_gr[[1]], annotations = annotations)

summarize_annotations(annotated_regions = bedfiles_gr_ann, quiet = TRUE)

plot_annotation(annotated_regions = bedfiles_gr_ann,
                plot_title = '# of sites tested for all chromosomes',
                x_label = 'knownGene Annotations', y_label = 'Count')



#-----------------------------------------------------------------------------------------------------#
#                           RRMS BED file regions per chromosome
#-----------------------------------------------------------------------------------------------------#

library(chromoMap)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

bsg <- BSgenome.Hsapiens.UCSC.hg38

# Calculate length of each chromosome
length_chromosome <- function(chromosome) {
  seq <- getSeq(bsg, names = chromosome)
}

chromosomes <- seqnames(bsg)[1:24]
chrom_length <- sapply(chromosomes, length_chromosome)

# Create data frame with start and end of chromosomes
length <- data.frame(chromosomes, "start" = 1, "end" = sapply(chrom_length, length))

# Create anno data frame with RRMS bed regions
RRMS_bed <- as.data.frame(read.table(paste0(s_ROOT_dir, "Reference-files/RRMS_human_hg38.bed") ,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
colnames(RRMS_bed) <- c("Chr", "Pos_start", "Pos_end")

anno.data <- data.frame(paste0("region",rownames(RRMS_bed)))
anno.data <- cbind(anno.data, RRMS_bed)

# Plot RRMS BED file regions on top of chromosomes
chromoMap(list(length),list(anno.data))

dev.off()


#-----------------------------------------------------------------------------------------------------#
#                 Genome wide significant CpGs per chromosome
#-----------------------------------------------------------------------------------------------------#

# Obtain all genome wide significant CpGs 
DM_Results_subset <- DM_Results[DM_Results$adj.P.Val < 0.05,]
anno.data <- data.frame(rownames(DM_Results_subset))
anno.data <- cbind(anno.data, "CHR" = paste0("chr", DM_Results_subset$CHR), "Start" = DM_Results_subset$MAPINFO, "End" = DM_Results_subset$MAPINFO + 1, "logFC" = DM_Results_subset$logFC)

# recode chr23 as chrX
anno.data$CHR[which(anno.data$CHR=="chr23")] <- "chrX"

# Plot methylation value changes on top of chromosomes (per CHR)
anno.data.chr <- list(unique(anno.data$CHR))
length.chr <- list(unique(anno.data$CHR))

for (i in unique(anno.data$CHR)) {
  length.chr[[i]] <- length[length$chromosomes == i,]
  anno.data.chr[[i]] <- anno.data[anno.data$CHR == i,]
  
  # store image
  image <- chromoMap(list(length.chr[[i]]),list(anno.data.chr[[i]]),
          data_based_color_map = TRUE,
          data_type = "numeric",
          plots = "scatter",
          plot_filter = list(c("gt",0,"red")),
          ref_line = TRUE,
          refl_pos = 0,
          chr_width = 8,
          chr_length = 8,
          plot_height = 100,
          export.options = TRUE)
  
  # save image as html file
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DMP_10_5mC_Chrom", i, ".html", sep="")
  save_html(image, destination)
}


# # Plot methylation value changes on top of chromosomes (all CHR together)
image <- chromoMap(list(length),list(anno.data),
          data_based_color_map = TRUE,
          data_type = "numeric",
          plots = "scatter",
          plot_filter = list(c("lt",0,"red")),
          ref_line = TRUE,
          refl_pos = 0,
          chr_width = 8,
          chr_length = 8,
          plot_height = 100,
          export.options = TRUE)

destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DMP_10_5mC_allCHR.html", sep="")
save_html(image, destination)


# CHROM PLOT



#-----------------------------------------------------------------------------------------------------#
#                       GenomicDistributions
#-----------------------------------------------------------------------------------------------------#

queryList = GRangesList(position_cpg_gr)



TSSdist = calcFeatureDistRefTSS(queryList, "hg38")
p = plotFeatureDist(TSSdist, featureName="TSS", tile=TRUE, nbin=200)

plotFeatureDist(TSSdist, featureName="TSS")

print(p)


TSS <- GenomicDistributionsData::TSS_hg38()
featuredist = calcFeatureDist(position_cpg_gr, TSS)
featuredist


position <- cbind(position, "TSSdist" = unlist(Featuredist))
head(position)
perc2 = calcPartitionsRef(position_cpg_gr, "hg38")
plotPartitions(perc2)


#-----------------------------------------------------------------------------------------------------#
#                       Volcano plot (Enhanced volcano)
#-----------------------------------------------------------------------------------------------------#

volcanoplot <- function(mod_values) {
  controls <- paste0("sample_", pheno_subset$CaseID[which(pheno_subset$Cohort =="CTL")])
  controls <- which(colnames(mod_values) %in% controls)
  cases <- paste0("sample_", pheno$CaseID[which(pheno$Cohort =="AD")])
  cases <- which(colnames(mod_values) %in% cases)
  
  mod_values_x <- rowMeans(as.matrix(mod_values[, cases])) - rowMeans(as.matrix(mod_values[, controls]))
  mod_values_x <- as.data.frame(mod_values_x)
  
  
  cpg_site <- rownames(mod_values)
  mod_values_x[,2] <- cpg_site
  mod_values_x <- filter(mod_values_x, cpg_site %in% rownames(DM_Results))
  colnames(mod_values_x) <- c("DeltaBeta", "cpg_site")
  
  dm <- DM_Results[match(rownames(mod_values_x), rownames(DM_Results)),c("P.Value", "adj.P.Val", "Name")]
  mod_values_x <- cbind(mod_values_x, dm)
  
  mod_values_x$Name[which(mod_values_x$adj.P.Val> 0.05)] <- NA
  
  
  bonfP<-0.05/nrow(mod_values_x)
  
  keyvals.colour <- ifelse(
    mod_values_x$DeltaBeta < -0.1 & mod_values_x$adj.P.Val <0.05, 'darkblue',
    ifelse(mod_values_x$DeltaBeta > 0.1 & mod_values_x$adj.P.Val <0.05, 'darkred',
           '#E5E6EB'))
  keyvals.colour[is.na(keyvals.colour)] <- '#E5E6EB'
  names(keyvals.colour)[keyvals.colour == 'darkblue'] <- 'Hypomethylated'
  names(keyvals.colour)[keyvals.colour == 'darkred'] <- 'Hypermethylated'
  names(keyvals.colour)[keyvals.colour == '#E5E6EB'] <- 'NS'


  #destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/GO_Volcanoplot.pdf", sep="")
  #pdf(file=destination)
  
  EnhancedVolcano(mod_values_x, 
                  lab = mod_values_x$Name,
                  selectLab = mod_values_x$Name[which(mod_values_x$adj.P.Val < 0.05)],
                  x = "DeltaBeta", 
                  y = "P.Value", 
                  title = "", 
                  subtitle = "",
                  pCutoff = bonfP, 
                  FCcutoff = 0.1, 
                  pointSize = 3,
                  colCustom = keyvals.colour,
                  ylim = c(0, -log10(10e-9)),
                  xlim = c(-0.5, 0.5),
                  gridlines.major = FALSE,
                  gridlines.minor = FALSE,
                  xlab = "Delta beta value",
                  boxedLabels = TRUE,
                  drawConnectors = TRUE)
  #dev.off()
}

dev.off()


volcanoplot(methyl_Metas)


#-----------------------------------------------------------------------------------------------------#
#                       Volcano plot (Bonefrro)
#-----------------------------------------------------------------------------------------------------#


bonfP<-0.05/nrow(DM_Results)

DM_Results_subset <- DM_Results[which(DM_Results$P.Value < 0.01),]

keyvals.colour <- ifelse(
  DM_Results_subset$logFC < 0 & DM_Results_subset$P.Value <bonfP, 'darkblue',
  ifelse(DM_Results_subset$logFC > 0 & DM_Results_subset$P.Value <bonfP, 'darkred',
         '#E5E6EB'))
keyvals.colour[is.na(keyvals.colour)] <- '#E5E6EB'
names(keyvals.colour)[keyvals.colour == 'darkblue'] <- 'Hypomethylated'
names(keyvals.colour)[keyvals.colour == 'darkred'] <- 'Hypermethylated'
names(keyvals.colour)[keyvals.colour == '#E5E6EB'] <- 'NS'

DM_Results_subset <- cbind(DM_Results_subset, "change" = names(keyvals.colour))
head(DM_Results_subset)


ggplot(DM_Results_subset, aes(logFC, -log10(P.Value), color = change)) +
  geom_point(color = keyvals.colour) + 
  geom_hline(yintercept=-log10(bonfP), linetype = "dashed", color = "black") +
  theme_classic()


dev.off()
#-----------------------------------------------------------------------------------------------------#
#                       Gene set enrichment analysis (GSEA)
#-----------------------------------------------------------------------------------------------------#





#-----------------------------------------------------------------------------------------------------#
#                       Pathway enrichment analysis
#-----------------------------------------------------------------------------------------------------#













