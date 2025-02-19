#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#	Downsampling.R
#
#	Purpose 
#		This code was made for simulating coverage and downsampling
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
library(dplyr)
library("viridis")           

library(purrr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(readxl)
library(annotatr)
library(envnames)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


#-----------------------------------------------------------------------------------------------------#
#					            	Load bedmethyl data
#-----------------------------------------------------------------------------------------------------#
# Create list of bedmethyl files to import
temp = list.files(s_bedfiles_folder, pattern="\\.bed$")


bedfile <- as.data.frame(read.table(paste0(s_bedfiles_folder,"/",temp[1]), header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

bedfile <- setNames(bedfile, c("chrom", "start_pos0", "end_pos0", "mod_base_code","score",
                               "strand", "start_pos","end_pos", "color", "Nvalidcov", "percent_mod", 
                               "N_mod", "N_canon", "N_other_mod", "N_delete", "N_fail", "N_diff", "N_nocall"))

head(bedfile)

N_methyl <- bedfile[bedfile$mod_base_code=="m",c("chrom", "start_pos0", "end_pos0", "score", "N_mod")]
N_hydroxy <- bedfile[bedfile$mod_base_code=="h",c("chrom", "start_pos0", "end_pos0", "score", "N_mod")]
N_unmod <- bedfile[bedfile$mod_base_code=="m",c("chrom", "start_pos0", "end_pos0", "score", "N_canon")]

N_total <- cbind(N_methyl, "N_hydroxy" = N_hydroxy$N_mod, "N_canon" = N_unmod$N_canon)

dev.off()



#-----------------------------------------------------------------------------------------------------#
#					        Density plots of all CpGs vs. >different coverages CpGs 
#-----------------------------------------------------------------------------------------------------#

# set colors
library(RColorBrewer)
viridis_colors <- viridis(30)
brewer_colors <- brewer.pal(n = 10, name = "RdBu")

cov_threshold_interest <- c(4, 6, 8, 10, 12, 14, 16)

par(mfrow=c(1,3))

# unmod
plot(density(N_total$N_canon/N_total$score, na.rm=TRUE), col = "darkgrey", main = "all CpGs",
     xlab="Modifications (%)")     

for (c in 1:length(cov_threshold_interest)) {
  cov_threshold <- cov_threshold_interest[c] 
  N_total_10x <- filter(N_total, score >=cov_threshold)
  lines(density(N_total_10x$N_canon/N_total_10x$score, na.rm=TRUE), col = brewer_colors[c])
}


# methyl
plot(density(N_total$N_mod/N_total$score, na.rm=TRUE), col = "darkgrey", main = "5mC CpGs",
     xlab="Modifications (%)")     

for (c in 1:length(cov_threshold_interest)) {
  cov_threshold <- cov_threshold_interest[c] 
  N_total_10x <- filter(N_total, score >=cov_threshold)
  
  lines(density(N_total_10x$N_mod/N_total_10x$score, na.rm=TRUE), col = brewer_colors[c])
}

# hydroxymethyl
plot(density(N_total$N_hydroxy/N_total$score, na.rm=TRUE), col = "darkgrey", main = "5hmC CpGs",
     xlab="Modifications (%)")     

for (c in 1:length(cov_threshold_interest)) {
  cov_threshold <- cov_threshold_interest[c] 
  N_total_10x <- filter(N_total, score >=cov_threshold)
  
  lines(density(N_total_10x$N_hydroxy/N_total_10x$score, na.rm=TRUE), col = brewer_colors[c])
}

dev.off()

#-----------------------------------------------------------------------------------------------------#
#					        Density plots of all CpGs vs. >10X CpGs 
#-----------------------------------------------------------------------------------------------------#


# Plot densities of all values (unmod, 5mC, and 5hmC)
par(mfrow=c(1,2))
plot(density(N_total$N_canon/N_total$score, na.rm=TRUE), col = "darkgrey", main = "all CpGs",
     xlab="Modifications (%)")
lines(density(N_total$N_mod/N_total$score, na.rm=TRUE), col = "darkblue")
lines(density(N_total$N_hydroxy/N_total$score, na.rm=TRUE), col = "orange")

plot(density(N_total_10x$N_canon/N_total_10x$score, na.rm=TRUE), col = "darkgrey", 
     main = "Above 10X CpGs",
     xlab="Modifications (%)")
lines(density(N_total_10x$N_mod/N_total_10x$score, na.rm=TRUE), col = "orange")
lines(density(N_total_10x$N_hydroxy/N_total_10x$score, na.rm=TRUE), col = "darkblue")
#legend("topright", legend=c("unmodified", "5mC", "5hmC"), lty=1, col=c("darkgrey","darkblue", "orange"))


# Plot densities of all values per modification type (unmod, 5mC, and 5hmC)
par(mfrow=c(1,3))
plot(density(N_total$N_canon/N_total$score, na.rm=TRUE), col = "darkgrey", main = "canonical CpGs",
     xlab="Modifications (%)")
lines(density(N_total_10x$N_canon/N_total_10x$score, na.rm=TRUE), col = "darkblue")

plot(density(N_total$N_mod/N_total$score, na.rm=TRUE), col = "darkgrey", main = "5mC CpGs",
     xlab="Modifications (%)")
lines(density(N_total_10x$N_mod/N_total_10x$score, na.rm=TRUE), col = "darkblue")

plot(density(N_total$N_hydroxy/N_total$score, na.rm=TRUE), col = "darkgrey", main = "5hmC CpGs",
     xlab="Modifications (%)")
lines(density(N_total_10x$N_hydroxy/N_total_10x$score, na.rm=TRUE), col = "darkblue")

dev.off()


#-----------------------------------------------------------------------------------------------------#
#					        Downsampling per CpG
#-----------------------------------------------------------------------------------------------------#

mods <- list(nrow(N_total))

N_total <- cbind(N_total, "repeats" = NA)
head(N_total)

N_total$repeats <- apply(N_total, 1, function(x) {
  unlist(c(replicate(as.integer(x[5]), "h"), 
           replicate(as.integer(x[6]), "m"), 
           replicate(as.integer(x[7]), "u")))})

save(N_total, file = paste0(s_OUT_dir, "Bedmethyl/N_total.rdata")) 

apply(N_total, 1, function(x) x[5])


N_total_10x <- cbind(N_total_10x, "repeats" = NA)

N_total_10x$repeats <- apply(N_total_10x, 1, function(x) {
  unlist(c(replicate(as.integer(x[5]), "h"), 
           replicate(as.integer(x[6]), "m"), 
           replicate(as.integer(x[7]), "u")))})


for (c in 1:10) {
  
  cov_threshold <- 5
  # Indicate how many times the list should be sampled 
  number_repeats <- 10
  
  
  
  # paste(sample(unlist(N_total_10x$repeats[3]), 5))
sampling_results <- apply(N_total_10x,1, function(position){
    
    list <- unlist(position[8])
  
  actual_h <- as.numeric(position[6])/as.numeric(position[4])*100
  actual_m <- as.numeric(position[5])/as.numeric(position[4])*100
  actual_u <- as.numeric(position[7])/as.numeric(position[4])*100
  
  results <- matrix(NA, ncol = 3, nrow = number_repeats)
  results <- data.frame(results)
  colnames(results) <- c("h", "m", "u")
  
  sampled_strings <- replicate(number_repeats, paste(sample(list, cov_threshold), collapse = " "))
  
  methyl 
  results[,"m"] <- str_count(sampled_strings, "m")/cov_threshold*100)
  results[,"h"] <- str_count(sampled_strings, "h")/cov_threshold*100)
  results[,"u"] <- str_count(sampled_strings, "u")/cov_threshold*100)
  
# somethng with mean 
  results
})

head(sampling_results)
  # Density plots
  
  
  max_density <- max(c(unlist(density(results[,"u"])["y"]), unlist(density(results[,"m"])["y"]),
                       unlist(density(results[,"h"])["y"])))
  
  plot(density(results[,"u"], na.rm=TRUE), 
       main = c("Downsampling", paste(total_calls, "total calls,", 
                                      number_repeats, "repetitions,",
                                      cov_threshold, "coverage threshold")), col = "darkgrey", 
       xlim = c(0, 100), ylim=c(0,(max_density)),
       xlab="Modifications (%)")
  lines(density(results[,"m"], na.rm=TRUE), col = "darkblue")
  lines(density(results[,"h"], na.rm=TRUE), col = "orange")
  abline(v = actual_h, col = "orange", lty ="dashed" )
  abline(v = actual_m, col = "darkblue", lty="dashed")
  abline(v=actual_u, col = "darkgrey", lty="dashed")
  legend("topright", legend=c("unmodified", "5mC", "5hmC"), lty=1, col=c("darkgrey","darkblue", "orange"))
  
}





# Take bedmethyl files N counts 

# Take 10 sa

head(hydroxymethyl)

library(stringr)

list <- c(replicate(3, "h"), replicate(10, "m"), replicate (7, "u"))
total_calls <- as.numeric(length(list))

# Set how many reads to sample (threshold high confidence calls) 
dev.off()



for (c in 1:10) {

cov_threshold <- c
# Indicate how many times the list should be sampled 
number_repeats <- 10000

actual_h <- str_count(paste(list, collapse= " "), "h")/total_calls*100
actual_m <- str_count(paste(list, collapse= " "), "m")/total_calls*100
actual_u <- str_count(paste(list, collapse= " "), "u")/total_calls*100

results <- matrix(NA, ncol = 3, nrow = number_repeats)
results <- data.frame(results)
colnames(results) <- c("h", "m", "u")

for (i in 1:number_repeats) {
  sampling <- paste(sample(list, cov_threshold), collapse = " ")
  results[i,"m"] <- str_count(sampling, "m")/cov_threshold*100
  results[i,"h"] <- str_count(sampling, "h")/cov_threshold*100
  results[i,"u"] <- str_count(sampling, "u")/cov_threshold*100
}

# Density plots


max_density <- max(c(unlist(density(results[,"u"])["y"]), unlist(density(results[,"m"])["y"]),
                     unlist(density(results[,"h"])["y"])))

plot(density(results[,"u"], na.rm=TRUE), 
     main = c("Downsampling", paste(total_calls, "total calls,", 
                                              number_repeats, "repetitions,",
                                              cov_threshold, "coverage threshold")), col = "darkgrey", 
     xlim = c(0, 100), ylim=c(0,(max_density)),
     xlab="Modifications (%)")
lines(density(results[,"m"], na.rm=TRUE), col = "darkblue")
lines(density(results[,"h"], na.rm=TRUE), col = "orange")
abline(v = actual_h, col = "orange", lty ="dashed" )
abline(v = actual_m, col = "darkblue", lty="dashed")
abline(v=actual_u, col = "darkgrey", lty="dashed")
legend("topright", legend=c("unmodified", "5mC", "5hmC"), lty=1, col=c("darkgrey","darkblue", "orange"))

}




#max_frequency <- max(c(results[,"u"], results[,"m"], results[,"h"]))

#hist(results[,"u"], xlim=c(0,100), ylim=c(0,max_density), breaks = 10, freq=FALSE, col = "darkgrey")
#hist(results[,"m"], xlim=c(0,100), ylim=c(0,max_density), breaks = 10, freq=FALSE, col = "darkblue", alpha = 0.1, add = TRUE)
#hist(results[,"m"], xlim=c(0,100), ylim=c(0,max_density), breaks = 10, freq=FALSE, col = rgb(1, 0, 0, 0.5), add = TRUE)
#lines(density(results[,"u"], na.rm=TRUE), col = "darkblue")


