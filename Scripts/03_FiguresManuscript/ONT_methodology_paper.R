
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
library(car)
library(emmeans)

library(ggprism)


#-----------------------------------------------------------------------------------------------------#
#						Load pheno data                                                           ----
#-----------------------------------------------------------------------------------------------------#
# Set cutoff coverage
cutoff = 5

# Pheno file
pheno_subset_pass <- readRDS(file = paste0(s_ROOT_dir, s_out_folder,"QC/Pheno_new_",cutoff,"X.rds"))

#-----------------------------------------------------------------------------------------------------#
#					FANS AND SEQUENCING                                                      ----
#-----------------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------------------------#
#						1. Nuclei RNA DNA yield                                                         ----
#-----------------------------------------------------------------------------------------------------#

# Load FANS data 
fans <- read_excel(paste0(s_ROOT_dir, "Anno/FANS_DNA_RNA.xlsx"), na = c("NA", "na", "Na", "nA", "?"),
                   col_types = c("numeric", "text", "numeric","numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric", "skip", "numeric"), sheet = 1)

fans$CaseID <- sub("_", "-", fans$CaseID)

# Load DNA isolation data (nanodrop)
fans_dna <- read_excel(paste0(s_ROOT_dir, "Anno/FANS_DNA_RNA.xlsx"), 
                       col_types = c("text", "numeric", "numeric", "numeric", "numeric"),
                       sheet = 2)

fans_dna$CaseID[duplicated(fans_dna$CaseID)]

# Combine data frames
fans <- left_join(fans, fans_dna, by = "CaseID")

# Combine with phenotype data (metadata)
load(file = paste0(s_ROOT_dir, s_out_folder, "Pheno/Pheno.Rdata"))
fans_ext <- left_join(fans, pheno[,c("CaseID", "Gender", "Age", "PMI", "Braak_score", "Group", "DNA_input")], by = "CaseID")
fans_ext <- fans_ext[fans_ext$CaseID != "16-51",]

# Check for duplicated cases
fans$CaseID[duplicated(fans$CaseID)]


# Load OLIG fragmentation and DNA 
fans_dna_olig <- read_excel(paste0(s_ROOT_dir, "Anno/Sample_fragmentation_qubit_Olig.xlsx"),
                   col_types = c("text", "text", "skip","numeric", "numeric", "numeric",
                                 "numeric"), sheet = 2)

colnames(fans_dna_olig) <- c("Group", "CaseID", "DNA_conc", "Mean_frag_size", "Total_DNA_Olig", "N_Flowcell")
fans_dna_olig$CaseID <- sub("_", "-", fans_dna_olig$CaseID)


# Load OLIG DIN
fans_DIN_olig <- read_excel(paste0(s_ROOT_dir, "Anno/Sample_fragmentation_qubit_Olig.xlsx"),
                            col_types = c("text", "text", "skip","skip", "skip", "skip", "numeric",
                                          "skip", "skip", "skip", "skip", "skip", "skip", "skip", "skip"), sheet = 1)

#-----------------------------------------------------------------------------------------------------#
#						2. Extra variables                                                     ----
#-----------------------------------------------------------------------------------------------------#

# Total number of nuclei
Total_nuclei <- apply(fans_ext, 1, function(x) sum(as.numeric(x["Total_nuclei_NeuN"]), as.numeric(x["Total Olig NeuN"]),as.numeric(x["Total_rest_nuclei"] ), na.rm = TRUE))
fans_ext <- cbind(fans_ext, Total_nuclei)

# Ratio of nuclei per population
Ratio_nuclei_NeuN <- fans_ext$Total_nuclei_NeuN / fans_ext$Total_nuclei * 100
Ratio_nuclei_Olig2 <- fans_ext$`Total Olig NeuN` / fans_ext$Total_nuclei * 100
Ratio_nuclei_Rest <- fans_ext$Total_rest_nuclei / fans_ext$Total_nuclei * 100

# Total amount of DNA
total_dna <- apply(fans_ext, 1, function(x) sum(as.numeric(x["Yield_NeuN"]), as.numeric(x["Yield_Olig2"]),as.numeric(x["Yield_Rest"] ), na.rm = TRUE))

# Create variable with dates 
year <- as.numeric(str_split_fixed(fans_ext$CaseID, "-",n = 2)[,1])
sample_date <- ifelse(year > 90, year + 1900, year + 2000)

# Add variables to fans_ext
fans_ext <- cbind(fans_ext, Ratio_nuclei_NeuN, Ratio_nuclei_Olig2, Ratio_nuclei_Rest, total_dna, sample_date)


#-----------------------------------------------------------------------------------------------------#
#					3. Scatterplot and boxplot function                                              ----
#-----------------------------------------------------------------------------------------------------#
# my_colors <- brewer.pal(4, "Blues")[4:2]
# scale_color_discrete(labels = c("CTL", "MCI", "AD")) +
# scale_color_paletteer_d("lisa::MarcelDuchamp") +
# scale_color_paletteer_d("fishualize::Acanthisthius_brasilianus", labels = c("CTL", "MCI", "AD")) 


my_colors <- c("#084594","#2171B5","#6BAED6", "#718FAB", "#343434", "#C6DBEF")
# 
# my_colors <- c(	"#343434", "#084594","#54A1E4", "#5EB373")
# 
# my_colors <- c(	"#161A53", "#084594","#2171B5", "#6BAED6")


# my_colors <- c(	"#343434", "#084594","#2171B5","#78B38C")

# my_colors <- c("#084594","#54A1E4","#446750", "#78B38C")



my_shapes <- c(16, 15, 17, 18)

# Scatterplot function
scatterplot_pheno <- function(data, x_var, y_var, group_var, label_text, xlab_text = "", ylab_text, facet_var, facet_var2) {
  
  destination <- file.path(s_OUT_dir, paste0("Plots/Pheno_",x_var,"_",y_var, "_", group_var,"_corr.svg"))
  svglite(file=destination)
  
  if (is.null(group_var)) {
    
    image <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
      geom_smooth(method = "lm", se = TRUE, colour = "#084594", fill="#C6DBEF", alpha = 0.5, size = 1) +  # will be colored by group
      geom_point(colour = "#084594", size = 3, alpha = 0.8) + 
      theme(legend.position = "right") +
      theme_bw() +
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      # stat_cor(method = "pearson") +  # per group correlation   
      sm_statCorr(fit.params = list(color = "#084594", linetype = 1, alpha = 0.8)) +  
      xlab(xlab_text) + 
      ylab(ylab_text)

    
    
  } else if (!is.null(facet_var2)) {
    
    image <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], 
                              color = .data[[group_var]], shape = .data[[group_var]])) +
      geom_smooth(method = "lm", se = TRUE, fill="#C6DBEF", alpha = 0.5, size = 1) +  # will be colored by group
      geom_point(alpha = 0.8) + 
      theme_bw() +
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      sm_statCorr() +  
      theme(legend.position = "right") +
      scale_color_manual(name = group_var, values = my_colors, labels = label_text) +
      scale_shape_manual(name = group_var, values = my_shapes, labels = label_text) +
      
      xlab(xlab_text) + 
      ylab(ylab_text) + 
      # facet_wrap(~.data[[facet_var]] + .data[[facet_var2]]) + 
      facet_grid(rows = vars(.data[[facet_var]]),
                 cols = vars(.data[[facet_var2]]), scales = "free")
    
  } else if (!is.null(facet_var)) {
    
    image <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], 
                              color = .data[[group_var]], shape = .data[[group_var]])) +
      geom_smooth(method = "lm", se = TRUE, fill="#C6DBEF", alpha = 0.5, size = 1) +  # will be colored by group
      geom_point(alpha = 0.8) + 
      theme_bw() +
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      sm_statCorr() +  
      theme(legend.position = "none") +
      scale_color_manual(name = group_var, values = my_colors, labels = label_text) +
      scale_shape_manual(name = group_var, values = my_shapes, labels = label_text) +
      
      xlab(xlab_text) + 
      ylab(ylab_text) + 
      facet_wrap(~.data[[facet_var]])
    
  } else {
    
    image <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
      theme_bw() +
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
           ) +
      sm_statCorr() +  
      theme(legend.position = "right") +
      scale_color_manual(name = group_var, values = c("#084594","#F8A058",  "#5EB373", "#DE9151", "#54A1E4"), labels = label_text) +
      scale_shape_manual(name = group_var, values = my_shapes, labels = label_text) +
      geom_smooth(method = "lm", se = TRUE, colour = alpha("#3B3B3B", 1.0), fill="#C6DBEF", alpha = 0.4, linewidth=0.8, n =5) +  # will be colored by group
      geom_point(aes(color = .data[[group_var]], shape = .data[[group_var]]), size =3, alpha = 1.0) + 
      
      xlab(xlab_text) + 
      ylab(ylab_text)
    
    
  }
    
  
  print(image)
  dev.off()
  
  cat("Saved plot to: ", destination, "\n")
  print(image)
  
  return(image)
}

# Boxplot function
boxplot_pheno <- function(data, x_var, y_var, group_var, label_text, xlab_text = "", ylab_text) {
  
  
  destination <- file.path(s_OUT_dir, paste0("Plots/Pheno_boxplot_",x_var,"_",y_var,".svg"))
  svglite(file=destination)
  
  if (is.null(group_var)) {
    
    df_p_val <- data %>%
      rstatix::t_test(as.formula(paste0(y_var, " ~ ", x_var))) %>%
      rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
      rstatix::add_significance(p.col = "p.adj") %>%
      rstatix::add_xy_position(x = x_var)
    
    set.seed(32)
    
    image <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
      geom_boxplot(fill = "white") + 
      geom_jitter(width=0.2, colour = "#2171B5")+
      theme_gray()+
      theme(legend.position = "right") +
      xlab(xlab_text) + 
      ylab(ylab_text) + 
      add_pvalue(df_p_val, y.position = "y.position", colour = "black", tip.length= 0)
    
    
  } else  {
    
    # Get significance values
    if (x_var == group_var) {
      
      df_p_val2 <- data %>%
        rstatix::wilcox_test(as.formula(paste0(y_var, " ~ ", group_var))) %>%
        rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
        rstatix::add_significance(p.col = "p.adj") %>%
        rstatix::add_xy_position(x = x_var, dodge=0.8)
      
      print(df_p_val2)
      
      
    } else {
      
      df_p_val <- data %>%
        rstatix::group_by(.data[[group_var]]) %>%
        rstatix::t_test(as.formula(paste0(y_var, " ~ ", x_var))) %>%
        rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
        rstatix::add_significance(p.col = "p.adj") %>%
        rstatix::add_xy_position(x = group_var, dodge = 0.8) # important for positioning!
      
      df_p_val2 <- data %>%
        rstatix::group_by(.data[[x_var]]) %>%
        rstatix::t_test(as.formula(paste0(y_var, " ~ ", group_var))) %>%
        rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
        rstatix::add_significance(p.col = "p.adj") %>%
        rstatix::add_xy_position(x = x_var, dodge = 0.8)
      
      print(df_p_val)
      print(df_p_val2)
    }
    


   
    names(df_p_val2)[3] <- group_var
    
    # max_dep <- max(data[,y_var], na.rm=TRUE)
    # seq <- round(max_dep / 10, digits = 0)
    
    # df_p_val <- df_p_val[df_p_val$p.adj.signif != "ns",]
    # df_p_val <- cbind(df_p_val, "y.position" = seq(max_dep, by = seq, length.out = nrow(df_p_val)))
    


    
    set.seed(32)
    

    image <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], 
                              colour = .data[[group_var]], shape = .data[[group_var]])) +
      geom_boxplot(fill = "white") + 
      # geom_jitter(width=0.2, size = 2)+
      geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.8) +
      theme_gray()+
      # theme(axis.line = element_line(color='black'),
      #       plot.background = element_blank(),
      #       panel.grid.major = element_blank(),
      #       panel.grid.minor = element_blank()) +
      theme(legend.position = "right") +
   
      xlab(xlab_text) + 
      ylab(ylab_text) +  
      # add_pvalue(df_p_val, y.position = "y.position", colour = "black",
      #            tip.length= 0, bracket.nudge.y = 2,
      #            xmin = "xmin",
      #            xmax = "xmax") +
      # 
      add_pvalue(df_p_val2, y.position = "y.position", colour = "black",
                 xmin = "xmin", xmax = "xmax", tip.length = 0, bracket.nudge.y = 2)+
      
      scale_color_manual(name = group_var, values = my_colors, labels = label_text) +
      scale_shape_manual(name = group_var, values = my_shapes, labels = label_text) 

    
  }
  
  print(image)
  dev.off()
  
  cat("Saved plot to: ", destination, "\n")
  
  print(image)
  return(image)
}



# Boxplot function
barplot_pheno <- function(data, x_var, y_var, group_var, label_text, xlab_text = "", ylab_text) {
  
  
  destination <- file.path(s_OUT_dir, paste0("Plots/Pheno_barplot_",x_var,"_",y_var,".svg"))
  svglite(file=destination)
  
  if (is.null(group_var)) {
    
    # Calculates mean, sd, se and IC
    my_sum <- data %>%
      group_by(.data[[x_var]]) %>%
      summarise( 
        N_set = length(.data[[y_var]]),
        mean = round(mean(.data[[y_var]], na.rm = TRUE), digits = 3),
        sd = round(sd(.data[[y_var]], na.rm = TRUE), digits = 3),
        se = round(sd(.data[[y_var]], na.rm = TRUE)/sqrt(length((.data[[y_var]]))), digits= 3)
        ) 


    df_p_val <- data %>%
      rstatix::t_test(as.formula(paste0(y_var, " ~ ", x_var))) %>%
      rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
      rstatix::add_significance(p.col = "p.adj") %>%
      rstatix::add_xy_position(x = x_var)
    
    set.seed(32)
    
    image <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
      geom_bar(stat = "identity", fill = "white") + 
      geom_jitter(width=0.2, colour = "#2171B5")+
      theme_gray()+
      theme(legend.position = "right") +
      xlab(xlab_text) + 
      ylab(ylab_text) + 
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.4, colour="black", alpha=0.9, size=1.5) +
      add_pvalue(df_p_val, y.position = "y.position", colour = "black", tip.length= 0)
    
    
  } else  {
    
    # Calculates mean, sd, se and IC
    my_sum <- data %>%
      group_by(.data[[group_var]]) %>%
      summarise( 
        N_set = length(.data[[y_var]]),
        mean = round(mean(.data[[y_var]], na.rm = TRUE), digits = 3),
        sd = round(sd(.data[[y_var]], na.rm = TRUE), digits = 3),
        se = round(sd(.data[[y_var]], na.rm = TRUE)/sqrt(length((.data[[y_var]]))), digits= 3),
        .groups = "drop"
      )
   
   
    # Get significance values
    if (x_var == group_var) {
      
      df_p_val2 <- data %>%
        rstatix::t_test(as.formula(paste0(y_var, " ~ ", group_var))) %>%
        rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
        rstatix::add_significance(p.col = "p.adj") %>%
        rstatix::add_xy_position(x = x_var, dodge=0.8)
      
      print(df_p_val2)
      
      
    } else {
      
      # Calculates mean, sd, se and IC
      my_sum <- data %>%
        group_by(.data[[x_var]], .data[[group_var]]) %>%
        summarise( 
          N_set = length(.data[[y_var]]),
          mean = round(mean(.data[[y_var]], na.rm = TRUE), digits = 3),
          sd = round(sd(.data[[y_var]], na.rm = TRUE), digits = 3),
          se = round(sd(.data[[y_var]], na.rm = TRUE)/sqrt(length((.data[[y_var]]))), digits= 3),
          .groups = "drop"
        )
      
      print(my_sum)
      
      df_p_val <- data %>%
        rstatix::group_by(.data[[group_var]]) %>%
        rstatix::t_test(as.formula(paste0(y_var, " ~ ", x_var))) %>%
        rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
        rstatix::add_significance(p.col = "p.adj") %>%
        rstatix::add_xy_position(x = group_var, dodge = 0.8) # important for positioning!
      
      df_p_val2 <- data %>%
        rstatix::group_by(.data[[x_var]]) %>%
        rstatix::t_test(as.formula(paste0(y_var, " ~ ", group_var))) %>%
        rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
        rstatix::add_significance(p.col = "p.adj") %>%
        rstatix::add_xy_position(x = x_var, dodge = 0.8)
      
      print(df_p_val)
      print(df_p_val2)
    }
    
    
    names(df_p_val2)[3] <- group_var
    
    set.seed(32)
    
    image <- ggplot() +
      geom_bar(data = my_sum, aes(x = .data[[x_var]], y = mean, color = .data[[group_var]]), fill = "white",
                                  stat= "identity",position=position_dodge(width = 0.8), width = 0.6, linewidth = 0.5) + 
      # geom_jitter(width=0.2, size = 2)+
      geom_point(data = data, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[group_var]], shape = .data[[group_var]]),
                                  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +
      theme_gray()+
      # theme(axis.line = element_line(color='black'),
      #       plot.background = element_blank(),
      #       panel.grid.major = element_blank(),
      #       panel.grid.minor = element_blank()) +
      theme(legend.position = "right") +
      xlab(xlab_text) + 
      ylab(ylab_text) +  
      geom_errorbar(data = my_sum, aes(x = .data[[x_var]], ymin=mean-se, ymax=mean+se, group = .data[[group_var]], color = .data[[group_var]]), 
                    width=0.2, linewidth = 1.0, alpha=1.0, position = position_dodge((width = 0.8))) +
      # add_pvalue(df_p_val, y.position = "y.position", colour = "black",
      #            tip.length= 0, bracket.nudge.y = 2,
      #            xmin = "xmin",
      #            xmax = "xmax") +

      add_pvalue(df_p_val2, y.position = "y.position", colour = "black",
                 xmin = "xmin", xmax = "xmax", tip.length = 0, bracket.nudge.y = 2) +
      scale_fill_manual(name = group_var, values = my_colors, labels = label_text) +
      scale_shape_manual(name = group_var, values = my_shapes, labels = label_text) +
      scale_color_manual(name = group_var, values = my_colors, labels = label_text)
      

  }
  
  print(image)
  dev.off()
  
  cat("Saved plot to: ", destination, "\n")
  
  print(image)
  return(image)
}

#-----------------------------------------------------------------------------------------------------#
#					3b. Statistics function                                            ----
#-----------------------------------------------------------------------------------------------------#

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
    write.table(summary_df, file = file_output, append = TRUE, quote = FALSE, sep = "\t")
    
    
    # Pairwise comparisons
    cat(paste0("Pairwise comparisons (Sidak): ~ ", indep_var1, "\n"))
    formula <- as.formula(paste("~", indep_var1))
    em <- emmeans(mod, formula)
    em_pairs <- pairs(em, adjust = "sidak")
    
    write.table(em_pairs, file = file_output, append = TRUE, quote = FALSE, sep = "\t")
    
    # Data for plot
    em_pairs <- pairs(em, adjust = "sidak")
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

cor_stats_txt <- function(data, y_var, x_var, group_var, group_var2) {
  
  if(is.null(group_var) & is.null(group_var2)) {
    
    cor_df <- data %>%
      summarise(
        cor_res <- broom::tidy(cor.test(.data[[x_var]], .data[[y_var]], method = "pearson")),
        .groups = "drop"
      ) %>% 
      mutate(
        estimate = round(estimate, 3),
        statistic = round(statistic, 3),
        p.value = round(p.value, 3),
        parameter = parameter,
        conf.low = round(conf.low, 3),
        conf.high = round(conf.high, 3),
        method = NULL,
        alternative = NULL
      )
    
    print(cor_df)
    
    file_output = file.path(s_OUT_dir, "QC", paste0("QC_stats_", y_var, "_", x_var, ".txt"))
    write.table(cor_df, file = file_output, quote = FALSE, sep = "\t")
    
  } else if(is.null(group_var2)) {
    
    cor_df <- data %>%
      group_by(.data[[group_var]]) %>%
      summarise(
        cor_res <- broom::tidy(cor.test(.data[[x_var]], .data[[y_var]], method = "pearson")),
        .groups = "drop"
      ) %>% 
      mutate(
        estimate = round(estimate, 3),
        statistic = round(statistic, 3),
        p.value = round(p.value, 3),
        parameter = parameter,
        conf.low = round(conf.low, 3),
        conf.high = round(conf.high, 3),
        method = NULL,
        alternative = NULL
      )
    
    print(cor_df)
    
    file_output = file.path(s_OUT_dir, "QC", paste0("QC_stats_", y_var, "_", x_var, "_",group_var, ".txt"))
    write.table(cor_df, file = file_output, quote = FALSE, sep = "\t")
    
  } else {
    
    cor_df <- data %>%
      group_by(.data[[group_var]], .data[[group_var2]]) %>%
      summarise(
        cor_res <- broom::tidy(cor.test(.data[[x_var]], .data[[y_var]], method = "pearson")),
        .groups = "drop"
      ) %>% 
      mutate(
        estimate = round(estimate, 3),
        statistic = round(statistic, 3),
        p.value = round(p.value, 3),
        parameter = parameter,
        conf.low = round(conf.low, 3),
        conf.high = round(conf.high, 3),
        method = NULL,
        alternative = NULL
      )
    
    print(cor_df)
    
    file_output = file.path(s_OUT_dir, "QC", paste0("QC_stats_", y_var, "_", x_var, "_",group_var, "_", group_var2, ".txt"))
    write.table(cor_df, file = file_output, quote = FALSE, sep = "\t")
  }
  
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



#-----------------------------------------------------------------------------------------------------#
#					4. Nuclei yield per type vs group                                            ----
#-----------------------------------------------------------------------------------------------------#

# Data with nuclei yield: total, NeuN+, Olig2+ and double negative 
total_nuclei_data <- list(as.numeric(fans_ext$Total_nuclei),
                          as.numeric(fans_ext$Total_nuclei_NeuN),
                          as.numeric(fans_ext$`Total Olig NeuN`),
                         as.numeric(fans_ext$Total_rest_nuclei))

names(total_nuclei_data) <- c("Total","NeuN+", "Olig2+", "NeuN-/Olig2-")

# Create data frame with gene and contrast
total_nuclei_data <- bind_rows(
  lapply(names(total_nuclei_data), function(nuclei_type) {
    tibble(Nuclei = total_nuclei_data[[nuclei_type]], 
           Set = nuclei_type, 
           Group= fans_ext$Group, 
           PMI = fans_ext$PMI, 
           sample_date = fans_ext$sample_date)}))

# Omit rows with NA values
total_nuclei_data <- total_nuclei_data[!is.na(total_nuclei_data$Nuclei),]

# Factorize the nuclei type variable
total_nuclei_data$Set <- factor(total_nuclei_data$Set, levels = c("Total", "NeuN+", "Olig2+", "NeuN-/Olig2-"))

# Plot nuclei yield per type per group
p_Nuclei_vs_Group <- boxplot_pheno(total_nuclei_data, 
                                x_var = "Set",
                                y_var="Nuclei", 
                                group_var = "Group", 
                                label_text = c("CTL", "MCI", "AD"), 
                                xlab_text = "",
                                ylab_text = "Total isolated nuclei (x 1000)")

# Plot nuclei yield per type per group
p_Nuclei_vs_Group <- boxplot_pheno(total_nuclei_data, 
                                   x_var = "Group",
                                   y_var="Nuclei", 
                                   group_var = "Set", 
                                   label_text = c("Total", "NeuN+", "Olig2+", "NeuN-/Olig2-"), 
                                   xlab_text = "",
                                   ylab_text = "Total isolated nuclei (x 1000)")



p_Nuclei_vs_Group_bar <- barplot_pheno(total_nuclei_data, 
                                   x_var = "Set",
                                   y_var="Nuclei", 
                                   group_var = "Group", 
                                   label_text = c("CTL", "MCI", "AD"), 
                                   xlab_text = "",
                                   ylab_text = "Total isolated nuclei (x 1000)")


p_Nuclei_vs_Group_bar <- barplot_pheno(total_nuclei_data, 
                                       x_var = "Group",
                                       y_var="Nuclei", 
                                       group_var = "Set", 
                                       label_text = c("Total", "NeuN+", "Olig2+", "NeuN-/Olig2-"), 
                                       xlab_text = "",
                                       ylab_text = "Total isolated nuclei (x 1000)")

# p_Nuclei_vs_Group_set <- boxplot_pheno(total_nuclei_data, 
#                                    x_var = "Group",
#                                    y_var="Nuclei", 
#                                    group_var = "Set", 
#                                    label_text = c("NeuN+", "Olig2+", "NeuN-/Olig2-"), 
#                                    xlab_text = "",
#                                    ylab_text = "Total isolated nuclei (x 1000)")


# Statistics
get_stats_txt(total_nuclei_data, depend_var = "Nuclei", "Group", "Set")



#-----------------------------------------------------------------------------------------------------#
#					5. Nuclei yield vs PMI and vs. sample acquisition                                           ----
#-----------------------------------------------------------------------------------------------------#


# Plot the correlation between nuclei yield and PMI
total_nuclei_data_PMI <- total_nuclei_data[total_nuclei_data$PMI < 10, ] 


p_Nuclei_vs_PMI <- scatterplot_pheno(total_nuclei_data_PMI, 
                                          x_var = "PMI",
                                          y_var="Nuclei", 
                                          group_var = "Set", 
                                          label_text = c("NeuN+", "Olig2+", "NeuN-/Olig2-"), 
                                          xlab_text ="Post-mortem interval (PMI, in hours)",
                                          ylab_text = "Total isolated nuclei",
                                          facet_var = "Set",
                                          facet_var2 = NULL)

cor_stats_txt(total_nuclei_data, x_var = "PMI",
              y_var="Nuclei", 
              group_var = "Set",
              group_var2 = NULL)

# Plot the correlation between nuclei yield and sample acquisition
p_Nuclei_vs_sampledate <- scatterplot_pheno(total_nuclei_data, 
                                     x_var = "sample_date",
                                     y_var="Nuclei", 
                                     group_var = "Set", 
                                     label_text = c("NeuN+", "Olig2+", "NeuN-/Olig2-"), 
                                     xlab_text ="Year of sample acquisition",
                                     ylab_text = "Total isolated nuclei",
                                     facet_var = "Set",
                                     facet_var2 = NULL)


cor_stats_txt(total_nuclei_data, x_var = "sample_date",
              y_var="Nuclei", 
              group_var = "Set",
              group_var2 = NULL)


#-----------------------------------------------------------------------------------------------------#
#					6. Fraction of nuclei type vs group                                            ----
#-----------------------------------------------------------------------------------------------------#

# Fraction nuclei
fraction_nuclei_data <- list(as.numeric(fans_ext$Ratio_nuclei_NeuN),
                          as.numeric(fans_ext$Ratio_nuclei_Olig2),
                          as.numeric(fans_ext$Ratio_nuclei_Rest))

names(fraction_nuclei_data) <- c("NeuN+", "Olig2+", "NeuN-/Olig2-")

# Create data frame with gene and contrast
fraction_nuclei_data <- bind_rows(
  lapply(names(fraction_nuclei_data), function(nuclei_type) {
    tibble(Fraction = fraction_nuclei_data[[nuclei_type]], Set = nuclei_type, Group= fans_ext$Group)
  }))

# Omit NA values
fraction_nuclei_data <- fraction_nuclei_data[!is.na(fraction_nuclei_data$Fraction),]

# Factorize nuclei type variable 
fraction_nuclei_data$Set <- factor(fraction_nuclei_data$Set, levels = c("NeuN+", "Olig2+", "NeuN-/Olig2-"))

# Plot fraction of nuclei type vs group
p_Fraction_vs_Group <- boxplot_pheno(fraction_nuclei_data, 
                                   x_var = "Set",
                                   y_var="Fraction", 
                                   group_var = "Group", 
                                   label_text = c("CTL", "MCI", "AD"), 
                                   xlab_text = "",
                                   ylab_text = "Percentage of total nuclei (%)")

p_Fraction_vs_Group <- boxplot_pheno(fraction_nuclei_data, 
                                     x_var = "Group",
                                     y_var="Fraction", 
                                     group_var = "Set", 
                                     label_text = c("NeuN+", "Olig2+", "NeuN-/Olig2-"), 
                                     xlab_text = "",
                                     ylab_text = "Percentage of total nuclei (%)")



my_colors <- c("#2171B5","#6BAED6", "#718FAB", "#343434", "#C6DBEF")

p_Fraction_vs_Group_bar <- barplot_pheno(fraction_nuclei_data, 
                                     x_var = "Set",
                                     y_var="Fraction", 
                                     group_var = "Group", 
                                     label_text = c("CTL", "MCI", "AD"), 
                                     xlab_text = "",
                                     ylab_text = "Percentage of total nuclei (%)")


p_Fraction_vs_Group_bar <- barplot_pheno(fraction_nuclei_data, 
                                         x_var = "Group",
                                         y_var="Fraction", 
                                         group_var = "Set", 
                                         label_text = c("NeuN+", "Olig2+", "NeuN-/Olig2-"), 
                                         xlab_text = "",
                                         ylab_text = "Percentage of total nuclei (%)")



my_colors <- c("#084594","#2171B5","#6BAED6", "#718FAB", "#343434", "#C6DBEF")


# mean and sd by group
get_stats_txt(fraction_nuclei_data, depend_var = "Fraction", indep_var1 = "Set", indep_var2 = "Group")


# Combine
destination <- file.path(s_OUT_dir, paste0("Plots/FANS_Nuclei_Group3.svg"))
svglite(file=destination)

ggarrange(p_Nuclei_vs_Group, p_Fraction_vs_Group,
          labels = c("B", "c"),
          ncol = 1, nrow=2)
dev.off()


destination <- file.path(s_OUT_dir, paste0("Plots/FANS_Nuclei_Group4.svg"))
svglite(file=destination)

ggarrange(p_Nuclei_vs_Group, p_Nuclei_vs_PMI, p_Fraction_vs_Group, p_Nuclei_vs_sampledate,
          labels = c("B", "D", "C", "E"),
          ncol = 2, nrow=2)
dev.off()


destination <- file.path(s_OUT_dir, paste0("Plots/FANS_Nuclei_Group_bar3.svg"))
svglite(file=destination)

ggarrange(p_Nuclei_vs_Group_bar, p_Nuclei_vs_PMI, p_Fraction_vs_Group_bar, p_Nuclei_vs_sampledate,
          labels = c("B", "D", "C", "E"),
          ncol = 2, nrow=2)
dev.off()


#-----------------------------------------------------------------------------------------------------#
#					7. DNA yield vs nuclei vs group   (nanodrop)                                       ----
#-----------------------------------------------------------------------------------------------------#

# Create data with DNA yield vs nuclei yield
total_nuclei_dna <- list(data.frame("Nuclei" = as.numeric(fans_ext$NeuN_nuclei_DNA), 
                                    "DNA" = as.numeric(fans_ext$Yield_NeuN), 
                                    "Type" = "NeuN+", "Group" = fans_ext$Group),
                          data.frame("Nuclei" = as.numeric(fans_ext$Olig2_nuclei_DNA), 
                                     "DNA" = as.numeric(fans_ext$Yield_Olig2), 
                                     "Type" = "Olig2+", "Group" = fans_ext$Group),
                          data.frame("Nuclei" = as.numeric(fans_ext$Total_rest_nuclei), 
                                     "DNA" = as.numeric(fans_ext$Yield_Rest), 
                                     "Type" = "NeuN-/Olig2-", "Group" = fans_ext$Group))

names(total_nuclei_dna) <- c("NeuN+", "Olig2+", "NeuN-/Olig2-")
total_nuclei_dna <- do.call("rbind", total_nuclei_dna)

# Factorize nuclei type variable
total_nuclei_dna$Type <- factor(total_nuclei_dna$Type, levels = c("NeuN+", "Olig2+", "NeuN-/Olig2-"))

# Plot the correlation between DNA yield per nuclei type vs Group
p_DNA_vs_nuclei_group <- scatterplot_pheno(total_nuclei_dna,
                                   x_var = "Nuclei",
                                   y_var="DNA",
                                   group_var = "Group",
                                   label_text = c("CTL", "MCI", "AD"),
                                   xlab_text ="Total isolated nuclei",
                                   ylab_text = "Total DNA yield (ng)",
                                   facet_var = "Group",
                                   facet_var2 = "Type")

cor_stats_txt(total_nuclei_dna, x_var = "Nuclei",
              y_var="DNA", 
              group_var = "Group",
              group_var2 = "Type")

# Plot the correlation between DNA yield per nuclei type
p_DNA_vs_nuclei_type <- scatterplot_pheno(total_nuclei_dna, 
                                     x_var = "Nuclei",
                                     y_var="DNA", 
                                     group_var = "Type", 
                                     label_text = c("NeuN+", "Olig2+", "NeuN-/Olig2-"), 
                                     xlab_text ="Total isolated nuclei",
                                     ylab_text = "Total DNA yield (ng)",
                                     facet_var = "Type",
                                     facet_var2 = NULL)

# DNA vs nuclei Group (including total)
total_nuclei_dna <- list(data.frame("Nuclei" = as.numeric(fans_ext$Total_nuclei), 
                                    "DNA" = as.numeric(fans_ext$total_dna),
                                    "Type" = "Total", "Group" = fans_ext$Group), 
                         data.frame("Nuclei" = as.numeric(fans_ext$NeuN_nuclei_DNA), 
                                    "DNA" = as.numeric(fans_ext$Yield_NeuN), 
                                    "Type" = "NeuN+", "Group" = fans_ext$Group),
                         data.frame("Nuclei" = as.numeric(fans_ext$Olig2_nuclei_DNA), 
                                    "DNA" = as.numeric(fans_ext$Yield_Olig2), 
                                    "Type" = "Olig2+", "Group" = fans_ext$Group),
                         data.frame("Nuclei" = as.numeric(fans_ext$Rest_nuclei_DNA), 
                                    "DNA" = as.numeric(fans_ext$Yield_Rest), 
                                    "Type" = "NeuN-/Olig2-", "Group" = fans_ext$Group))

names(total_nuclei_dna) <- c("Total", "NeuN+", "Olig2+", "NeuN-/Olig2-")

total_nuclei_dna <- do.call("rbind", total_nuclei_dna)

total_nuclei_dna$Type <- factor(total_nuclei_dna$Type, levels = c("Total", "NeuN+", "Olig2+", "NeuN-/Olig2-"))



p_DNA_vs_type_Group <- boxplot_pheno(total_nuclei_dna, 
                                   x_var = "Type",
                                   y_var="DNA", 
                                   group_var = "Group", 
                                   label_text = c("CTL", "MCI", "AD"), 
                                   xlab_text = "",
                                   ylab_text = "Total DNA yield (ng)")

get_stats_txt(total_nuclei_dna, "DNA", "Type", "Group")


# DNA vs nuclei Group (excluding total)
total_nuclei_dna_nototal <- total_nuclei_dna[total_nuclei_dna$Type!="Total",]
colnames(total_nuclei_dna_nototal) <- c("Nuclei", "DNA", "Type2", "Group")


p_DNA_vs_type_Group_nototal <- boxplot_pheno(total_nuclei_dna_nototal, 
                                     x_var = "Type2",
                                     y_var="DNA", 
                                     group_var = "Group", 
                                     label_text = c("CTL", "MCI", "AD"), 
                                     xlab_text = "",
                                     ylab_text = "Total DNA yield (ng)")

get_stats_txt(total_nuclei_dna, "DNA", "Type", "Group")

library(devtools)
library(nlcor)

# x <- data.frame("NeuN_nuclei_DNA" = fans_ext$NeuN_nuclei_DNA, "DNA_input" = fans_ext$DNA_input, "Group"= fans_ext$Group)
# x <- x[!is.na(x$NeuN_nuclei_DNA) & !is.na(x$DNA_input),]
# 
# 
# nlcor(x$NeuN_nuclei_DNA[x$Group=="AD"], x$DNA_input[x$Group=="AD"])
# 
# start_values <- c(a=0, b=0)
# 
# fit <- nls(fans_ext$DNA_input ~ a * exp(b * fans_ext$NeuN_nuclei_DNA),
#            start = start_values,
#            algorithm = "port",
#            control = nls.control(maxiter = 1000))


#-----------------------------------------------------------------------------------------------------#
#					7. DNA yield vs nuclei vs group   (Qubit)                                       ----
#-----------------------------------------------------------------------------------------------------#

## QUBITTTTT
# Plot the correlation between DNA yield per nuclei type
# p_DNA_vs_nuclei_NeuN <- scatterplot_pheno(fans_ext, 
#                                           x_var = "NeuN_nuclei_DNA",
#                                           y_var="DNA_input", 
#                                           group_var = "Group", 
#                                           label_text = c("CTL", "MCI", "AD"), 
#                                           xlab_text ="Number of NeuN+ nuclei",
#                                           ylab_text = "Post-fragmentation DNA yield (ng)",
#                                           facet_var = "Group",
#                                           facet_var2 = NULL)
# 
# # Correlation stats
# cor_stats_txt(fans_ext, x_var = "NeuN_nuclei_DNA",
#               y_var="DNA_input", 
#               group_var = "Group",
#               group_var2 = NULL)

fans_ext<- left_join(fans_ext, fans_dna_olig[,c("CaseID", "Mean_frag_size", "Total_DNA_Olig")], by = "CaseID")

# p_DNA_vs_nuclei_Olig <- scatterplot_pheno(fans_ext, 
#                                           x_var = "Olig2_nuclei_DNA",
#                                           y_var="Total_DNA_Olig", 
#                                           group_var = "Group", 
#                                           label_text = c("CTL", "MCI", "AD"), 
#                                           xlab_text ="Number of Olig2+ nuclei",
#                                           ylab_text = "Post-fragmentation DNA yield (ng)",
#                                           facet_var = "Group",
#                                           facet_var2 = NULL)

# cor_stats_txt(total_nuclei_dna, x_var = "Nuclei",
#               y_var="DNA", 
#               group_var = "Group",
#               group_var2 = "Type")

# Create data with DNA yield vs nuclei yield (Qubit)
total_nuclei_dna_qubit <- list(data.frame("Nuclei" = as.numeric(fans_ext$NeuN_nuclei_DNA), 
                                    "DNA_post" = as.numeric(fans_ext$DNA_input), 
                                    "Type" = "NeuN+", "Group" = fans_ext$Group),
                         data.frame("Nuclei" = as.numeric(fans_ext$Olig2_nuclei_DNA), 
                                    "DNA_post" = as.numeric(fans_ext$Total_DNA_Olig), 
                                    "Type" = "Olig2+", "Group" = fans_ext$Group))

names(total_nuclei_dna_qubit) <- c("NeuN+", "Olig2+")
total_nuclei_dna_qubit <- do.call("rbind", total_nuclei_dna_qubit)

# Factorize nuclei type variable
total_nuclei_dna_qubit$Type <- factor(total_nuclei_dna_qubit$Type, levels = c("NeuN+", "Olig2+"))

# Plot the correlation between DNA yield per nuclei type vs Group
p_DNA_vs_nuclei_group_post <- scatterplot_pheno(total_nuclei_dna_qubit,
                                           x_var = "Nuclei",
                                           y_var="DNA_post",
                                           group_var = "Group",
                                           label_text = c("CTL", "MCI", "AD"),
                                           xlab_text ="Total isolated nuclei",
                                           ylab_text = "Total DNA yield (ng) post-fragmentation",
                                           facet_var = "Type",
                                           facet_var2 = "Group")

my_colors <- c("#2171B5","#6BAED6", "#718FAB", "#343434", "#C6DBEF")

# Plot the correlation between DNA yield per nuclei type vs Group
p_DNA_vs_nuclei_group_post <- scatterplot_pheno(total_nuclei_dna_qubit,
                                                x_var = "Nuclei",
                                                y_var="DNA_post",
                                                group_var = "Type",
                                                label_text = c("NeuN+", "Olig2+"),
                                                xlab_text ="Total isolated nuclei",
                                                ylab_text = "Total DNA yield (ng) post-fragmentation",
                                                facet_var = "Type",
                                                facet_var2 = "Group")

cor_stats_txt(total_nuclei_dna_qubit, x_var = "Nuclei",
              y_var="DNA_post",
              group_var = "Type",
              group_var2 = "Group")


p_DNA_vs_type_Group_post <- boxplot_pheno(total_nuclei_dna_qubit, 
                                     x_var = "Type",
                                     y_var="DNA_post", 
                                     group_var = "Group", 
                                     label_text = c("CTL", "MCI", "AD"), 
                                     xlab_text = "",
                                     ylab_text = "Total DNA yield (ng) post-fragmentation")


p_DNA_vs_type_Group_post_bar <- barplot_pheno(total_nuclei_dna_qubit, 
                                          x_var = "Type",
                                          y_var="DNA_post", 
                                          group_var = "Group", 
                                          label_text = c("CTL", "MCI", "AD"), 
                                          xlab_text = "",
                                          ylab_text = "Total DNA yield (ng) post-fragmentation")

my_colors <- c("#2171B5","#6BAED6", "#718FAB", "#343434", "#C6DBEF")

p_DNA_vs_type_Group_post_bar <- barplot_pheno(total_nuclei_dna_qubit, 
                                              x_var = "Group",
                                              y_var="DNA_post", 
                                              group_var = "Type", 
                                              label_text = c("NeuN+", "Olig2+"), 
                                              xlab_text = "",
                                              ylab_text = "Total DNA yield (ng) post-fragmentation")

get_stats_txt(total_nuclei_dna_qubit, "DNA_post", "Type", "Group")

# # Combine plots
# destination <- file.path(s_OUT_dir, paste0("Plots/FANS_DNA_Group5.svg"))
# svglite(file=destination)
# # 
# # top_row <- p_DNA_vs_type_Group + p_DNA_vs_nuclei_NeuN + plot_layout(ncol =2, widths = c(1.5, 2))
# # 
# # bottom_row <- p_DNA_vs_nuclei_group + plot_layout(ncol = 1, widths = c(1, 1))
# 
# # layout <- top_row / bottom_row + plot_layout(heights= c(1.5,2)) + 
# #   plot_annotation(tag_levels = 'A')
# # layout
# 
# layout <- p_DNA_vs_type_Group_post / p_DNA_vs_nuclei_group_post + plot_layout(heights= c(1.5,2)) + 
#   plot_annotation(tag_levels = 'A')
# layout
# 
# 
# dev.off()




# bar plot version
# Combine plots
destination <- file.path(s_OUT_dir, paste0("Plots/FANS_DNA_Group_bar3.svg"))
svglite(file=destination)
# 
# top_row <- p_DNA_vs_type_Group + p_DNA_vs_nuclei_NeuN + plot_layout(ncol =2, widths = c(1.5, 2))
# 
# bottom_row <- p_DNA_vs_nuclei_group + plot_layout(ncol = 1, widths = c(1, 1))

# layout <- top_row / bottom_row + plot_layout(heights= c(1.5,2)) + 
#   plot_annotation(tag_levels = 'A')
# layout

layout <- p_DNA_vs_nuclei_group_post / p_DNA_vs_type_Group_post_bar  + plot_layout(heights= c(2, 1.5)) + 
  plot_annotation(tag_levels = 'A')
layout


dev.off()

# Free up memory
# rm(p_Fraction_vs_Group,fraction_nuclei_data, p_Nuclei_vs_Group, p_Nuclei_vs_Group_set, p_Nuclei_vs_PMI, p_Nuclei_vs_sampledate, total_nuclei_data)


#-----------------------------------------------------------------------------------------------------#
#					7. DNA integrity                                 ----
#-----------------------------------------------------------------------------------------------------#


# Create data with DNA yield vs nuclei yield (Qubit)
# total_nuclei_dna_qubit <- list(data.frame("DIN" = as.numeric(fans_ext$Mean_frag_size), 
#                                           "DNA_post" = as.numeric(fans_ext$DNA_input), 
#                                           "Type" = "NeuN+", "Group" = fans_ext$Group),
#                                data.frame("Nuclei" = as.numeric(fans_ext$Olig2_nuclei_DNA), 
#                                           "DNA_post" = as.numeric(fans_ext$Total_DNA_Olig), 
#                                           "Type" = "Olig2+", "Group" = fans_ext$Group))

# FANS DIN combined


# Create data with DNA yield vs nuclei yield (Qubit)
FANS_DIN_combined <- list(data.frame("DIN" = as.numeric(fans_ext$DIN_NeuN), 
                                          "Type" = "NeuN+", "Group" = fans_ext$Group),
                               data.frame("DIN" = as.numeric(fans_DIN_olig$DIN), 
                                          "Type" = "Olig2+", "Group" = fans_DIN_olig$Group))


names(FANS_DIN_combined) <- c("NeuN+", "Olig2+")

FANS_DIN_combined <- do.call("rbind", FANS_DIN_combined)

FANS_DIN_combined$Type <- factor(FANS_DIN_combined$Type, levels = c("NeuN+", "Olig2+"))
FANS_DIN_combined$Group <- factor(FANS_DIN_combined$Group, levels = c("CTL", "MCI", "AD"))

FANS_DIN_combined <- FANS_DIN_combined[complete.cases(FANS_DIN_combined),]


p_DIN_vs_Group <- boxplot_pheno(FANS_DIN_combined, 
                                     x_var = "Type",
                                     y_var="DIN", 
                                     group_var = "Group", 
                                     label_text = c("CTL", "MCI", "AD"), 
                                     xlab_text = "",
                                     ylab_text = "DNA integrity (DIN)")


p_DIN_vs_Group_bar <- barplot_pheno(FANS_DIN_combined, 
                                x_var = "Type",
                                y_var="DIN", 
                                group_var = "Group", 
                                label_text = c("CTL", "MCI", "AD"), 
                                xlab_text = "",
                                ylab_text = "DNA integrity (DIN)")


# DIN vs Group
fans_DIN_olig$Group <- factor(fans_DIN_olig$Group, levels = c("CTL", "MCI", "AD"))

fans_DIN_olig <- fans_DIN_olig[complete.cases(fans_DIN_olig),]

p_DIN_vs_Group_Olig <- boxplot_pheno(fans_DIN_olig, 
                                     x_var = "Group",
                                     y_var="DIN", 
                                     group_var = "Group", 
                                     label_text = c("CTL", "MCI", "AD"), 
                                     xlab_text = "",
                                     ylab_text = "DNA integrity (DIN)")

p_DIN_vs_Group_Olig_bar <- barplot_pheno(fans_DIN_olig, 
                                     x_var = "Group",
                                     y_var="DIN", 
                                     group_var = "Group", 
                                     label_text = c("CTL", "MCI", "AD"), 
                                     xlab_text = "",
                                     ylab_text = "DNA integrity (DIN)")

get_stats_txt(FANS_DIN_combined, "DIN", "Type", "Group")


#-----------------------------------------------------------------------------------------------------#
#					8. Library prep plots                                       ----
#-----------------------------------------------------------------------------------------------------#

# Create data table with date, coverage and cohort 
sample_dates <- data.table(date, "Coverage" = pheno_subset_pass$Coverage, "Group" = as.vector(pheno_subset_pass$Group), "DNA_input" = as.numeric(pheno_subset_pass$DNA_input))
sample_dates$Group <- factor(sample_dates$Group, levels = c("CTL", "MCI", "AD"))


# DNA vs Group
p_DNA_vs_Group_NeuN <- boxplot_pheno(pheno_subset_pass, 
                                x_var = "Group",
                                y_var="DNA_input", 
                                group_var = "Group", 
                                label_text = c("CTL", "MCI", "AD"), 
                                xlab_text = "",
                                ylab_text = "Total DNA input (ng)")

# DNA vs PMI
pheno_subset_pass_PMI <- pheno_subset_pass[!pheno_subset_pass$sample_names=="sample_11-71",]

p_DNA_vs_PMI <- scatterplot_pheno(pheno_subset_pass_PMI, 
                                  x_var = "PMI",
                                  y_var="DNA_input", 
                                  group_var = "Group", 
                                  label_text = c("CTL", "MCI", "AD"), 
                                  xlab_text ="Post-mortem interval (PMI)",
                                  ylab_text = "Total DNA input (ng)",
                                  facet_var = "Group",
                                  facet_var2 = NULL)

# DNA vs year
p_DNA_vs_year <- scatterplot_pheno(sample_dates, 
                                   x_var = "date",
                                   y_var="DNA_input", 
                                   group_var = "Group", 
                                   label_text = c("CTL", "MCI", "AD"), 
                                   xlab_text ="",
                                   ylab_text = "Total DNA input (ng)",
                                   facet_var = "Group",
                                   facet_var2 = NULL)


# Combine
destination <- file.path(s_OUT_dir, paste0("Plots/Pheno_Group_DNA_year2.svg"))
svglite(file=destination)

ggarrange(p_DNA_vs_Group_NeuN, p_DNA_vs_PMI, p_DNA_vs_year,
          labels = c("A", "B", "C"),
          ncol = 2, nrow=2)
dev.off()

# Free up memory
rm(p_DNA_vs_Group, p_DNA_vs_PMI, p_DNA_vs_year, p_DNA_vs_Group_NeuN, p_DNA_vs_nuclei_group,
   p_DNA_vs_nuclei_NeuN, p_DNA_vs_nuclei_type, p_DNA_vs_type_Group, pheno_subset_pass_PMI, total_nuclei_dna)





#-----------------------------------------------------------------------------------------------------#
#						10. Flowcell metadata                                                ----
#-----------------------------------------------------------------------------------------------------#

flowcell_metadata = read_excel(paste0(s_ROOT_dir, "Anno/FlowcellMetadata.xlsx"), sheet=1)

names(flowcell_metadata)



# Effect of initital pore count on data output
p_Output_vs_Pore <- scatterplot_pheno(flowcell_metadata, 
                                      x_var = "InitialPore",
                                      y_var="Data_output", 
                                      group_var = "Library", 
                                      label_text = c(">45 fmol", "≤45 fmol"), 
                                      xlab_text ="Initial pore count",
                                      ylab_text = "Total data output (Gb)",
                                      facet_var = NULL,
                                      facet_var2 = NULL)

p_Output_vs_Pore_notbatch <- scatterplot_pheno(flowcell_metadata, 
                                      x_var = "InitialPore",
                                      y_var="Data_output", 
                                      group_var = NULL, 
                                      label_text = NULL, 
                                      xlab_text ="Initial pore count",
                                      ylab_text = "Total data output (Gb)",
                                      facet_var = NULL,
                                      facet_var2 = NULL)


# Correlation stats
cor_stats_txt(flowcell_metadata, x_var = "InitialPore",
              y_var="Data_output", 
              group_var = NULL,
              group_var2 = NULL)


# Effect of number of multiplexed samples on data ouput
# p_Output_vs_Barcodes <- scatterplot_pheno(flowcell_metadata, 
#                                           x_var = "N_Barcodes",
#                                           y_var="Data_output", 
#                                           group_var = NULL,
#                                           label_text = NULL, 
#                                           xlab_text ="Number of multiplexed samples",
#                                           ylab_text = "Total data output (Gb)")

flowcell_metadata_short <- flowcell_metadata[!flowcell_metadata$N_Barcodes==10,]

flowcell_metadata_short$N_Barcodes_fact <- as.factor(flowcell_metadata_short$N_Barcodes)

flowcell_metadata$N_Barcodes_fact <- as.factor(flowcell_metadata$N_Barcodes)


# p_Output_vs_Barcodes <- scatterplot_pheno(flowcell_metadata, 
#                                           x_var = "N_Barcodes",
#                                           y_var="Data_output", 
#                                           group_var = NULL,
#                                           label_text = NULL, 
#                                           xlab_text ="Number of multiplexed samples",
#                                           ylab_text = "Total data output (Gb)")

### NOT ENOUGH OBSERVATIONS

p_Output_vs_Barcodes <- boxplot_pheno(flowcell_metadata_short, 
                                      x_var = "N_Barcodes_fact",
                                      y_var="Data_output", 
                                      group_var = "N_Barcodes_fact",
                                      label_text = c("7 samples", "8 samples", "9 samples"), 
                                      xlab_text ="Number of multiplexed samples",
                                      ylab_text = "Total data output (Gb)")


p_Output_vs_Barcodes_long <- boxplot_pheno(flowcell_metadata, 
                                      x_var = "N_Barcodes_fact",
                                      y_var="Data_output", 
                                      group_var = "N_Barcodes_fact",
                                      label_text = c("7 samples", "8 samples", "9 samples", "10 samples"), 
                                      xlab_text ="Number of multiplexed samples",
                                      ylab_text = "Total data output (Gb)")


p_Output_vs_Barcodes_bar <- barplot_pheno(flowcell_metadata_short, 
                                      x_var = "N_Barcodes_fact",
                                      y_var="Data_output", 
                                      group_var = "N_Barcodes_fact",
                                      label_text = c("7 samples", "8 samples", "9 samples"), 
                                      xlab_text ="Number of multiplexed samples",
                                      ylab_text = "Total data output (Gb)")

# Correlation stats
# get_stats_txt(flowcell_metadata, "Data_output", "N_Barcodes_fact", indep_var2 = NULL)

# Kruskal Wallis
kruskal.test(Data_output ~ N_Barcodes_fact, data = flowcell_metadata)
pairwise.wilcox.test(flowcell_metadata$Data_output, flowcell_metadata$N_Barcodes_fact,
                     p.adjust.method = "BH")




# Association between washing step and data oupput
p_Output_vs_Washing <- scatterplot_pheno(flowcell_metadata, 
                                         x_var = "Time",
                                         y_var="Data_output", 
                                         group_var = "Batch",
                                         label_text = c("Batch 1", "Batch 2"), 
                                         xlab_text ="Sequencing time (h) before flowcell washing",
                                         ylab_text = "Total data output (Gb)",
                                         facet_var = NULL,
                                         facet_var2 = NULL)

# Correlation stats
cor_stats_txt(flowcell_metadata, x_var = "Time",
              y_var="Data_output", 
              group_var = NULL,
              group_var2 = NULL)

# Association between the fmol library input and data ouput
p_Output_vs_fmol <- scatterplot_pheno(flowcell_metadata, 
                                      x_var = "fmol",
                                      y_var="Data_output", 
                                      group_var = "Library",
                                      label_text = c(">45 fmol", "≤45 fmol"), 
                                      xlab_text ="Total fmol input",
                                      ylab_text = "Total data output (Gb)",
                                      facet_var = NULL,
                                      facet_var2 = NULL)


# Association between the fmol library input and data ouput
p_Output_vs_fmol_nobatch <- scatterplot_pheno(flowcell_metadata, 
                                      x_var = "fmol",
                                      y_var="Data_output", 
                                      group_var = NULL,
                                      label_text = NULL, 
                                      xlab_text ="Total fmol input",
                                      ylab_text = "Total data output (Gb)",
                                      facet_var = NULL,
                                      facet_var2 = NULL)

# Correlation stats
cor_stats_txt(flowcell_metadata, x_var = "fmol",
              y_var="Data_output", 
              group_var = NULL,
              group_var2 = NULL)

# Coverage vs. number of multiplxed 
cov_flowcell_data <- left_join(pheno[,c("CaseID","Coverage", "Flowcell")], flowcell_metadata, by = "Flowcell", relationship = "many-to-many")

# Omit NA values
cov_flowcell_data <- cov_flowcell_data[!is.na(cov_flowcell_data$N_Barcodes_fact),]

p_Coverage_vs_Barcodes <- boxplot_pheno(cov_flowcell_data, 
                                        x_var = "N_Barcodes_fact",
                                        y_var="Coverage", 
                                        group_var = "N_Barcodes_fact",
                                        label_text = c("7 samples", "8 samples", "9 samples", "10 samples"), 
                                        xlab_text ="Number of multiplexed samples",
                                        ylab_text = "Mean sequencing depth (X) in target region")

p_Coverage_vs_Barcodes_bar <- barplot_pheno(cov_flowcell_data, 
                                        x_var = "N_Barcodes_fact",
                                        y_var="Coverage", 
                                        group_var = "N_Barcodes_fact",
                                        label_text = c("7 samples", "8 samples", "9 samples", "10 samples"), 
                                        xlab_text ="Number of multiplexed samples",
                                        ylab_text = "Mean sequencing depth (X) in target region")



get_stats_txt(cov_flowcell_data, "Coverage", "N_Barcodes_fact", indep_var2 = NULL)





# Kruskal Wallis
kruskal.test(Data_output ~ N_Barcodes_fact, data = cov_flowcell_data)
pairwise.wilcox.test(cov_flowcell_data$Data_output, cov_flowcell_data$N_Barcodes_fact,
                     p.adjust.method = "BH")



#### Data output per sample

flowcell_metadata_sample <- read_excel(paste0(s_ROOT_dir, "Anno/FlowcellMetadata.xlsx"), sheet =2)
flowcell_metadata_sample$TotalMb <- as.numeric(flowcell_metadata_sample$TotalMb)




# remove the pilots 
flowcell_metadata_sample<-flowcell_metadata_sample[flowcell_metadata_sample$N_barcodes != 5 & flowcell_metadata_sample$N_barcodes != 14,]
flowcell_metadata_sample$N_Barcodes_fact <- as.factor(flowcell_metadata_sample$N_barcodes)

get_stats_txt(flowcell_metadata_sample, "TotalMb", "N_Barcodes_fact", indep_var2= NULL)

# Kruskal Wallis
kruskal.test(TotalMb ~ N_Barcodes_fact, data = flowcell_metadata_sample)
# pairwise.wilcox.test(cov_flowcell_data$Data_output, cov_flowcell_data$N_Barcodes_fact,
#                      p.adjust.method = "BH")


p_TotalMb_vs_Barcodes <- boxplot_pheno(flowcell_metadata_sample, 
                                        x_var = "N_Barcodes_fact",
                                        y_var="TotalMb", 
                                        group_var = "N_Barcodes_fact",
                                        label_text = c("7 samples", "8 samples", "9 samples", "10 samples"), 
                                        xlab_text ="Number of multiplexed samples",
                                        ylab_text = "Total data ouput (Mb)")

p_TotalMb_vs_Barcodes_bar <- barplot_pheno(flowcell_metadata_sample, 
                                       x_var = "N_Barcodes_fact",
                                       y_var="TotalMb", 
                                       group_var = "N_Barcodes_fact",
                                       label_text = c("7 samples", "8 samples", "9 samples", "10 samples"), 
                                       xlab_text ="Number of multiplexed samples",
                                       ylab_text = "Total data ouput (Mb)")




p_CovDNA_vs_Barcodes <- boxplot_pheno(flowcell_metadata_sample, 
                                       x_var = "N_Barcodes_fact",
                                       y_var="Coverage_per_DNA", 
                                       group_var = "N_Barcodes_fact",
                                       label_text = c("7 samples", "8 samples", "9 samples", "10 samples"), 
                                       xlab_text ="Number of multiplexed samples",
                                       ylab_text = "Coverage per ng DNA")


# Combine
destination <- file.path(s_OUT_dir, paste0("Plots/Pheno_Flowcell_Metadata4.svg"))
svglite(file=destination)

ggarrange(p_Output_vs_Pore, p_Output_vs_Barcodes, p_Output_vs_fmol, p_Coverage_vs_Barcodes,
          labels = c("A", "C", "B", "D"),
          ncol = 2, nrow=2)


dev.off()


# Combine
destination <- file.path(s_OUT_dir, paste0("Plots/Pheno_Flowcell_Metadata6.svg"))
svglite(file=destination, width = 15, height = 11)


emptyplot <- ggplot() + theme_void()

ggarrange(plotlist = list(p_Output_vs_Pore, p_Output_vs_fmol, p_Output_vs_Barcodes_long, 
          p_TotalMb_vs_Barcodes_bar, emptyplot,emptyplot,
          p_Coverage_vs_Barcodes_bar, emptyplot, emptyplot), nrow = 3, ncol = 3,
          labels = c("A", "B", "C", "D", "", "", "E","", "" ), align = "v") 



dev.off()


# Combine
destination <- file.path(s_OUT_dir, paste0("Plots/Pheno_Flowcell_Metadata6_nobatch.svg"))
svglite(file=destination, width = 15, height = 11)


emptyplot <- ggplot() + theme_void()

ggarrange(plotlist = list(p_Output_vs_Pore_notbatch, p_Output_vs_fmol_nobatch, p_Output_vs_Barcodes_long, 
                          p_TotalMb_vs_Barcodes_bar, emptyplot,emptyplot,
                          p_Coverage_vs_Barcodes_bar, emptyplot, emptyplot), nrow = 3, ncol = 3,
          labels = c("A", "B", "C", "D", "", "", "E","", "" ), align = "v") 



dev.off()
# rm(p_Output_vs_Pore, p_Output_vs_fmol, p_Output_vs_Barcodes, p_Output_vs_Washing, 
#    flowcell_metadata, p_Coverage_vs_Barcodes, cov_flowcell_data)


#-----------------------------------------------------------------------------------------------------#
#						11. Scores - problem regions                                              ----
#-----------------------------------------------------------------------------------------------------#

coverage <- readRDS(paste0(s_OUT_dir, "Bedmethyl_all/bedfiles_scores_all.rds"))
gc()

# Calculate median coverage per CpG site
cov_dens <- apply(coverage, 1, median)

# Create data frame of median coverage
cov_dens_df <- as.data.frame(cov_dens)

# 7164108 positions 
# 10%, 1%, 0.1%

# Top positions - show median coverage + confidence interval 
top_pos_id <- cov_dens_df$cov_dens >= 32

high_cov_cps <- rownames(cov_dens_df)[top_pos_id]

coverage_high <- coverage[high_cov_cps,]

# Free up memory
rm(coverage)
gc()

# Plot densities coverage
destination <- file.path(s_OUT_dir, paste0("Plots/QC_CoverageDistribution_density2.svg"))
svglite(file=destination)

p_cov_distribution <- ggplot(cov_dens_df, aes(x=cov_dens)) +
  geom_density(fill="#08519C", size = 1)+
  geom_vline(aes(xintercept=quantile(cov_dens, probs = c(0.99))), color="magenta",
             linetype="dashed", size = 1.2)+
  geom_vline(aes(xintercept=quantile(cov_dens, probs = c(0.999))), color="green",
             linetype="dashed", size = 1.2)+
  geom_vline(aes(xintercept=quantile(cov_dens, probs = c(0.9999))), color="cyan",
             linetype="dashed", size = 1.2)+
  labs(x="Sequencing depth (X)", y = "Density")+
  theme_gray()

print(p_cov_distribution)

dev.off()



# Calculate stats coverages
cov_median <- apply(coverage_high, 1, median)
cov_min <- apply(coverage_high, 1, min)
cov_max <- apply(coverage_high, 1, max)

# coverage_high_stats <- data.frame(cbind(cov_median, cov_min, cov_max))

cov_40 <- apply(coverage_high, 1, function(x) quantile(x, probs = 0.40))
cov_60 <- apply(coverage_high, 1, function(x) quantile(x, probs = 0.60))

coverage_high_stats <- data.frame(cbind(cov_median, cov_40, cov_60))

# Obtain positions
position <- rownames(coverage_high_stats)
position <- data.frame("chr" = sub("chr", "", str_split_fixed(position, "_", 2)[,1]),
                       "start" = as.numeric(str_split_fixed(position, "_", 2)[,2]))

coverage_high_stats <- cbind(coverage_high_stats, position)
coverage_high_stats$chr<- factor(coverage_high_stats$chr, levels = unique(position$chr))

# Order and group regions
coverage_high_stats_ord <- coverage_high_stats %>%
  arrange(chr, start) %>%
  mutate(row_id = row_number())  # or use as factor for discrete x-axis


# coverage_high_stats_ord %>% 
#   ggplot(aes(x = row_id, y = cov_median, group =1 )) +
#   geom_ribbon(aes(ymin = cov_min, ymax = cov_max), color = "#C6DBEF", alpha = 0.3) +
#   geom_line(color = "#08519C", linewidth = 0.7) +
#   labs(
#     x = "Position",
#     y = "Sequencing deoth (X)") +
#   theme_bw() +
#   theme(legend.position = "top")

# Example graph
coverage_high_stats_ord %>% 
  ggplot(aes(x = row_id, y = cov_median, group =1)) +
  geom_ribbon(aes(ymin = cov_40, ymax = cov_60), fill = "#6BAED6", alpha = 0.5) +
  geom_line(color = "#08519C", linewidth = 0.7) +
  theme_bw() +
  labs(
    x = "Position",
    y = "Sequencing deoth (X)") +
  theme(legend.position = "top")


# Arrange and compute blocks (regions)
df_blocks <- coverage_high_stats_ord %>%
  group_by(chr) %>%
  mutate(gap = c(0, diff(start)) > 500,
         block = cumsum(gap)) %>%
  group_by(chr, block) %>%
  summarise(
    block_start = data.table::first(start),
    block_end = data.table::last(start),
    block_length = data.table::last(start) - data.table::first(start),
    cov_40 = mean(cov_40),
    median = mean(cov_median),
    cov_60 = mean(cov_60),
    CG = length(cov_median),
    .groups = "drop"
  )

print(df_blocks)

# chr   block block_start block_end block_length
# <fct> <int>       <dbl>     <dbl>        <dbl>
#   1 7         0    56371558  56372856         1298
# 2 8         0    57208589  57210471         1882
# 3 8         1    85658942  85660227         1285
# 4 17        0    43280217  43287969         7752
# 5 18        0    47022674  47027934         5260
# 6 21        0     8215695   8227903        12208
# 7 21        1     8232323   8234927         2604
# 8 Y         0    11312483  11313216          733
# 9 Y         1    11313795  11317063         3268


df_blocks_gap <- coverage_high_stats_ord %>%
  arrange(chr, start) %>%
  mutate(
    is_new_block = c(TRUE, abs(diff(start)) > 500),     # Start of new block
    block = cumsum(is_new_block)                # Block number
  ) %>%
  ungroup()


# Plot the regions
destination <- file.path(s_OUT_dir, paste0("Plots/QC_CoverageDistribution_perregion.svg"))
svglite(file=destination)

p_highcov_region <- df_blocks_gap %>% 
  ggplot(aes(x = start, y = cov_median, group =1)) +
  geom_ribbon(aes(ymin = cov_40, ymax = cov_60), fill = "#6BAED6", alpha = 0.3) +
  geom_line(color = "#08519C", linewidth = 0.8) +
  theme_bw() +
  labs(
    x = "Position",
    y = "Sequencing deoth (X)") +
  theme(legend.position = "top") + 
  facet_wrap(~ block, scales = "free", 
             labeller = labeller(block = c("1" = "chr 7",
                                           "2" = "chr 8",
                                           "3" = "chr 8",
                                           "4" = "chr 17",
                                           "5" = "chr 18",
                                           "6" = "chr 21",
                                           "7" = "chr 21",
                                           "8" = "chr Y",
                                           "9" = "chr Y")))

print(p_highcov_region)

dev.off()


# Combine
destination <- file.path(s_OUT_dir, paste0("Plots/Pheno_Cov_distribution_combined.svg"))
svglite(file=destination)

ggarrange(p_cov_distribution, p_highcov_region,
          labels = c("A", "B"),
          ncol = 1, nrow=2, heights = c(1,2))


dev.off()


#-----------------------------------------------------------------------------------------------------#
#						12. Filtering options                                              ----
#-----------------------------------------------------------------------------------------------------#

# Pre vs post filter for three independent ones
filter_data <- read_excel(paste0(s_ROOT_dir, "Anno/Coverage_high_flowcell.xlsx"), sheet = 1)
filter_data$Filter <- factor(filter_data$Filter, levels=c("Pre-filter", "Post-filter"))

p_Coverage_vs_filter_flowcell <- boxplot_pheno(filter_data, 
                                               x_var = "Flowcell_batch",
                                               y_var="Cov", 
                                               group_var = "Filter",
                                               label_text = c("Pre-filter", "Post-filter"), 
                                               xlab_text ="Flowcell",
                                               ylab_text = "Mean sequencing depth (X) in target region")

get_stats_txt(filter_data, "Cov", "Flowcell_batch", "Filter")

# example of sequencing depths high coverage regions 

filter_data_highregion <- read_excel(paste0(s_ROOT_dir, "Anno/Coverage_high_flowcell.xlsx"), sheet = 2)
colnames(filter_data_highregion) <- c("Prefilter", "Postfilter")

filter_data_highregion <- pivot_longer(filter_data_highregion, 
                                       cols = c(Prefilter, Postfilter), 
                                       names_to = "Filter", 
                                       values_to = "Coverage")

filter_data_highregion$Filter <- factor(filter_data_highregion$Filter, levels = c("Prefilter", "Postfilter"))

p_Coverage_vs_filter <- boxplot_pheno(filter_data_highregion, 
                                        x_var = "Filter",
                                        y_var="Coverage", 
                                        group_var = "Filter",
                                        label_text = c("Pre-filter", "Post-filter"), 
                                        xlab_text =NULL,
                                        ylab_text = "Mean sequencing depth (X) in HC region")

get_stats_txt(filter_data_highregion, "Coverage", "Filter", NULL)


# example of MAPQ and Read length high coverage regions 
filter_data_NCvsHC <- read_excel(paste0(s_ROOT_dir, "Anno/Coverage_high_flowcell.xlsx"), sheet = 3)

## Read length
filter_data_NCvsHC_readlength <- pivot_longer(filter_data_NCvsHC[,c(1,3)], 
                                       cols = c(ReadLength_normal, ReadLength_high), 
                                       names_to = "Region", 
                                       values_to = "ReadLength")

filter_data_NCvsHC_readlength$Region <- factor(filter_data_NCvsHC_readlength$Region , levels = c("ReadLength_normal", "ReadLength_high"))
filter_data_NCvsHC_readlength$ReadLength <- filter_data_NCvsHC_readlength$ReadLength / 1000

p_readlength_region <- ggplot(filter_data_NCvsHC_readlength, aes(x = ReadLength, fill = Region)) +
  geom_histogram(binwidth = 1) + 
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major =element_line(color='lightgrey'),
        panel.grid.minor = element_line(color='lightgrey')) +
  theme(legend.position = "top") +
  # scale_color_manual(name = group_var, values = my_colors, labels = label_text) +
  scale_fill_manual(name = "Region", values = my_colors[c(1,3)], labels = c("Normal coverage (0 - 99.9 perc.)", "High coverage (top 0.1 perc.)")) +
  xlab("Read length (kb)") + 
  ylab("Number of reads") + facet_wrap(~Region, nrow = 2, scales = "free_y") +
  scale_y_continuous(labels = label_comma())

p_readlength_region




## MAPQ and filter
filter_data_NCvsHC_MAPQ <- pivot_longer(filter_data_NCvsHC[,c(2,4)], 
                                              cols = c(MAPQ_normal, MAPQ_high), 
                                              names_to = "Region", 
                                              values_to = "MAPQ")

filter_data_NCvsHC_MAPQ$Region <- factor(filter_data_NCvsHC_MAPQ$Region , levels = c("MAPQ_normal", "MAPQ_high"))

p_MAPQ_region <- ggplot(filter_data_NCvsHC_MAPQ, aes(x = MAPQ, fill = Region)) +
  geom_histogram(binwidth = 1) + 
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major =element_line(color='lightgrey'),
        panel.grid.minor = element_line(color='lightgrey')) +
  theme(legend.position = "top") +
  # scale_color_manual(name = group_var, values = my_colors, labels = label_text) +
  scale_fill_manual(name = "Region", values = my_colors[c(1,3)], labels = c("Normal coverage (0 - 99.9 perc.)", "High coverage (top 0.1 perc.)")) +
  xlab("Mapping quality (MAPQ)") + 
  ylab("Number of reads") + facet_wrap(~Region, nrow = 2, scales = "free_y") + scale_y_continuous(labels = label_comma())

p_MAPQ_region


## MAPQ and readlength stats

filter_data_stats <- read_excel(paste0(s_ROOT_dir, "Anno/Coverage_high_flowcell.xlsx"), sheet = 4)

names(filter_data_stats)

filter_data_stats$Region <- factor(filter_data_stats$Region, levels = c("NormalCov", "HighCov"))

filter_data_stats$interaction <- interaction(filter_data_stats$Region,filter_data_stats$Filter )

p_ReadLength_vs_filter_Region <- boxplot_pheno(filter_data_stats, 
                                      x_var = "Flowcell",
                                      y_var="ReadLength", 
                                      group_var = "interaction",
                                      label_text = c("NC", "HC", "NC + filter", "HC+filter"), 
                                      xlab_text ="Flowcell",
                                      ylab_text = "Median read length (kb)")

p_MAPQ_vs_filter_Region <- boxplot_pheno(filter_data_stats, 
                                               x_var = "Flowcell",
                                               y_var="MAPQ", 
                                               group_var = "interaction",
                                               label_text = c("NC", "HC", "NC + filter", "HC+filter"), 
                                               xlab_text ="Flowcell",
                                               ylab_text = "Median mapping quality (MAPQ)")

get_stats_txt(filter_data_stats, "ReadLength", "Flowcell", "interaction")
get_stats_txt(filter_data_stats, "ReadLength", "Region", "Filter")

get_stats_txt(filter_data_stats, "MAPQ", "Flowcell", "interaction")
get_stats_txt(filter_data_stats, "MAPQ", "Region", "Filter")


# Combine
destination <- file.path(s_OUT_dir, paste0("Plots/FilterThresholds_readlength_MAPQ.svg"))
svglite(file=destination)

ggarrange( p_readlength_region,p_MAPQ_region,
           p_Coverage_vs_filter_flowcell, p_Coverage_vs_filter,
          p_ReadLength_vs_filter_Region,p_MAPQ_vs_filter_Region,
         
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow=3, align = "v")
dev.off()

# Free up workspace
rm( p_readlength_region,p_MAPQ_region,
    p_Coverage_vs_filter_flowcell, p_Coverage_vs_filter,
    p_ReadLength_vs_filter_Region,p_MAPQ_vs_filter_Region, filter_data_stats, 
    filter_data_highregion, filter_data_NCvsHC_readlength, filter_data_NCvsHC_MAPQ,
    filter_data_NCvsHC, filter_data)



#-----------------------------------------------------------------------------------------------------#
#						13. Effects on coverage, reads and calls                                                   ----
#-----------------------------------------------------------------------------------------------------#

# Coverage vs Group
p_Cov_vs_Group <- boxplot_pheno(pheno, 
                                x_var = "Group",
                                y_var="Coverage", 
                                group_var = "Group", 
                                label_text = c("CTL", "MCI", "AD"), 
                                xlab_text = "",
                                ylab_text = "Mean sequencing depth (X) in target region")


p_Cov_vs_Group_bar <- barplot_pheno(pheno, 
                                    x_var = "Group",
                                    y_var="Coverage", 
                                    group_var = "Group", 
                                    label_text = c("CTL", "MCI", "AD"), 
                                    xlab_text = "",
                                    ylab_text = "Mean sequencing depth (X) in target region")

get_stats_txt(pheno, "Coverage", "Group", indep_var2 = NULL)

# Cov vs DNA
p_Cov_vs_DNA <- scatterplot_pheno(pheno, 
                                  x_var = "DNA_input",
                                  y_var="Coverage", 
                                  group_var = "Group", 
                                  label_text = c("CTL", "MCI", "AD"), 
                                  xlab_text ="Total DNA input (ng)",
                                  ylab_text = "Mean sequencing depth (X) in target region",
                                  facet_var = "Group",
                                  facet_var2 = NULL)

cor_stats_txt(pheno, "Coverage", "DNA_input",  "Group", group_var2 = NULL)

# Reads vs Group
pheno$Reads_div <- pheno$Reads /1000

p_Reads_vs_Group <- boxplot_pheno(pheno, 
                                  x_var = "Group",
                                  y_var="Reads_div", 
                                  group_var = "Group", 
                                  label_text = c("CTL", "MCI", "AD"), 
                                  xlab_text ="",
                                  ylab_text = "Number of reads in target region (x 1000)")

p_Reads_vs_Group_bar<- barplot_pheno(pheno, 
                                     x_var = "Group",
                                     y_var="Reads_div", 
                                     group_var = "Group", 
                                     label_text = c("CTL", "MCI", "AD"), 
                                     xlab_text ="",
                                     ylab_text = "Number of reads in target region (x 1000)")

get_stats_txt(pheno, "Reads_div", "Group", indep_var2 = NULL)

# Number of calls per group
number_calls <- read.csv(file.path(s_ROOT_dir, s_out_folder, "QC/number_calls_5X.csv"))
number_calls <- data.frame(t(number_calls[,-1]))
colnames(number_calls) <- c("Total_calls", "Passed_calls")

# Rename sample names
number_calls <- cbind(number_calls, "CaseID" = rownames(number_calls))
number_calls$CaseID <- sub("sample_", "", sub("\\.", "-", number_calls$CaseID))

number_calls <- left_join(number_calls, fans_ext[,c("CaseID", "Group")])


# Long data frame of number calls
number_calls_long <- pivot_longer(number_calls, 
                                  cols = c(Total_calls, Passed_calls), 
                                  names_to = "Type", 
                                  values_to = "Calls")


number_calls_long$Type <- factor(number_calls_long$Type, levels = c("Total_calls", "Passed_calls"))


p_calls_vs_Group <- boxplot_pheno(number_calls_long, 
                                  x_var = "Group",
                                  y_var="Calls", 
                                  group_var = "Type", 
                                  label_text = c("Total", ">5X"), 
                                  xlab_text ="",
                                  ylab_text = "Number of CpG sites (x 1,000,000)")


p_calls_vs_Group_bar <- barplot_pheno(number_calls_long, 
                                      x_var = "Group",
                                      y_var="Calls", 
                                      group_var = "Type", 
                                      label_text = c("Total", ">5X"), 
                                      xlab_text ="",
                                      ylab_text = "Number of CpG sites (x 1,000,000)")


get_stats_txt(number_calls_long, "Calls", "Group", "Type")


# Combine plots
destination <- file.path(s_OUT_dir, paste0("Plots/Pheno_Group_DNA_cov2.svg"))
svglite(file=destination)

top_row <- p_Reads_vs_Group + p_Cov_vs_Group + plot_layout(ncol = 2)
bottom_row <- p_calls_vs_Group +  p_Cov_vs_DNA + plot_layout(ncol = 2, widths = c(1, 2.5))

layout <- top_row / bottom_row + plot_layout(heights= c(1,1)) + 
  plot_annotation(tag_levels = 'A')
layout

dev.off()

# Free up memory
rm(p_Reads_vs_Group, p_Cov_vs_Group, p_calls_vs_Group, number_calls, number_calls_long,
   p_Cov_vs_DNA, p_Cov_vs_Year, top_row, bottom_row, layout)
gc()



# Combine plots
destination <- file.path(s_OUT_dir, paste0("Plots/Pheno_Group_DNA_cov_bar2.svg"))
svglite(file=destination)

top_row <- p_Reads_vs_Group_bar + p_Cov_vs_Group_bar + plot_layout(ncol = 2)
bottom_row <- p_calls_vs_Group_bar + p_Cov_vs_DNA + plot_layout(ncol = 2, widths = c(1, 2.5))

layout <- top_row / bottom_row + plot_layout(heights= c(1,1)) + 
  plot_annotation(tag_levels = 'A')
layout

dev.off()

# Free up memory
rm(p_Reads_vs_Group, p_Cov_vs_Group, p_calls_vs_Group, number_calls, number_calls_long,
   p_Cov_vs_DNA, p_Cov_vs_Year, top_row, bottom_row, layout)
gc()


## On vs Off target
on_off_target <- read_excel(paste0(s_ROOT_dir, "Anno/Off_vs_Ontarget_Coverage.xlsx"), 
                       col_types = c("text", "text", "text", "numeric", "numeric", "numeric", "numeric",
                                     "numeric"),
                       sheet = 1)

# Add groups 
on_off_target <- left_join(on_off_target, pheno[,c("Barcode", "Flowcell", "Group")], by = c("Barcode", "Flowcell"))

# Divide by 1000
on_off_target$Count_global_div <- on_off_target$Count_global / 1000
on_off_target$Count_inside_mapped_div <- on_off_target$Count_inside_mapped / 1000
on_off_target$Count_off_target_div <- on_off_target$Count_off_target / 1000

index <- on_off_target$Filtered=="No"

# Reformat data (n = 142)
on_off_target_reads <- list(data.frame("Reads" = as.numeric(on_off_target$Count_off_target_div)[index],
                                          "Region" = "Off_target", "Group" = on_off_target$Group[index]),
                               data.frame("Reads" = as.numeric(on_off_target$Count_inside_mapped_div[index]),
                                          "Region" = "On_target", "Group" = on_off_target$Group[index]))

names(on_off_target_reads) <- c("Off_target", "On_target")

# Calculate factor on vs off target (reads)
on_off_target_reads <- on_off_target[on_off_target$Filtered=="No",]
ratio_reads_target <- on_off_target_filtered$Count_inside_mapped_div / on_off_target_filtered$Count_off_target_div
on_off_target_reads <- cbind(on_off_target_reads, ratio_reads_target)
get_stats_txt(on_off_target_reads, "ratio_reads_target", "Group", NULL)

# Calculate factor on vs off target (cov)

on_off_target_cov <- on_off_target[on_off_target$Filtered=="Yes",]
ratio_cov_target <- on_off_target_cov$On_target / on_off_target_cov$Off_target
on_off_target_cov <- cbind(on_off_target_cov, ratio_cov_target)
get_stats_txt(on_off_target_cov, "ratio_cov_target", "Group", NULL)


mean(ratio_reads_target)
mean(ratio_cov_target, na.rm=TRUE)


# Create data frame with gene and contrast
on_off_target_reads <- bind_rows(
  lapply(names(on_off_target_reads), function(Region) {
    tibble(Reads = on_off_target_reads[[Region]]$Reads, Set = Region, Group= on_off_target$Group[index])
  }))

# Remove NA
on_off_target_reads <- on_off_target_reads[!is.na(on_off_target_reads$Group),]

# Factorize nuclei type variable 
on_off_target_reads$Region <- factor(on_off_target_reads$Set, levels = c("Off_target", "On_target"))

p_on_off_target_bar_reads <- barplot_pheno(on_off_target_reads, 
                                              x_var = "Group",
                                              y_var="Reads", 
                                              group_var = "Region", 
                                              label_text = c("Off-target", "On-target"), 
                                              xlab_text = "",
                                              ylab_text = "Number of reads in target region (x 1000)")

get_stats_txt(on_off_target_reads, "Reads", "Region", "Group")





# Coverage 

index <- on_off_target$Filtered=="Yes"

# Reformat data (n = 142)
on_off_target_cov <- list(data.frame("Cov" = as.numeric(on_off_target$Off_target)[index],
                                       "Region" = "Off_target", "Group" = on_off_target$Group[index]),
                            data.frame("Cov" = as.numeric(on_off_target$On_target[index]),
                                       "Region" = "On_target", "Group" = on_off_target$Group[index]))

names(on_off_target_cov) <- c("Off_target", "On_target")





# Create data frame with gene and contrast
on_off_target_cov <- bind_rows(
  lapply(names(on_off_target_cov), function(Region) {
    tibble(Cov = on_off_target_cov[[Region]]$Cov, Set = Region, Group= on_off_target$Group[index])
  }))

# Remove NA
on_off_target_cov <- on_off_target_cov[!is.na(on_off_target_cov$Group),]

# Factorize nuclei type variable 
on_off_target_cov$Region <- factor(on_off_target_cov$Set, levels = c("Off_target", "On_target"))

p_on_off_target_bar_cov <- barplot_pheno(on_off_target_cov, 
                                     x_var = "Group",
                                     y_var="Cov", 
                                     group_var = "Region", 
                                     label_text = c("Off-target", "On-target"), 
                                     xlab_text = "",
                                     ylab_text = "Mean sequencing depth (X)")

get_stats_txt(on_off_target_cov, "Cov", "Region", "Group")



# Combine plots
destination <- file.path(s_OUT_dir, paste0("Plots/Pheno_Group_DNA_cov_bar3.svg"))
svglite(file=destination)

top_row <- p_on_off_target_bar_reads + p_on_off_target_bar_cov + plot_layout(ncol = 2)
bottom_row <- p_calls_vs_Group_bar + p_Cov_vs_DNA + plot_layout(ncol = 2, widths = c(1, 2.5))

layout <- top_row / bottom_row + plot_layout(heights= c(1,1)) + 
  plot_annotation(tag_levels = 'A')
layout

dev.off()


#-----------------------------------------------------------------------------------------------------#
#					HANDLING COVERAGE - REAL DATA                                                    ----
#-----------------------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------------------#
#					1. Paths to different project folders                                                          ----
#-----------------------------------------------------------------------------------------------------#
# Load libraries
library(VennDiagram)
library(data.table)
library(pbapply)
library(ComplexUpset)
library(tidyverse)
library(ggVennDiagram)
library(ggvenn)
library(UpSetR)
library(tibble)
library(dplyr)


# Jaccard Similarity function
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


# Specify results folders
DSS_folder <- file.path(s_ROOT_dir, "Results", "RESEARCH_GSM0172RRMS_ALL", "13062025_analysis_DSS")

Limma_weights_folder <- file.path(s_ROOT_dir, "Results", "RESEARCH_GSM0172RRMS_ALL","13062025_analysis_Limma")


#-----------------------------------------------------------------------------------------------------#
#					2. Overlap between outcomes - Upset Plot                                                           ----
#-----------------------------------------------------------------------------------------------------#

# Load DSS DM results
DM_Results_list_DSS <- readRDS(file.path(DSS_folder, "DM", "DM_Results_list_all.rds"))

# Load Limma results
DM_Results_list_Limma <- readRDS(file.path(Limma_weights_folder, "DM", "DM_Results_list_all.rds"))

# Function to rename the contrasts
rename_contrasts_DM <- function(DM_Results_list, method) {
  names(DM_Results_list) <-
    c(paste0("hydroxymethyl_AD_vs_CTL_", method),
      paste0("hydroxymethyl_AD_vs_MCI_", method),
      paste0("hydroxymethyl_MCI_vs_CTL_", method),
      paste0("methyl_AD_vs_CTL_", method),
      paste0("methyl_AD_vs_MCI_", method),
      paste0("methyl_MCI_vs_CTL_", method),
      paste0("modified_AD_vs_CTL_", method),
      paste0("modified_AD_vs_MCI_", method),
      paste0("modified_MCI_vs_CTL_", method))

  return(DM_Results_list)
}

# Rename the contrasts
DM_Results_list_Limma <- rename_contrasts_DM(DM_Results_list_Limma, "Limma")
DM_Results_list_DSS <- rename_contrasts_DM(DM_Results_list_DSS, "DSS")

# Combine DM Results of all methods
DM_Results_list_all <- do.call(c, list(DM_Results_list_Limma, DM_Results_list_DSS))

# Exclude modified data
temp = str_detect(names(DM_Results_list_all), "modified", negate = TRUE)

# Get top 1000 CpG positions
cpg_ids <- pblapply(DM_Results_list_all[temp], function(DM_Results) {
  cpg_id <-DM_Results$cpg_id[1:1000]
  return(cpg_id)
})

# Prepare table output
results <- data.frame(Vector1 = character(),
                      Vector2 = character(),
                      SharedCpGs = character(),
                      stringsAsFactors = FALSE)

# Get number of vectors
names_list <- names(cpg_ids)
n <- length(cpg_ids)

# Loop through all unique pairs
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    common_cpgs <- intersect(cpg_ids[[i]], cpg_ids[[j]])
    if (length(common_cpgs) > 0) {
      results <- rbind(
        results,
        data.frame(
          Vector1 = names_list[i],
          Vector2 = names_list[j],
          SharedCpGs = paste(common_cpgs, collapse = ", "),
          stringsAsFactors = FALSE))}}
}

print(results)

# Create data frame with cpg ID and contrast (long data frame)
df_long_region <- bind_rows(
  lapply(names(cpg_ids), function(name) {
    tibble(CpG = cpg_ids[[name]], Set = name)
  })
)

# Function to create a matrix with binary values (0 = no overlap, 1 = overlap)
create_binary_matrix <- function(df_long, name) {
  
  # Turn to long format and spread to wide binary format
  binary_matrix <- df_long %>%
    dplyr::distinct() %>%
    dplyr::mutate(present = 1) %>%
    tidyr::pivot_wider(names_from = Set, values_from = present, values_fill = 0)
  
  # Add rowsums to binary_matrix (total CpG sites that overlap)
  binary_matrix_rowsums <- binary_matrix %>%
    dplyr::mutate(overlap_count = rowSums(across(-CpG)))
  
  # Rename column names
  binary_matrix <- binary_matrix[,-1]
  colnames(binary_matrix) <- gsub("_", " ", colnames(binary_matrix))
  # binary_matrix <- binary_matrix[,order(colnames(binary_matrix))]
  
  # Create Upsetplot
  image <- ComplexUpset::upset(binary_matrix, intersect = colnames(binary_matrix),
                               set_sizes = FALSE, min_size = 2, 
                               themes=list(default=theme(axis.text.x=element_blank()))) + 
    theme_classic() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
                          
  print(image)
  
  # Save UpsetPlot
  library(ggpubr)
  destination = file.path(s_OUT_dir, "Plots", paste("DM17_UpsetPlot_",name,".svg"))
  svglite(destination, width = 20, height = 8)
  
  print(image)
  dev.off()
 
  return(binary_matrix_rowsums)
}

binary_matrix_region <- create_binary_matrix(df_long_region, "top1000cpgs_nomod")


#-----------------------------------------------------------------------------------------------------#
#					3. Overlap between outcomes - Jaccard Similarity                                                           ----
#-----------------------------------------------------------------------------------------------------#

# Jaccard similarity (0 = no similarity, 1 = perfect similarity)

# Get top 1000 CpG positions
cpg_ids <- pblapply(DM_Results_list_all, function(DM_Results) {
  cpg_id <-DM_Results$cpg_id[1:1000]
  return(cpg_id)
})

# Create matrix with all contrast pairs
jaccard_matrix <- matrix(data = NA, nrow = length(names(DM_Results_list_all)), ncol = length(names(DM_Results_list_all)))
rownames(jaccard_matrix) <- names(DM_Results_list_all)
colnames(jaccard_matrix) <- names(DM_Results_list_all)

# Calculate jaccard value and add to matrix
# Loop through all unique pairs
for (first in 1:(n)) {
  for (second in (first):n) {
    jaccard_value <- jaccard(cpg_ids[[first]], cpg_ids[[second]])
    jaccard_matrix[first,second] <- jaccard_value
    jaccard_matrix[second,first] <- jaccard_value
    
    }
}

print(jaccard_matrix)

# Heatmap of Jaccard matrix
destination = file.path(s_OUT_dir, "Plots", paste("DM17_HeatmapJaccard_blue.svg"))
svglite(destination, height = 10, width = 10)

labels = sub("methyl", "5-mC", sub("hydroxymethyl", "5-hmC", colnames(jaccard_matrix)))
labels = sub("modified", "mod", labels)
labels = gsub("_", "  ",sub("_vs_", "-", labels))
labels

coul <- colorRampPalette(brewer.pal(9, "Blues"))(30)
heatmap(jaccard_matrix, Colv = NA, Rowv = NA, scale="column", col = coul, 
        labRow = labels, labCol = labels, cexRow = 1, cexCol = 1, margins = c(10,10)) 

dev.off()


#-----------------------------------------------------------------------------------------------------#
#				4. Correlation test statistics  (effect sizes, p-values)                                                        ----
#-----------------------------------------------------------------------------------------------------#

# All contrasts 
contrasts <- c("hydroxymethyl_AD_vs_CTL",
               "hydroxymethyl_AD_vs_MCI",
               "hydroxymethyl_MCI_vs_CTL",
               "methyl_AD_vs_CTL",
               "methyl_AD_vs_MCI",
               "methyl_MCI_vs_CTL",
               "modified_AD_vs_CTL",
               "modified_AD_vs_MCI",
               "modified_MCI_vs_CTL")

# Name methods of differential methylation analysis
method_one = "LimmaW"
method_two = "DSS"

# Calculate and plot correlation statistics of all contrast pairs
corr_list <- pblapply(contrasts, function(contrast){
  
  DM_Results_method_one <- DM_Results_list_all[[paste0(contrast, "_", method_one)]]
  
  DM_Results_method_two <- DM_Results_list_all[[paste0(contrast, "_", method_two)]]
  
  compare_df <- left_join(DM_Results_method_one[, c("stat", "pvals", "cpg_id")],
                          DM_Results_method_two[, c("stat", "pvals", "cpg_id")], by = "cpg_id")
  
  head(compare_df)
  
  cat(paste0("Contrast: ", contrast, "\n"))
  
  corr_stat <- stats::cor.test(compare_df$stat.x, compare_df$stat.y)
  corr_pvals <- stats::cor.test(compare_df$pvals.x, compare_df$pvals.y)
  
  corr_df <- data.frame("stat_t" = corr_stat$statistic, "stat_p.value" = corr_stat$p.value, "stat_estimate" = corr_stat$estimate,
             "pvals_t" = corr_pvals$statistic, "pvals_p.value" = corr_pvals$p.value, "pvals_estimate" = corr_pvals$estimate)
  

  # destination = file.path(s_OUT_dir, "Plots", paste("CorrPlot_Tstats_",contrast,".png"))
  # png(destination)
  
  image <- ggplot(compare_df, aes(x = stat.x, y = stat.y)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "lm", color = "#2171B5", fill = "grey") +
    labs(
      title = paste("Correlation = ", round(corr_stat$estimate, 3), " p-value = ", corr_stat$p.value),
      x = "T-statistic using Limma",
      y =  "T-statistic using DSS") +
    theme_bw()
  
  # print(image)
  # 
  # dev.off()
  
  # Free up memory
  rm(compare_df,DM_Results_method_one, DM_Results_method_two )
  gc()
  
  return(image)
  
})

names(corr_list) <- contrasts
# corr_list <- do.call("rbind", corr_list)

# Combine the plots
destination = file.path(s_OUT_dir, "Plots", paste("CorrPlot_Tstats_combined.png"))
png(destination, width = 1080, height = 1080)

ggarrange(plotlist=corr_list, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
          ncol=3, nrow=3)

dev.off()

# ggqqplot(compare_df$pvals.x, ylab = "Limma") 


#-----------------------------------------------------------------------------------------------------#
#					5. Plot association between coverage and mod %                                                      ----
#-----------------------------------------------------------------------------------------------------#

# Subset scores and turn into vector 
scores_subset <- scores[idx_cpgs, idx_samples]
saveRDS(scores_subset, file = file.path(s_OUT_dir, "Bedmethyl_all/bedfiles_scores_all_5X_pass_nosexchr_subset.rds"))

# turn scores into vector
scores <- as.vector(scores)
scores_subset <- as.vector(scores_subset)

gc()


temp = c("modified", "methyl", "hydroxymethyl")


pblapply(temp, function(mod_type){
  
  # Load M values
  # Metas <- as.matrix(readRDS(file.path(s_OUT_dir, paste0("QC/", mod_type, "_all_", cutoff, "X_pass_Metas_nosexchr.rds"))))
  betas <- as.matrix(readRDS(file.path(s_OUT_dir, paste0("Bedmethyl_all/bedfiles_", mod_type,"_all_",cutoff, "X_pass_nosexchr.rds"))))
  gc()
  
  # Long data format
  # data <- data.table(cov = scores, mod = as.vector(betas))
  # data <- na.omit(data)
  # gc()
  # 
  # # calculate correlation
  # corr <- stats::cor(data$cov, data$mod, method = "pearson")
  # cat(paste0("Correlation all elements (N = ", nrow(data), ": ", corr))
  
  # Subset matrix and convert to vector
  mod_subset <- betas[idx_cpgs,idx_samples]
  mod_subset <- as.vector(mod_subset)
  
  rm(betas)
  gc()
  
  data_subset <- data.table(cov = scores_subset, mod = mod_subset)
  data_subset <- na.omit(data_subset)
  
  # Calculate correlation for subset 
  corr_subset <- stats::cor(data_subset$cov, data_subset$mod, method = "pearson")
  cat(paste0("Correlation subset (N = ", length(idx_cpgs), ": ", corr_subset, "\n"))
  
  # Plot correlation plot
  image <- ggplot(data_subset, aes(x = cov, y = mod)) +
    geom_point(alpha = 0.5) +
    # geom_smooth(method = "lm", color = "#2171B5") +
    labs(
      title = paste("Correlation = ", round(corr_subset, 3)),
      x = "Sequencing depth (X)",
      y =  paste0(mod_type, " %")
    ) +
    theme_bw() +
    coord_cartesian(xlim = c(0,100))
  
  cat(paste0("X - Y plot for ", mod_type, "\n"))
  
  # ggexport(image, filename= file.path(s_ROOT_dir, s_out_folder,paste0("Plots/QC_Corr_Cov_vs_",mod_type,"_subset_zoom.png")))
  
  destination = file.path(s_ROOT_dir, s_out_folder,paste0("Plots/QC_Corr_Cov_vs_",mod_type,"_subset_zoom.svg"))
  svglite(filename = destination)
  
  print(image) 
  
  dev.off()
  
  rm(image)
  
})


rm(scores)

#-----------------------------------------------------------------------------------------------------#
#					HANDLING COVERAGE - SUBSET DATA                                                       ----
#-----------------------------------------------------------------------------------------------------#

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

# # Get random samples per group 
#
# n_samples = 30 # per group
# 
# idx_samples <- pblapply(group_samples, function(group) {
#   set.seed(36)
# 
#   idx <- sample(group, size = n_samples)
#   return(idx)
# })
#
# idx_samples_CTL <- idx_samples$CTL
# idx_samples_AD <- idx_samples$AD
#
# idx_samples <- unlist(idx_samples)

rm(groups, group_samples, sample_names)

#-----------------------------------------------------------------------------------------------------#
#						2. Select CpG sites                                                       ----
#-----------------------------------------------------------------------------------------------------#
# Load scores file
scores <- as.matrix(readRDS(file = paste0(s_ROOT_dir, s_out_folder,"Bedmethyl_all/bedfiles_scores_all_5X_pass_nosexchr.rds")))


# Randomly select cpg sites
n_cpg_sites = 200000

set.seed(72)
idx_cpgs <- sample(nrow(scores), n_cpg_sites)

rm(scores)
gc()




# -----------------------------------------------------------------------------------------------------#
# 				5. Create subset QC passed data to compare DM analyses                                                          ----
# -----------------------------------------------------------------------------------------------------#

# Passed samples 5 mC and 5-hmc data
temp = list.files(file.path(s_OUT_dir, "QC/"), pattern=paste0(cutoff, "X_pass_Metas_nosexchr"),
                  full.names = TRUE)

pblapply(temp, function(file) {
  
  # Read M values
  Metas <- readRDS(file)

  # Vector with file name
  file_base_name <- sub(".rds", "", basename(file))

  # Subset M value matrix
  Metas <- Metas[idx_cpgs, idx_samples]

  # Save new matrix
  saveRDS(Metas, file = file.path(s_OUT_dir, "QC/", paste0(file_base_name, "_subset.rds")))
})


# Count and beta data
temp = list.files(file.path(s_OUT_dir, "Bedmethyl_all/"), pattern=paste0(cutoff, "X_pass_nosexchr.rds"),
                  full.names = TRUE)

# temp = temp[str_detect(temp, "count", negate = FALSE)]

pblapply(temp, function(file){

  # Read count values
  counts <- readRDS(file)

  # Extract the file base name
  file_base_name <- sub(".rds", "", basename(file))

  # Subset the count data
  counts <- counts[idx_cpgs, idx_samples]

  # Save count data
  saveRDS(counts, file = file.path(s_OUT_dir, "Bedmethyl_all/", paste0(file_base_name, "_subset.rds")))
})


#-----------------------------------------------------------------------------------------------------#
#					3. Impute differences                                                        ----
#-----------------------------------------------------------------------------------------------------#

temp = list.files(file.path(s_OUT_dir, "Bedmethyl_all/"), pattern=paste0(cutoff, "X_pass_nosexchr_subset"),
                  full.names = TRUE)

temp = str_subset(temp, "count|scores", negate=TRUE)


betas <- readRDS(temp[2])

n_samples <- ncol(betas)
n_true_cpgs <- 500  # Number of true DMPs to spike in
effect_size <- 20  # e.g., 20% methylation difference
# group <- rep(c(0, 1), each = n_samples)  # binary phenotype

# Randomly pick CpGs to spike
set.seed(476)

true_diff_cpgs <- sample(1:length(idx_cpgs), n_true_cpgs)

# Inject effect
betas_sim <- betas

for (i in true_diff_cpgs) {
  if (runif(1) > 0.5) {
    betas_sim[i, idx_samples_AD] <- pmin(betas[i, idx_samples_AD] + effect_size, 100)
  } else {
    betas_sim[i, idx_samples_AD] <- pmax(betas[i, idx_samples_AD] - effect_size, 0)
  }
}


# x <- rbind(betas_sim[true_diff_cpgs[1],], betas[true_diff_cpgs[1],])


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

# Convert to M values
Metas_sim <- mod2M(as.data.frame(imputezerohundred(betas_sim)))
summary(Metas_sim)

# Convert to count values
coverage = readRDS(file.path(s_OUT_dir, "Bedmethyl_all", "bedfiles_scores_all_5X_pass_nosexchr_subset.rds"))

counts_sim <- round(betas_sim * coverage / 100)
counts_sim[is.na(counts_sim)] <- 0

rm(betas)

# Save objects
saveRDS(betas_sim, file = file.path(s_OUT_dir, "Bedmethyl_all/", "bedfiles_methyl_all_5X_pass_nosexchr_subset_sim.rds"))

saveRDS(Metas_sim, file = file.path(s_OUT_dir, "QC/", "methyl_all_5X_pass_Metas_nosexchr_subset_sim.rds"))

saveRDS(counts_sim, file = file.path(s_OUT_dir, "Bedmethyl_all/", "bedfiles_methyl_count_all_5X_pass_nosexchr_subset_sim.rds"))

# Loading objects
betas_sim <- readRDS(file = file.path(s_OUT_dir, "Bedmethyl_all/", "bedfiles_methyl_all_5X_pass_nosexchr_subset_sim.rds"))
Metas_sim <- readRDS(file = file.path(s_OUT_dir, "QC/", "methyl_all_5X_pass_Metas_nosexchr_subset_sim.rds"))
counts_sim <- readRDS(file = file.path(s_OUT_dir, "QC/", "methyl_all_5X_pass_Metas_nosexchr_subset_sim.rds"))


#-----------------------------------------------------------------------------------------------------#
#					6. Load subset data to compare                                                           ----
#-----------------------------------------------------------------------------------------------------#

# temp = list.files(file.path(s_OUT_dir, "QC/"), pattern=paste0(cutoff, "X_pass_Metas_nosexchr_subset"),
#                   full.names = TRUE)

# Pheno subset for benchmarking
pheno_subset_benchmark <- pheno_subset_pass[idx_samples,]
pheno_subset_benchmark$Group <- factor(pheno_subset_benchmark$Group, levels = c("CTL", "AD"))
pheno_subset_benchmark$Gender <- factor(pheno_subset_benchmark$Gender)
pheno_subset_benchmark$Flowcell <- factor(pheno_subset_benchmark$Flowcell)


#-----------------------------------------------------------------------------------------------------#
#					7. Perform statistical analyses using different methods                                                         ----
#-----------------------------------------------------------------------------------------------------#

# Perform analysis - check run time and memory (profvis)

#-----------------------------------------------------------------------------------------------------# 
# ~ Limma regular ----

# Make design matrix formula; correcting for Subject, determining the effect of Tissue
s_CovFormulaImportantSV = ~ Group + Age + factor(Gender) + PMI + factor(Flowcell)

# make design matrix
temp_design = model.matrix(as.formula(s_CovFormulaImportantSV),data=pheno_subset_benchmark)
colnames(temp_design) = make.names(colnames(temp_design))
head(temp_design)


# Limma without weights function
DM_analysis_single_limma <- function(Metas_file, mod_type) {
  
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
  tryCatch(saveRDS(Limmaefit, file = file.path(paste0(s_ROOT_dir,s_out_folder,"DM/Limmaefit_Limma_", mod_type, ".rds"))), error = function(e) warning(e$message))
  
  loop_contrast = "GroupAD"
    s_coeff = loop_contrast 
    
    cat(paste0("Extract table of the top-ranked genes", "\n\n"))
    DM_Results = limma::topTable(Limmaefit, coef=s_coeff, num=Inf, sort.by = "P", adjust="BH") # Data of specified variables
    
    # Order based on p value
    DM_Results_sign=DM_Results[DM_Results$P.Value<0.01,]
    print(head(DM_Results,20))
    
    destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DM_hist_Limma_",mod_type, "_",loop_contrast, ".pdf")
    pdf(file=destination)
    
    # histogram of p values and adjusted p values
    hist(DM_Results$P.Value ,xlim=c(0,1),main=s_coeff)
    hist(DM_Results$adj.P.Val ,xlim=c(0,1),main=s_coeff)
    min(DM_Results$adj.P.Val)
    
    
    # Descriptives
    cat(paste0("Minimum adjusted p value of ", s_coeff, ": " ,round(min(DM_Results$fdrs),digits = 3), "\n\n"))
    
    cat(paste0("Contrast: ",s_coeff,".\n  number of features below TH (0.01)\n  Pvalue:",
               sum(DM_Results$P.Value<0.01)," (",round((sum(DM_Results$P.Value<0.01)/dim(DM_Results)[1])*100,0),"%)",
               "\n  AdjPvalue:",sum(DM_Results$adj.P.Val<0.05)," (",round((sum(DM_Results$adj.P.Val<0.05)/dim(DM_Results)[1])*100,0),"%)","\n\n"))
    
    dev.off()
    
    cat(paste0("Save results"),"\n\n")
    
    tryCatch(saveRDS(DM_Results, file = file.path(paste0(s_ROOT_dir,s_out_folder,"DM/DM_Results_Limma_", mod_type, "_",loop_contrast, ".rds"))), error = function(e) warning(e$message))
    # tryCatch(fwrite(DM_Results,file = paste0(s_ROOT_dir,s_out_folder,"DM/DM_Results_",mod_type, "_", loop_contrast,".csv")), error = function(e) warning(e$message))
    
    rm(Limmaefit)
    gc()
    
    return(DM_Results)
}


# Running analysis
limma_profile <- profvis(
  DM_analysis_single_limma(Metas_file = file.path(s_OUT_dir, "QC", "methyl_all_5X_pass_Metas_nosexchr_subset_sim.rds"),
                           mod_type = "methyl")
)

# Save profile
saveRDS(limma_profile, file.path(s_OUT_dir, "DM", "limma_profile.rds"))
htmlwidgets::saveWidget(limma_profile, file = file.path(s_OUT_dir, "DM", "limma_profile.html"))


#-----------------------------------------------------------------------------------------------------# 
# ~ Limma with weights ----

# Limma with weights function
DM_analysis_single_limmaWeights <- function(Metas_file, mod_type, scores_file, group_factor, coverage_file, pheno_file) {
  
  # Load Metas
  Metas <- readRDS(Metas_file)
  Metas <- as.matrix(Metas)
  
  # Calculate weights
  scores_weights <- weights_global_sample(scores_file, group_factor, coverage_file, pheno_file)
  scores_weights <- as.matrix(scores_weights)
  
  # Make LMmodel with weights
  cat(paste0("Fitting linear model"),"\n\n")
  temp_fit = limma::lmFit(as.matrix(Metas),temp_design, weights=scores_weights)
  
  # Show the coefs for first gene
  data.frame(temp_fit$coef[1,])
  
  cat(paste0("Rank genes using empirical Bayes method"),"\n\n")
  Limmaefit = limma::eBayes(temp_fit)
  
  # Save results limma
  tryCatch(saveRDS(Limmaefit, file = file.path(paste0(s_ROOT_dir,s_out_folder,"DM/Limmaefit_LimmaW_", mod_type, ".rds"))), error = function(e) warning(e$message))
  
  loop_contrast = "GroupAD"
  s_coeff = loop_contrast 
  
  cat(paste0("Extract table of the top-ranked genes", "\n\n"))
  DM_Results = limma::topTable(Limmaefit, coef=s_coeff, num=Inf, sort.by = "P", adjust="BH") # Data of specified variables
  
  # Order based on p value
  DM_Results_sign=DM_Results[DM_Results$P.Value<0.01,]
  print(head(DM_Results,20))
  
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DM_hist_LimmaW_",mod_type, "_",loop_contrast, ".pdf")
  pdf(file=destination)
  
  # histogram of p values and adjusted p values
  hist(DM_Results$P.Value ,xlim=c(0,1),main=s_coeff)
  hist(DM_Results$adj.P.Val ,xlim=c(0,1),main=s_coeff)
  min(DM_Results$adj.P.Val)
  
  
  # Descriptives
  cat(paste0("Minimum adjusted p value of ", s_coeff, ": " ,round(min(DM_Results$fdrs),digits = 3), "\n\n"))
  
  cat(paste0("Contrast: ",s_coeff,".\n  number of features below TH (0.01)\n  Pvalue:",
             sum(DM_Results$P.Value<0.01)," (",round((sum(DM_Results$P.Value<0.01)/dim(DM_Results)[1])*100,0),"%)",
             "\n  AdjPvalue:",sum(DM_Results$adj.P.Val<0.05)," (",round((sum(DM_Results$adj.P.Val<0.05)/dim(DM_Results)[1])*100,0),"%)","\n\n"))
  
  dev.off()
  
  cat(paste0("Save results"),"\n\n")
  
  tryCatch(saveRDS(DM_Results, file = file.path(paste0(s_ROOT_dir,s_out_folder,"DM/DM_Results_LimmaW_", mod_type, "_",loop_contrast, ".rds"))), error = function(e) warning(e$message))
  # tryCatch(fwrite(DM_Results,file = paste0(s_ROOT_dir,s_out_folder,"DM/DM_Results_",mod_type, "_", loop_contrast,".csv")), error = function(e) warning(e$message))
  
  rm(Limmaefit)
  gc()
  
  return(DM_Results)
}

# Creating weights function
weights_global_sample <- function(scores_file, group_factor, coverage_file, pheno_file) {
  
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
    group_weights <- round(sweep(group_scores, 1, rowSums(coverage, na.rm=TRUE), `/`), digits = 3)
    
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
  saveRDS(scores_weights, file = paste0(s_ROOT_dir,s_out_folder,"DM/scores_group_weights_", cutoff,"X_pass_benchmarking.rds")) 
  
  return(scores_weights)
  
  # Free up memory
  rm(scores_weights,group_scores, group_weights)
  gc()
  
}



# Running analysis
limma_weights_profile <- profvis(
  
  DM_analysis_single_limmaWeights(file.path(s_OUT_dir, "QC", "methyl_all_5X_pass_Metas_nosexchr_subset_sim.rds"),
                                  mod_type = "methyl", 
                                  scores_file = file.path(s_OUT_dir, "Bedmethyl_all", "bedfiles_scores_all_5X_pass_nosexchr_subset.rds"), 
                                  group_factor = pheno_subset_benchmark$Group, 
                                  coverage_file = file.path(s_OUT_dir, "Bedmethyl_all", "bedfiles_scores_all_5X_pass_nosexchr_subset.rds"), 
                                  pheno_file = pheno_subset_benchmark)
)

# Save Results
saveRDS(limma_weights_profile, file.path(s_OUT_dir, "DM", "limma_weights_profile.rds"))
htmlwidgets::saveWidget(limma_weights_profile, file = file.path(s_OUT_dir, "DM", "limma_weights_profile.html"))




# table(DM_Results$P.Value<1)
# 
# FALSE   TRUE 
# 473 199527 



#-----------------------------------------------------------------------------------------------------# 
# ~ Methylkit  ----

DM_analysis_single_methylKit <- function(mod_type, counts_file, coverage, treatment_factor){
  
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
  
  saveRDS(methylRaw_list_subset_sim, file = file.path(s_OUT_dir, "QC", paste0("methylRaw_list_",mod_type,".rds")))
  
  # Unite into one Base object
  cat(paste0("Uniting methylRawList", "\n\n"))
  methylBase_obj_sim <- unite(methylRaw_list_subset_sim)
  
  saveRDS(methylBase_obj_sim, file = file.path(s_OUT_dir, "QC", paste0("methylBase_obj_",mod_type,".rds")))
  
  rm(methylRaw_list_subset_sim)
  gc()

  
  # ---------------------------------------------------------------------------- # 
  
  # Modelling without correcting for overdispersion
  cat(paste0("Running model: performing calculateDiffMeth", "\n\n"))

  DM_Results = calculateDiffMeth(methylBase_obj_sim, covariates = covariates, overdispersion = "none")
  colnames(DM_Results) <- c("chr", "start", "end", "strand", "pvals", "fdrs", "meth.diff")
  
  
  # Histogram of p values and adjusted p values
  destination <- file.path(s_OUT_dir, "Plots", paste0("DM_hist_Methylkit_", mod_type, ".pdf"))
  pdf(file=destination)
  
  hist(DM_Results$pvals, xlim=c(0,1), main =s_coeff)
  hist(DM_Results$fdrs, xlim=c(0,1), main =s_coeff)
  
  dev.off()
  
  # Descriptives
  cat(paste0("Minimum adjusted p value of ", s_coeff, ": " ,round(min(DM_Results$fdrs),digits = 3), "\n\n"))
  
  cat(paste0("Contrast: ",s_coeff,".\n  number of features below TH (0.01)\n  Pvalue:",
             sum(DM_Results$pvals<0.01, na.rm =TRUE)," (",round((sum(DM_Results$pvals<0.01, na.rm=TRUE)/dim(DM_Results)[1])*100,0),"%)",
             "\n  AdjPvalue:",sum(DM_Results$fdrs<0.05, na.rm = TRUE)," (",round((sum(DM_Results$fdrs<0.05, na.rm = TRUE)/dim(DM_Results)[1])*100,0),"%)","\n\n"))
  
  # Saving DM Results
  saveRDS(DM_Results, file = file.path(s_OUT_dir, "DM", paste0("DM_Results_methylkit_",mod_type,".rds")))
}

# WITH CORRECTION
DM_analysis_single_methylKit_MNcor <- function(mod_type, counts_file, coverage, treatment_factor){
  
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
  
  saveRDS(methylRaw_list_subset_sim, file = file.path(s_OUT_dir, "QC", paste0("methylRaw_list_",mod_type,".rds")))
  
  # Unite into one Base object
  cat(paste0("Uniting methylRawList", "\n\n"))
  methylBase_obj_sim <- unite(methylRaw_list_subset_sim)
  
  saveRDS(methylBase_obj_sim, file = file.path(s_OUT_dir, "QC", paste0("methylBase_obj_",mod_type,".rds")))
  
  rm(methylRaw_list_subset_sim)
  gc()
  
  
  # ---------------------------------------------------------------------------- # 
  
  # Modelling with correction for overdispersion
  cat(paste0("Running model: performing calculateDiffMeth with MN correction", "\n\n"))
  
  DM_Results_overd = calculateDiffMeth(methylBase_obj_sim, covariates = covariates, overdispersion = "MN", test = "Chisq")
  colnames(DM_Results_overd) <- c("chr", "start", "end", "strand", "pvals", "fdrs", "meth.diff")
  
  # Histogram of p values and adjusted p values
  destination <- file.path(s_OUT_dir, "Plots", paste0("DM_hist_Methylkit_overd_", mod_type, ".pdf"))
  pdf(file=destination)
  
  hist(DM_Results_overd$pvals, xlim=c(0,1), main =s_coeff)
  hist(DM_Results_overd$fdrs, xlim=c(0,1), main =s_coeff)
  
  dev.off()
  
  # Descriptives
  cat(paste0("Minimum adjusted p value of ", s_coeff, ": " ,round(min(DM_Results_overd$fdrs),digits = 3), "\n\n"))
  
  cat(paste0("Contrast: ",s_coeff,".\n  number of features below TH (0.01)\n  Pvalue:",
             sum(DM_Results_overd$pvals<0.01, na.rm =TRUE)," (",round((sum(DM_Results_overd$pvals<0.01, na.rm=TRUE)/dim(DM_Results_overd)[1])*100,0),"%)",
             "\n  AdjPvalue:",sum(DM_Results_overd$fdrs<0.05, na.rm = TRUE)," (",round((sum(DM_Results_overd$fdrs<0.05, na.rm = TRUE)/dim(DM_Results_overd)[1])*100,0),"%)","\n\n"))
  
  saveRDS(DM_Results_overd, file = file.path(s_OUT_dir, "DM", paste0("DM_Results_methylkit_",mod_type,"_MNcorrection.rds")))
  
}


# Treatment factor coded 0 to 1
treatment_factor <- factor(pheno_subset_benchmark$Group, labels = c(0,1), levels = c("CTL", "AD"))

# Covariates
covariates <- pheno_subset_benchmark[,c("Age", "Gender", "PMI", "Flowcell")]
covariates$Flowcell<- factor(covariates$Flowcell)
covariates$Gender <- factor(covariates$Gender)

# coverage = readRDS(file.path(s_OUT_dir, "Bedmethyl_all", "bedfiles_scores_all_5X_pass_nosexchr_subset.rds"))

methylKit_profile <- profvis(
  DM_analysis_single_methylKit(mod_type = "methyl", 
                       counts_file = file.path(s_OUT_dir, "Bedmethyl_all/", paste0("bedfiles_methyl_count_all_",cutoff,"X_pass_nosexchr_subset_sim.rds")),
                       coverage = coverage,
                       treatment_factor = treatment_factor)
)

saveRDS(methylKit_profile, file.path(s_OUT_dir, "DM", "methylKit_profile.rds"))

methylKit_MNcor_profile <- profvis(
  DM_analysis_single_methylKit_MNcor(mod_type = "methyl", 
                               counts_file = file.path(s_OUT_dir, "Bedmethyl_all/", paste0("bedfiles_methyl_count_all_",cutoff,"X_pass_nosexchr_subset_sim.rds")),
                               coverage = coverage,
                               treatment_factor = treatment_factor)
)

saveRDS(methylKit_profile_MNcor, file.path(s_OUT_dir, "DM", "methylKit_MNcor_profile.rds"))

# overdispersion "MN" is too conservative and strict - almost all p values are 1 

htmlwidgets::saveWidget(methylKit_profile, file = file.path(s_OUT_dir, "DM", "methylKitprofile.html"))
htmlwidgets::saveWidget(methylKit_MNcor_profile, file = file.path(s_OUT_dir, "DM", "methylKit_MNcor_profile.html"))

#-----------------------------------------------------------------------------------------------------# 
# ~ DSS ----

# Make design matrix formula; correcting for Subject, determining the effect of Tissue
s_CovFormulaImportantSV = ~ Group + Age + factor(Gender) + PMI + factor(Flowcell)

# make design matrix
temp_design = model.matrix(as.formula(s_CovFormulaImportantSV),data=pheno_subset_benchmark)
colnames(temp_design) = make.names(colnames(temp_design))
head(temp_design)


# Differential 5-mC and 5-hmC locus analysis
DM_analysis_single_DSS <- function(mod_type, counts_file, coverage) {
  
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
  
  tryCatch(saveRDS(BSobj, file = file.path(paste0(s_ROOT_dir,s_out_folder,"DM/BSobj_",mod_type, ".rds"))), error = function(e) warning(e$message))
  
  # Free up memory
  rm(counts, position)
  gc()
  
  
  # Run DML
  cat(paste0("Running DMLfit...", "\n\n"))

  DMLfit = DMLfit.multiFactor(BSobj, design=pheno_subset_benchmark, formula=s_CovFormulaImportantSV, smoothing = FALSE)
  tryCatch(saveRDS(DMLfit, file = file.path(paste0(s_ROOT_dir,s_out_folder,"DM/DMLfit_",mod_type,".rds"))), error = function(e) warning(e$message))
  
  # Contrast analysis
  cat(paste0("Contrast analysis on DMLfit"),"\n\n")

  loop_contrast <- "GroupAD"
    s_coeff = loop_contrast
    
    # Obtain statistics for contrasts
    
    DM_Results <- DMLtest.multiFactor(DMLfit, coef = s_coeff)
    print(head(DM_Results[order(DM_Results$pvals, decreasing = FALSE),]))
    
    # Histogram of p values and adjusted p values
    destination <- file.path(s_OUT_dir, paste0("Plots/DM_hist_DSS_",loop_contrast, mod_type, ".pdf"))
    pdf(file=destination)
    
    hist(DM_Results$pvals, xlim=c(0,1), main =s_coeff)
    hist(DM_Results$fdrs, xlim=c(0,1), main =s_coeff)
    
    dev.off()
    
    # Descriptives
    cat(paste0("Minimum adjusted p value of ", s_coeff, ": " ,round(min(DM_Results$fdrs),digits = 3), "\n\n"))
    
    cat(paste0("Contrast: ",s_coeff,".\n  number of features below TH (0.01)\n  Pvalue:",
               sum(DM_Results$pvals<0.01, na.rm =TRUE)," (",round((sum(DM_Results$pvals<0.01, na.rm=TRUE)/dim(DM_Results)[1])*100,0),"%)",
               "\n  AdjPvalue:",sum(DM_Results$fdrs<0.05, na.rm = TRUE)," (",round((sum(DM_Results$fdrs<0.05, na.rm = TRUE)/dim(DM_Results)[1])*100,0),"%)","\n\n"))
    
    tryCatch(saveRDS(DM_Results, file = file.path(paste0(s_ROOT_dir,s_out_folder,"DM/DM_Results_DSS_", mod_type, "_",loop_contrast, ".rds"))), error = function(e) warning(e$message))
    # tryCatch(fwrite(DM_Results,file = paste0(s_ROOT_dir,s_out_folder,"DM/DM_Results_",mod_type, "_", loop_contrast,".csv")), error = function(e) warning(e$message))
    

  cat(paste0("Save results"),"\n\n")
  
  rm(DMLfit)
  gc()
}

# Run the DM analysis
# coverage <- readRDS(file.path(s_OUT_dir, paste0("Bedmethyl_all/", "bedfiles_scores_all_", cutoff, "X_pass_nosexchr_subset.rds")))
# coverage <- as.matrix(coverage)

DSS_profile <- profvis(
  DM_analysis_single_DSS(mod_type = "methyl", 
                  counts_file = file.path(s_OUT_dir, paste0("Bedmethyl_all/", "bedfiles_methyl_count_all_", cutoff, "X_pass_nosexchr_subset.rds")),
                  coverage = coverage)
)

saveRDS(DSS_profile, file.path(s_OUT_dir, "DM", "DSS_profile.rds"))
htmlwidgets::saveWidget(DSS_profile, file = file.path(s_OUT_dir, "DM", "DSS_profile.html"))


# DM_analysis_single_DSS(mod_type = "hydroxymethyl", 
#                   counts_file = file.path(s_OUT_dir, paste0("Bedmethyl_all/", "bedfiles_hydroxymethyl_count_all_", cutoff, "X_pass_nosexchr_subset.rds")),
#                   coverage=coverage)
# 
# DM_analysis_single_DSS(mod_type = "modified", 
#                   counts_file = file.path(s_OUT_dir, paste0("Bedmethyl_all/", "bedfiles_modified_count_all_", cutoff, "X_pass_nosexchr_subset.rds")),
#                   coverage=coverage)

# NA values are probably because of design matrix not being full rank and model not being able to converge
# The internal EM algorithm fails to converge for reasons it doesn't report.


#-----------------------------------------------------------------------------------------------------#
#				8. List of DM Results                                               ----
#-----------------------------------------------------------------------------------------------------#
# List of DM Results
temp = list.files(file.path(s_OUT_dir, "DM/"), pattern="DM_Results_",
                  full.names = TRUE)



file_base_name <- sub(".rds", "", sub("DM_Results_", "", basename(temp)))

DM_Results_list = pblapply(temp, function(file){
  DM_Results <- readRDS(file)
  
  if(identical(names(DM_Results), c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B" ))) {
    colnames(DM_Results) <- c("logFC", "AveExpr", "stat", "pvals", "fdrs", "B" )
  } 
  
  return(DM_Results)
})


names(DM_Results_list) <- file_base_name


#-----------------------------------------------------------------------------------------------------#
#				9. QQ plot                                               ----
#-----------------------------------------------------------------------------------------------------#

# # function to calculate lambda
# calc_lambda <- function(p_values){
#   p_values <- na.omit(p_values)
#   median(qchisq(1-p_values,df= 1))/qchisq(0.5,df = 1)
# }

# destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DM1_qqplots_", DM_Results_method,".svg", sep="")
# svglite(destination)
# qqman::qq(DM_Results$pvals)
# text(0.5,4, paste("lambda","=",  signif(alpha, digits = 3)))
#
# dev.off()

# Citation: Kamil Slowikowski, https://slowkow.com/notes/ggplot2-qqplot/
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



temp = names(DM_Results_list)


qqplots <- pblapply(temp, function(DM_Results_method){
  
  # Load DM Results
  DM_Results <- DM_Results_list[[DM_Results_method]]
  
  DM_Results_method_name <- gsub("_", " ", sub("GroupAD", "", DM_Results_method))
  
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
  
  destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DM1_qqplots_", DM_Results_method_name,".svg", sep="")
  svglite(destination)
  
  print(image)
  dev.off()
  
  return(image)
})

destination <- paste0(s_ROOT_dir, s_out_folder, "Plots/DM1_qqplots_combined2.svg", sep="")
svglite(destination)

ggarrange(plotlist = qqplots, labels = c("A", "B", "C", "D", "E"))

dev.off()


#-----------------------------------------------------------------------------------------------------#
#			10. Check known methylation difference                                                        ----
#-----------------------------------------------------------------------------------------------------#

# True differentially methylated CpGs 
true_diff_cpgs_names <- rownames(coverage)[true_diff_cpgs]

true_neg_cpg_names <- rownames(coverage)[-c(true_diff_cpgs)]

temp = names(DM_Results_list)


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

fdr_stats <- pblapply(temp, function(DM_Results_method){
  
  # Load DM Results
  DM_Results <- DM_Results_list[[DM_Results_method]]
  
  na_values <- sum(is.na(DM_Results$fdrs))
  
  # Get Significant results
  DM_Results_sign <- DM_Results[DM_Results$fdrs < 0.05, ]
  cpg_ids_sign <- get_cpg_ids(DM_Results_sign)
  
  
  # Get not significant Cpg sites
  DM_Results_notsign <- DM_Results[DM_Results$fdrs >= 0.05, ]
  cpg_ids_notsign <- get_cpg_ids(DM_Results_notsign)
  
  
  
  # True Positive Rate (TPR)
  true_pos <- length(intersect(cpg_ids_sign, true_diff_cpgs_names))
  
  # False Discovery Rate (FDR)
  false_pos <- length(setdiff(cpg_ids_sign, true_diff_cpgs_names))
  
  true_neg <- length(intersect(cpg_ids_notsign, true_neg_cpg_names))
  
  false_neg <- length(setdiff(cpg_ids_notsign, true_neg_cpg_names))
  
  TPR <- true_pos / (true_pos + false_neg)
  
  FDR <- false_pos / (true_pos+ false_pos)
  

  stats <- data.frame(true_pos, false_pos, true_neg, false_neg, TPR, FDR, na_values)
  return(stats)
})

names(fdr_stats) <- names(DM_Results_list)
fdr_stats <- do.call("rbind", fdr_stats)

fdr_stats





# Plot coverages densities
# cat(paste0("Plotting coverage densities", "\n\n"))
# destination = file.path(s_OUT_dir, "Plots", "methylBase_obj_",mod_type,"_covdensities.svg")
# svglite(destination) 
# 
# plot(density(methylBase_obj_sim@.Data[[methylBase_obj_sim@coverage.index[1]]]), main = "")
# 
# for (cov_index in methylBase_obj_sim@coverage.index[-1]) {
#   lines(density(methylBase_obj_sim@.Data[[cov_index]]))
# }
# 
# dev.off()




df_p_val <- data %>%
  rstatix::group_by(.data[[group_var]]) %>%
  rstatix::t_test(as.formula(paste0(y_var, " ~ ", x_var))) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
  rstatix::add_significance(p.col = "p.adj") %>%
  rstatix::add_xy_position(x = group_var, dodge = 0.8) # important for positioning!

# df_p_val2 <- data %>%
#   # rstatix::group_by(.data[[x_var]]) %>% 
#   rstatix::t_test(as.formula(paste0(y_var, " ~ ", group_var))) %>%
#   rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
#   rstatix::add_significance(p.col = "p.adj") %>% 
#   rstatix::add_xy_position(x = x_var)


df_p_val2 <- data %>%
  rstatix::group_by(.data[[x_var]]) %>%
  rstatix::t_test(as.formula(paste0(y_var, " ~ ", group_var))) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "BH") %>%
  rstatix::add_significance(p.col = "p.adj") %>%
  rstatix::add_xy_position(x = x_var)

# names(df_p_val2)[1] <- 
names(df_p_val2)[3] <- group_var

# max_dep <- max(data[,y_var], na.rm=TRUE)
# seq <- round(max_dep / 10, digits = 0)

df_p_val <- df_p_val[df_p_val$p.adj.signif != "ns",]
# df_p_val <- cbind(df_p_val, "y.position" = seq(max_dep, by = seq, length.out = nrow(df_p_val)))



print(df_p_val)
print(df_p_val2)

set.seed(32)


image <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], 
                          colour = .data[[group_var]], shape = .data[[group_var]])) +
  geom_boxplot(fill = "white") + 
  # geom_jitter(width=0.2, size = 2)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.8) +
  theme_gray()+
  # theme(axis.line = element_line(color='black'),
  #       plot.background = element_blank(),
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank()) +
  theme(legend.position = "right") +
  
  xlab(xlab_text) + 
  ylab(ylab_text) +  
  # add_pvalue(df_p_val, y.position = "y.position", colour = "black",
  #            tip.length= 0, bracket.nudge.y = 2,
  #            xmin = "xmin",
  #            xmax = "xmax") +
  # 
  add_pvalue(df_p_val2, y.position = "y.position", colour = "black",
             xmin = "xmin", xmax = "xmax", tip.length = 0, bracket.nudge.y = 2)+
  
  scale_color_manual(name = group_var, values = my_colors, labels = label_text) +
  scale_shape_manual(name = group_var, values = my_shapes, labels = label_text) 
