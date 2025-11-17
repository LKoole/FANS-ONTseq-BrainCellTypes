#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#		Settings.R
#
# Author comment:
#	Lisa Koole

#-----------------------------------------------------------------------------------------------------#
#							Comment
#-----------------------------------------------------------------------------------------------------#
# Use '<<-' to assign to global envir. Do this here or settings might be lost


#-----------------------------------------------------------------------------------------------------#
#							SEED
#-----------------------------------------------------------------------------------------------------#
# set seed for consistency
s_seed <<- 42


#-----------------------------------------------------------------------------------------------------#
#							ROOT FOLDER
#-----------------------------------------------------------------------------------------------------#
# define ROOT here 
# >>> change this to your path! <<<
s_ROOT_dir <<- "G:/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/"

#-----------------------------------------------------------------------------------------------------#
#							Folder settings
#-----------------------------------------------------------------------------------------------------#
s_ROOT_current_folder_name <<- "RESEARCH_GSM0172RRMS_ALL" # Project

s_figure_folder_name <<- "Plots"

s_bedfiles_folder <- paste0(s_ROOT_dir,"Data/RESEARCH_GSM0172RRMS_ALL") # data

s_vcfs_folder <- paste0(s_ROOT_dir,"Data/RESEARCH_GSM0172RRMS_ALL/VCF_files/") # vcf files

s_analysis_run <- "31102025_analysis_Benchmarking"


#-----------------------------------------------------------------------------------------------------#
#							DEFINE OUTPUT FOLDER
#-----------------------------------------------------------------------------------------------------#

# Folder to store project
s_project_folder <- paste0(s_ROOT_dir, "Results/",s_ROOT_current_folder_name,"/")

# Change this to generate a different run, stored in new location
s_out_folder <<- paste0("Results/",s_ROOT_current_folder_name,"/", s_analysis_run, "/")

# define out location
s_OUT_dir <<- paste0(s_ROOT_dir,s_out_folder)

# define figure folder
s_figure_folder <<- paste0(s_ROOT_dir,s_out_folder,s_figure_folder_name,"/")


#-----------------------------------------------------------------------------------------------------#
#							MAKE OUTPUT FOLDERS IF NEEDED
#-----------------------------------------------------------------------------------------------------#
# Create project folder
if(!dir.exists(paste0(s_project_folder))){dir.create(file.path(paste0(s_project_folder)))}

# Create main output folder (per analysis)
if(!dir.exists(paste0(s_OUT_dir))){dir.create(file.path(paste0(s_OUT_dir)))}

# Bedmethyl folder
if(!dir.exists(paste0(s_project_folder, "Bedmethyl"))){dir.create(file.path(paste0(s_project_folder, "Bedmethyl")))}

# Create subfolders
temp_all_dirs = c(s_figure_folder_name, "Bedmethyl_all", "VCF_all", "QC","DM","Pheno")
lapply(temp_all_dirs,function(i){if(!dir.exists(paste0(s_OUT_dir,i))){dir.create(file.path(paste0(s_OUT_dir,i)))}})


#-----------------------------------------------------------------------------------------------------#
#							GET FUNCTIONS
#-----------------------------------------------------------------------------------------------------#
# source functions
# source(paste0(s_ROOT_dir,"Scripts\\.MAIN\\Functions.R"))


#-----------------------------------------------------------------------------------------------------#
#							Make SETTINGS file (IF NOT EXISTS)
#-----------------------------------------------------------------------------------------------------#
if(!file.exists(paste0(s_OUT_dir,"SETTINGS.txt"))){
  temp_settings_names = ls()[grep(ls(),pattern = "^s_")]
  
  temp_set = data.frame(variable = temp_settings_names, value = NA)
  
  for(i in 1:nrow(temp_set)){
    temp_set[i,"value"] =  get(temp_set[i,"variable"] )
  }
  
  write.table("SETTINGS:",file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
  write.table(temp_set,file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
  write.table("\nVERSION:",file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table(data.frame(name=names(version),value=unlist(version)),file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table("\nSYSTEM INFO:",file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table(data.frame(name=names(Sys.info()),value=unlist(Sys.info())),file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table("\nBASE PACKAGES:",file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table(sessionInfo()$basePkgs,file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table("\nOTHER PACKAGES:",file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table(names(sessionInfo()$otherPkgs),file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
  
}

