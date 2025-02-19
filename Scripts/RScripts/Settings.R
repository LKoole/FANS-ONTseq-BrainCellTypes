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
s_ROOT_current_folder_name <<- "RESEARCH_GSM0172RRMS_ALL" # Add date
s_figure_folder_name <<- "Plots"

s_bedfiles_folder <- paste0(s_ROOT_dir,"Data/RESEARCH_GSM0172RRMS_ALL")

s_vcfs_folder <- paste0(s_ROOT_dir,"Data/RESEARCH_GSM0172RRMS_ALL/VCF_files/")


#-----------------------------------------------------------------------------------------------------#
#							DEFINE OUTPUT FOLDER
#-----------------------------------------------------------------------------------------------------#
# change this to generate a different run, stored in new location
s_out_folder <<- paste0("Results/",s_ROOT_current_folder_name,"/")

# define out loc
s_OUT_dir <<- paste0(s_ROOT_dir,s_out_folder)

# define figure folder
s_figure_folder <<- paste0(s_ROOT_dir,s_out_folder,s_figure_folder_name,"/")


#-----------------------------------------------------------------------------------------------------#
#							MAKE OUTPUT FOLDERS IF NEEDED
#-----------------------------------------------------------------------------------------------------#
temp_all_dirs = c(s_figure_folder_name,"Qualimap", "Bedmethyl", "VCF", "QC","DM","Pheno", "Run Reports")
if(!dir.exists(paste0(s_OUT_dir))){dir.create(file.path(paste0(s_OUT_dir)))}
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
