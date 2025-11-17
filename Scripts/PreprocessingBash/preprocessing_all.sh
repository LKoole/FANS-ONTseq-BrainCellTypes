#!bin/bash

# sudo mount -t drvfs G: /mnt/g

#--------------------------------------------------------------------------------
#               MODKIT
#--------------------------------------------------------------------------------
# Add Modkit to PATH
#export PATH=$PATH:/mnt/c/modkit/bin
export PATH=$PATH:/add/path/to/user


#--------------------------------------------------------------------------------
#               DEFINE DIRECTORIES
#--------------------------------------------------------------------------------
# Define parent directory (root)
PARENT_DIR='/definy/path/to/root/directory'

# Define data directory containing all sequencing run folders
DATA_DIR="${PARENT_DIR}/Data"

# Define scripts directory
SCRIPTS_DIR="${PARENT_DIR}/Scripts"

# Define reference directory containing the reference genome
REFERENCE_DIR="${PARENT_DIR}/Reference-files"


#--------------------------------------------------------------------------------
#               LIST OF RUNS TO PREPROCESS
#--------------------------------------------------------------------------------

# Check if the run_list file exists 
RUN_LIST="${DATA_DIR}/run_list.txt"
if [[ ! -f "$RUN_LIST" ]]; then
    echo "Error: File ${RUN_LIST} not found!"
    exit 1
fi

dos2unix "${RUN_LIST}" # If first time using the txt file

# Remove empty lines
sed -i '/^$/d' "${RUN_LIST}"

#--------------------------------------------------------------------------------
#               PREPROCESSING
#--------------------------------------------------------------------------------
# Loop through each sequencing run listed in the file
while IFS= read -r FLOWCELL; do

    # Retrieve working directory based on sequencing runs listed
    WORKING_DIR="${DATA_DIR}/${FLOWCELL}"
    chmod +x "${WORKING_DIR}"

    if [[ -z "${FLOWCELL}" ]]; then
        echo "Skipping empty line in run_list.txt"
        continue
    fi

    if [[ ! -d "${WORKING_DIR}" ]]; then  # Ensure it's a directory
        
        echo "Warning: Directory ${WORKING_DIR} does not exist. Skipping..."
        continue
    fi

    echo "Processing sequencing run in: ${WORKING_DIR}"
        
    # Run the filtering and coverage script
    bash "${SCRIPTS_DIR}/Bash-scripts/index_filter_reads_coverage_v2.sh" "$WORKING_DIR" "$REFERENCE_DIR" "$FLOWCELL"
    
    # Run the Modkit script
    bash "${SCRIPTS_DIR}/Bash-scripts/bedmethyl_v7.sh" "$WORKING_DIR" "$REFERENCE_DIR" "$FLOWCELL"

    # Summarise the number of reads and coverage
    bash "${SCRIPTS_DIR}/Bash-scripts/summary_reads_coverage.sh" "$WORKING_DIR" "$REFERENCE_DIR" "$FLOWCELL"


done < "${RUN_LIST}"