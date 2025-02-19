#!bin/bash

# Install NanoPlot
# pip install NanoPlot

# Create variables
REFERENCE_DIR='/mnt/e/reference-files'
PROJECT_NAME='NBD_RRMS'
WORKING_DIR='/mnt/e/data/RESEARCH_GSM0172RRMS_16092024/'

# Move to working directory 
cd ${WORKING_DIR}

NanoPlot --summary $WORKING_DIR/sequencing_summary_PAW50995_392b1e61_a1e48a88.txt --loglength -o $WORKING_DIR/summary-plots-log-transformed --barcoded


# pip install NanoComp
NanoComp --summary sequencing_summary_PAW50995_392b1e61_a1e48a88.txt sequencing_summary_PAW51133_ffb0a6cd_fd75adc5.txt --outdir compare-runs --barcoded

