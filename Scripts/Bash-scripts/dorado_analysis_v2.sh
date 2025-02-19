#!bin/bash

# Install Dorado stand alone software, Mosdepth, and Modkit

# Create variables
REFERENCE_DIR='/mnt/g/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/Reference-files'
PROJECT_NAME='NBD_RRMS'
WORKING_DIR='/mnt/g/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/Data/RESEARCH_GSM0172RRMS_16092024'

# Move to working directory 
cd ${WORKING_DIR}

# Make BAM file directory 
mkdir bam 

# Basecalling 5mC/5hmC and trimming read ends
dorado basecaller hac,5mCG_5hmCG pod5/ --trim all --kit-name SQK-NBD114-24 --barcode-both-ends --emit-summary > bam/calls.bam

# Demultiplexing
dorado demux --output-dir bam/bam_demuxed --no-classify bam/calls.bam --emit-summary

# Alignment to reference
dorado aligner ${REFERENCE_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna bam/bam_demuxed --output-dir bam/bam_aligned --emit-summary

# Check coverage 
#mosdepth -x -t 8 -n -b ${REFERENCE_DIR}/RRMS_human_hg38.bed ${PROJECT_NAME} bam/bam_aligned


