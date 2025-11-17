#!bin/bash

# Install Dorado stand alone software, Mosdepth, and Modkit, Samtools

#################################################################################
#                   Create variables
#################################################################################

REFERENCE_DIR='/mnt/e/Reference-files'
WORKING_DIR='/mnt/e/add/path/to/data'
FLOWCELL='name_of_flowcell'


#################################################################################
#                  Dorado post run analysis
#################################################################################
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



