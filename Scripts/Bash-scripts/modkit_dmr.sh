#!bin/bash

# Install Dorado stand alone software, Mosdepth, and Modkit
# Install modkit: cargo install --git https://github.com/nanoporetech/modkit.git

# Add Modkit to PATH
#export PATH=$PATH:/mnt/c/modkit/bin

export PATH=$PATH:/mnt/c/Users/lucp13068/modkit/bin

# Create variables
REFERENCE_DIR='/mnt/g/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/Reference-files'
WORKING_DIR='/mnt/g/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/Data/RESEARCH_GSM0172RRMS_ALL'
BATCH='RESEARCH_GSM0172RRMS_ALL'

# Move to working directory 
cd "${WORKING_DIR}"

# Make bedmethyl directory 
mkdir dmr 

bgzip -o "${WORKING_DIR}/dmr/RESEARCH_GSM0172RRMS_04112024_filtered_SQK-NBD114-24_barcode11.bam.bed.gz" "${WORKING_DIR}/RESEARCH_GSM0172RRMS_04112024_filtered_SQK-NBD114-24_barcode11.bam.bed"
tabix -p bed "${WORKING_DIR}/dmr/RESEARCH_GSM0172RRMS_04112024_filtered_SQK-NBD114-24_barcode11.bam.bed.gz.tbi"

bgzip -o "${WORKING_DIR}/dmr/RESEARCH_GSM0172RRMS_04112024_filtered_SQK-NBD114-24_barcode12.bam.bed.gz" "${WORKING_DIR}/RESEARCH_GSM0172RRMS_04112024_filtered_SQK-NBD114-24_barcode12.bam.bed"
tabix -p bed "${WORKING_DIR}/dmr/RESEARCH_GSM0172RRMS_04112024_filtered_SQK-NBD114-24_barcode12.bam.bed.gz.tbi"

modkit dmr pair \
  -a "${WORKING_DIR}/dmr/RESEARCH_GSM0172RRMS_04112024_filtered_SQK-NBD114-24_barcode11.bam.bed.gz" \
  --index-a "${WORKING_DIR}/dmr/RESEARCH_GSM0172RRMS_04112024_filtered_SQK-NBD114-24_barcode11.bam.bed.gz.tbi" \
  -b "${WORKING_DIR}/dmr/RESEARCH_GSM0172RRMS_04112024_filtered_SQK-NBD114-24_barcode12.bam.bed.gz" \
  --index-b "${WORKING_DIR}/dmr/RESEARCH_GSM0172RRMS_04112024_filtered_SQK-NBD114-24_barcode12.bam.bed.gz.tbi" \
  -o "dmr_output_AD_CTL.bed" \
  -r "${REFERENCE_DIR}/RRMS_human_hg38.bed" \
  --ref "${REFERENCE_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" \
  --base C \
  --threads 10 \
  --log-filepath "${WORKING_DIR}/dmr/dmr.log"