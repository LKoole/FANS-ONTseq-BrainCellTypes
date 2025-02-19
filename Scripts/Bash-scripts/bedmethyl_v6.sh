#!bin/bash

# Install Dorado stand alone software, Mosdepth, and Modkit
# Install modkit: cargo install --git https://github.com/nanoporetech/modkit.git

# Add Modkit to PATH
#export PATH=$PATH:/mnt/c/modkit/bin

export PATH=$PATH:/mnt/c/Users/lucp13068/modkit/bin

# Create variables
REFERENCE_DIR='/mnt/g/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/Reference-files'
WORKING_DIR='/mnt/g/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/Data/RESEARCH_GSM0172RRMS_31072024'
BATCH='RESEARCH_GSM0172RRMS_31072024'

# Move to working directory 
cd "${WORKING_DIR}"

# Make bedmethyl directory 
mkdir bedmethyl

cd "${WORKING_DIR}/filtered_reads"

BAMFILES=($(ls *.bam))

# Create strand-aggregated methylation frequencies for all CpGs 
for f in ${BAMFILES[@]}
do 
    chmod +x $f 
    chmod +x $f.bai
    SAMPLE_NAME=$(basename $f)
    echo "Running modkit pileup for ${SAMPLE_NAME}..."
    modkit pileup "${WORKING_DIR}/filtered_reads/${SAMPLE_NAME}.bam" "${WORKING_DIR}/bedmethyl/${BATCH}_${SAMPLE_NAME}.bam.bed" --cpg --combine-strands --threads 10 --ref "${REFERENCE_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" --include-bed "${REFERENCE_DIR}/RRMS_human_hg38.bed" --no-filtering
    # modkit summary ${WORKING_DIR}/bam/bam_aligned/${SAMPLE_NAME}.bam --include-bed ${REFERENCE_DIR}/RRMS_human_hg38.bed --filter-threshold C:0.9 --mod-thresholds m:0.9 --mod-thresholds h:0.9 --tsv > modkit_summary_${SAMPLE_NAME}_threshold.txt
done
