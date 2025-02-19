#!bin/bash

# sudo mount -t drvfs E: /mnt/e
# Use dos2unix command to convert file to bash file

# Install modkit, samtools

#################################################################################
#                   Create variables
#################################################################################

REFERENCE_DIR='/mnt/g/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/Reference-files'
WORKING_DIR='/mnt/g/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/Data/RESEARCH_GSM0172RRMS_27062024'
FLOWCELL='RESEARCH_GSM0172RRMS_27062024'
ALIGNED_DIR=${WORKING_DIR}/bam/bam_aligned

#################################################################################
#                   Aligned BAM files
#################################################################################

# Move to working directory 
cd "${WORKING_DIR}"

# Make bedmethyl directory 
mkdir filtered_reads2

# Move to directory with aligned BAM files
cd "${ALIGNED_DIR}"

BAMFILES=($(ls *.bam))

for b in ${BAMFILES[@]}
do 
    # Index bam file
    echo "Running samtools index for $b"
    chmod +x "${ALIGNED_DIR}/$b"
    chmod +x "${ALIGNED_DIR}/$b.bai"

    # samtools index "${ALIGNED_DIR}/$b"

    # Make directory for specific bam file
    # BASENAME=($(basename $b .bam))
    # mkdir "${ALIGNED_DIR}/${BASENAME}"
    # cd "${ALIGNED_DIR}/${BASENAME}"

    # Calculate number of reads
    #    echo "Calculating global number of reads (mapped and unmapped) for $b"
    #    samtools view -c "${ALIGNED_DIR}/$b" --output count_global.txt --with-header
    #    echo "Calculating target region mapped number of reads for $b"
    #    samtools view -c -F 4 "${ALIGNED_DIR}/$b" --target-file "${REFERENCE_DIR}/RRMS_human_hg38.bed" --output count_inside_mapped.txt --with-header
   
    # Run mosdepth
    #    echo "Running mosdepth for $b"
    #    mosdepth -x -t 8 -n -b "${REFERENCE_DIR}/RRMS_human_hg38.bed" $b "${ALIGNED_DIR}/$b"

   # Filter reads
    echo "Filter reads of ${BASENAME}"
    samtools view -h -q 55 "${ALIGNED_DIR}/$b" | awk 'length($10) >= 400 || $1 ~ /^@/' | samtools view -bS > "${WORKING_DIR}/filtered_reads2/filtered_400bp_$b"

    samtools view -h -q 55 "${ALIGNED_DIR}/$b" | awk 'length($10) >= 2 || $1 ~ /^@/' | samtools view -bS > "${WORKING_DIR}/filtered_reads2/filtered_2bp_$b"
    echo "Running samtools index for filtered_$b"
    samtools index "${WORKING_DIR}/filtered_reads2/filtered_400bp_$b"
    samtools index "${WORKING_DIR}/filtered_reads2/filtered_2bp_$b"

done

#################################################################################
#                   Filter reads 
#################################################################################

# # Move to working directory 
# cd ${WORKING_DIR}

# # Make bedmethyl directory 
# mkdir filtered_reads

# # Move to aligned BAMs directory 
# cd ${ALIGNED_DIR}

# BAMFILES=($(ls *.bam))

# # Filter reads for mapping quality of 55 or higher and sequence length of at least 1kb 

# for f in ${BAMFILES[@]}
# do 
#     chmod +x $f 
#     chmod +x $f.bai
#     SAMPLE_NAME=$(basename $f .bam)
#     echo "Running samtools view for ${SAMPLE_NAME}..."
#     samtools view -h -q 55 $f | awk 'length($10) >= 1000 || $1 ~ /^@/' | samtools view -bS > ${WORKING_DIR}/filtered_reads/filtered_$f
#     echo "Running samtools index for filtered_$f"
#     samtools index ${WORKING_DIR}/filtered_reads/filtered_$f
# done

#################################################################################
#                 Filtered BAM files
#################################################################################

FILTERED_DIR="${WORKING_DIR}/filtered_reads2"

# Move to working directory 
cd "${FILTERED_DIR}"

BAMFILES=($(ls *.bam))

# Create strand-aggregated methylation frequencies for all CpGs 
for f in ${BAMFILES[@]}
do 
    BASENAME=($(basename $f .bam))
   mkdir "${FILTERED_DIR}/${BASENAME}"
   cd "${FILTERED_DIR}/${BASENAME}"
   chmod +x "${FILTERED_DIR}/$f"

   echo "Running mosdepth for $f"
   mosdepth -x -t 8 -n -b "${REFERENCE_DIR}/RRMS_human_hg38.bed" $f "${FILTERED_DIR}/$f"

   echo "Calculating global number of reads (mapped and unmapped) for $f"
   samtools view -c "${FILTERED_DIR}/$f" --output count_global.txt --with-header

   echo "Calculating target region mapped number of reads for $f"
   samtools view -c -F 4 "${FILTERED_DIR}/$f" --target-file "${REFERENCE_DIR}/RRMS_human_hg38.bed" --output count_inside_mapped.txt --with-header
done
