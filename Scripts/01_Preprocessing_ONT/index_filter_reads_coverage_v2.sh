#!bin/bash

# sudo mount -t drvfs E: /mnt/e
# Use dos2unix command to convert file to bash file

# Install modkit, samtools

#--------------------------------------------------------------------------------
#                   Create variables
#--------------------------------------------------------------------------------
WORKING_DIR=$1
REFERENCE_DIR=$2
FLOWCELL=$3
ALIGNED_DIR="${WORKING_DIR}/bam/bam_aligned"


#--------------------------------------------------------------------------------
#                   Aligned BAM files
#--------------------------------------------------------------------------------

# Move to working directory 
cd "${WORKING_DIR}"

# Make bedmethyl directory 
mkdir filtered_reads

# Move to directory with aligned BAM files
cd "${ALIGNED_DIR}"

BAMFILES=($(ls *.bam))

for b in ${BAMFILES[@]}
do 
    # Index bam file
    echo "Running samtools index for $b"
    chmod +x "${ALIGNED_DIR}/$b"
    chmod +x "${ALIGNED_DIR}/$b.bai"

    samtools index "${ALIGNED_DIR}/$b"

    # Make directory for specific bam file
    BASENAME=($(basename $b .bam))
    if [[ "${BASENAME}" == *_SQK-NBD114-24_barcode* ]]; then
        BASENAME=${BASENAME#*_} # Removes everything before the first "_"
        echo "$BASENAME"
    else
        echo "No prefix in: ${BASENAME}"
    fi

    # BASENAME=${BASENAME#$PREFIX}
    mkdir "${ALIGNED_DIR}/${BASENAME}"
    cd "${ALIGNED_DIR}/${BASENAME}"

    # Calculate number of reads
   echo "Calculating global number of reads (mapped and unmapped) for $b"
   samtools view -c "${ALIGNED_DIR}/$b" --output count_global.txt --with-header
   echo "Calculating target region mapped number of reads for $b"
   samtools view -c -F 4 "${ALIGNED_DIR}/$b" --target-file "${REFERENCE_DIR}/RRMS_human_hg38.bed" --output count_inside_mapped.txt --with-header
   
    # Run mosdepth
   echo "Running mosdepth for $b"
   mosdepth -x -t 8 -n -b "${REFERENCE_DIR}/RRMS_human_hg38.bed" ${BASENAME} "${ALIGNED_DIR}/$b"

   # Filter reads
    echo "Filter reads of ${BASENAME}"
    samtools view -h -q 55 "${ALIGNED_DIR}/$b" | awk 'length($10) >= 1000 || $1 ~ /^@/' | samtools view -bS > "${WORKING_DIR}/filtered_reads/filtered_${BASENAME}.bam"

    echo "Running samtools index for filtered_${BASENAME}.bam"
    samtools index "${WORKING_DIR}/filtered_reads/filtered_${BASENAME}.bam"

done


#--------------------------------------------------------------------------------
#                 Filtered BAM files
#--------------------------------------------------------------------------------

FILTERED_DIR="${WORKING_DIR}/filtered_reads"

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
