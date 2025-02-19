#!bin/bash

# Install bcftools

# wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
# bunzip2 bcftools-1.19.tar.bz2
# tar -xvf bcftools-1.19.tar
# rm bcftools-1.19.tar
# cd bcftools-1.19/
# ./configure --prefix=/usr/local/bin/bcftools-1.
# make
# sudo make install
# export PATH=$PATH:/usr/local/bin/bcftools-1.19/bin

# dos2unix

# Create variables
REFERENCE_DIR='/mnt/g/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/Reference-files'
WORKING_DIR='/mnt/g/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/Data/RESEARCH_GSM0172RRMS_31072024'
FLOWCELL='RESEARCH_GSM0172RRMS_31072024'

# Move to working directory 
cd "${WORKING_DIR}"

# Make SNP calling directory 
mkdir SNPcalling

# Move to directory with filtered BAM files
cd "${WORKING_DIR}/filtered_reads"

BAMFILES=($(ls *.bam))

# Calling SNPs with bcftools
for bam in ${BAMFILES[@]}
do 
    chmod +x $bam
    chmod +x $bam.bai
    SAMPLE_NAME=$(basename $bam)
    echo "Running bcftools SNP calling for ${SAMPLE_NAME}..."

    # creating a vcf file and calling SNPs
    bcftools mpileup -f "${REFERENCE_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" --config ont $bam | bcftools call --ploidy 2 -mv -Ob -o "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_var.vcf" -O v 
    
    # normalize the file and remove duplicates
    bcftools norm -d all -f "${REFERENCE_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" -o "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_norm.bcf" -O u "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_var.vcf"

    # remove low quality and low depth (<30) variants or high depth (>100, might represent variation in copy number repeats) variants 
    bcftools filter -e "QUAL<55 || DP<30 || DP>100" -Ob -o "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_norm_filter.bcf" "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_norm.bcf"

    # get a IGV compatible copy (uncompressed vcf) of bcf compressed file
    bcftools view -Ov -o "${WORKING_DIR}/SNPcalling/${FLOWCELL}_${SAMPLE_NAME}_igv.vcf" "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_norm_filter.bcf"

    # index the file (csi file)
    bcftools index "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_norm_filter.bcf"

    # stats 
    bcftools stats -f "${REFERENCE_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_norm_filter.bcf" > "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_norm_filter.stats"

    mkdir "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_plots"

    plot-vcfstats -T "Variants in ONT" -P -p "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_plots" "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_norm_filter.stats"
    # if above command says matplotlib not installed, then run following command
    # pip3 install matplotlib
    # if pip or pip3 not installed, install that first
done










#bcftools mpileup -f "${REFERENCE_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" --config ont -o "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_raw.vcf" -O v --threads 8 $bam
# getting only variants
#bcftools call --ploidy 2 -v -m -o "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_var.vcf" -O v "${WORKING_DIR}/SNPcalling/${SAMPLE_NAME}_raw.vcf" 


# Check https://angus.readthedocs.io/en/2013/snp_tutorial.html
# Check https://medium.com/@manabeel.vet/a-beginners-guide-to-genomic-data-analysis-variant-calling-ad8515eebddf 




# # Install DeepVariant

# # Follow instructions on: https://cloud.google.com/sdk/docs/install#deb
# git clone https://git.launchpad.net/ubuntu/+source/python3-stdlib-extensions

# git clone https://github.com/google/deepvariant.git
# cd deepvariant

# ./build-prereq.sh

# ./build_and_test.sh
# sudo apt install -y unzip
# curl -Lo bazelisk https://github.com/bazelbuild/bazelisk/releases/download/v1.12.0/bazelisk-linux-amd64
# chmod +x bazelisk
# sudo mv bazelisk /usr/local/bin/bazel




# Create a vcf file (pileup) and get the variants (call)
    #bcftools mpileup -f "${REFERENCE_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" --config ont $bam | bcftools call --ploidy 2 -mv -Ob -o "${WORKING_DIR}/SNPcalling/snp_calling_${SAMPLE_NAME}.vcf" -O v 
    # Filter for SNPs with read depth higher than 100 (might represent variation in copy number repeats)
    # And convert to VCF file
    # bcftools view ${WORKING_DIR}/${PROJECT_NAME}/SNPcalling/snp_calling_${SAMPLE_NAME}.bcf | vcfutils.pl varFilter -D100 > ${WORKING_DIR}/${PROJECT_NAME}/SNPcalling/snp_calling_${SAMPLE_NAME}.vcf
# bcftools view -vcg - > ${WORKING_DIR}/${PROJECT_NAME}/SNPcalling/snp_calling_${SAMPLE_NAME}.bcf

# -b for BCF format
# -c SNP calling
# -v Only output potential variant sites 
# -g Call genotypes for each sample