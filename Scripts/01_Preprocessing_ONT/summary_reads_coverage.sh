#!/bin/bash


#--------------------------------------------------------------------------------
#               DEFINE DIRECTORIES
#--------------------------------------------------------------------------------
# Create variables
WORKING_DIR="$1"
REFERENCE_DIR="$2"
FLOWCELL="$3"

ALIGNED_DIR="${WORKING_DIR}/bam/bam_aligned"
FILTERED_DIR="${WORKING_DIR}/filtered_reads"

# Move to working directory 
cd "${WORKING_DIR}" || exit 1


#--------------------------------------------------------------------------------
#               SUMMARY FILES
#--------------------------------------------------------------------------------
# Create summary files
output_file_aligned="summary_${3}_aligned_untarget.txt"
output_file_filtered="summary_${3}_filtered_untarget.txt"


# Write header (tab delimited)
echo -e "Flowcell\tFiltered\tBarcode\tCount_global\tCount_inside_mapped\tCoverage\tCoverage_untargeted" > "$output_file_aligned"
echo -e "Flowcell\tFiltered\tBarcode\tCount_global\tCount_inside_mapped\tCoverage\tCoverage_untargeted" > "$output_file_filtered"


#--------------------------------------------------------------------------------
#               BAM_ALIGNED FILES
#--------------------------------------------------------------------------------
# Check if bam_aligned path exists
if [ ! -d "$ALIGNED_DIR" ]; then
    echo "Warning: Directory ${ALIGNED_DIR} does not exist. Skipping..."
    continue
fi

# Extract information for each barcode directory
for barcode_dir in "$ALIGNED_DIR"/SQK-NBD114-24_barcode*; do
    
    # Check directory
    [ -d "$barcode_dir" ] || continue
    barcode=$(basename "$barcode_dir")

    if [[ "${barcode}" == *SQK-NBD114-24_barcode* ]]; then
     BASENAME=${barcode#*SQK-NBD114-24_} # Removes everything before the first "_"
        echo "$BASENAME"
    else
        echo "No prefix in: ${BASENAME}"
    fi

    echo "Extracting values for: $barcode"

    count_global_file="${barcode_dir}/count_global.txt"
    count_inside_file="${barcode_dir}/count_inside_mapped.txt"
    mosdepth_summary_file="${barcode_dir}/${barcode}.mosdepth.summary.txt" # check name 

    # sed -n '51p' "$barcode_aligned_dir/$barcode.mosdepth.summary.txt" # Printing line of coverages

    if [[ -f "$count_global_file" ]] && [[ -f "$count_inside_file" ]] && [[ -f "$mosdepth_summary_file" ]]; then

      read -r global_reads < "$count_global_file"
      read -r inside_mapped_reads < "$count_inside_file"
      coverage=$(awk -F$'\t' 'NR==51 {print $4}' "$mosdepth_summary_file")
      
      total_length=$(awk -F$'\t' 'NR==50 {print $2}' "$mosdepth_summary_file")
      total_bases=$(awk -F$'\t' 'NR==50 {print $3}' "$mosdepth_summary_file")
      total_region_length=$(awk -F$'\t' 'NR==51 {print $2}' "$mosdepth_summary_file")
      total_region_bases=$(awk -F$'\t' 'NR==51 {print $3}' "$mosdepth_summary_file")

      coverage_untargeted=$(echo "scale=3; ($total_bases-$total_region_bases)/($total_length-$total_region_length)" | bc)



      echo -e "Line to write: ${FLOWCELL}\tNo\t${BASENAME}\t${global_reads}\t${inside_mapped_reads}\t${coverage}\t${coverage_untargeted}"

      # Write to summary file
      echo -e "${FLOWCELL}\tNo\t${BASENAME}\t${global_reads}\t${inside_mapped_reads}\t${coverage}\t${coverage_untargeted}" >> "$output_file_aligned"
    fi
done

echo "Summary written to $output_file_aligned"

#--------------------------------------------------------------------------------
#              FILTERED_READS FILES
#--------------------------------------------------------------------------------
# Check if bfiltered reads path exists
if [ ! -d "$FILTERED_DIR" ]; then
    echo "Warning: Directory ${FILTERED_DIR} does not exist. Skipping..."
    continue
fi

# Extract information for each barcode directory
for barcode_dir in "$FILTERED_DIR"/filtered_SQK-NBD114-24_barcode*; do
    [ -d "$barcode_dir" ] || continue
    barcode=$(basename "$barcode_dir")

    if [[ "${barcode}" == *SQK-NBD114-24_barcode* ]]; then
     BASENAME=${barcode#*SQK-NBD114-24_} # Removes everything before the first "_"
        echo "$BASENAME"
    else
        echo "No prefix in: ${BASENAME}"
    fi

    echo "Extracting values for: $barcode"

    count_global_file="${barcode_dir}/count_global.txt"
    count_inside_file="${barcode_dir}/count_inside_mapped.txt"
    mosdepth_summary_file="${barcode_dir}/${barcode}.bam.mosdepth.summary.txt" # check name

    # sed -n '51p' "$barcode_aligned_dir/$barcode.mosdepth.summary.txt" # Printing line of coverages


    if [[ -f "$count_global_file" ]] && [[ -f "$count_inside_file" ]] && [[ -f "$mosdepth_summary_file" ]]; then

      read -r global_reads < "$count_global_file"
      read -r inside_mapped_reads < "$count_inside_file"
      coverage=$(awk -F$'\t' 'NR==51 {print $4}' "$mosdepth_summary_file")


      total_length=$(awk -F$'\t' 'NR==50 {print $2}' "$mosdepth_summary_file")
      total_bases=$(awk -F$'\t' 'NR==50 {print $3}' "$mosdepth_summary_file")
      total_region_length=$(awk -F$'\t' 'NR==51 {print $2}' "$mosdepth_summary_file")
      total_region_bases=$(awk -F$'\t' 'NR==51 {print $3}' "$mosdepth_summary_file")

      coverage_untargeted=$(echo "scale=3; ($total_bases-$total_region_bases)/($total_length-$total_region_length)" | bc)



      echo -e "Line to write: ${FLOWCELL}\tNo\t${BASENAME}\t${global_reads}\t${inside_mapped_reads}\t${coverage}\t${coverage_untargeted}"

      # Write to summary file
      echo -e "${FLOWCELL}\tYes\t${BASENAME}\t${global_reads}\t${inside_mapped_reads}\t${coverage}\t${coverage_untargeted}" >> "$output_file_filtered"
    fi
done

echo "Summary written to $output_file_filtered"
