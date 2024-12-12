#!/bin/bash

################################################################################
# Script Name: HiC_abc.sh
# Description: This script processes Hi-C data using juicebox_dump.py. 
#              It iterates over a list of Hi-C files and processes them 
#              using a specified Hi-C processing pipeline.
# Creation Date: August 2024
# Version: 1.0
# Usage: bash HiC_abc.sh
################################################################################

# Load necessary modules
module load java
module load tabix
conda activate basepengfei

# Define paths and parameters
OUT_DIR="/path/to/output_directory"
HIC_PATH="/path/to/hic_files"
CHROM_SIZE="/path/to/hg38_sizes"
JUICER_PATH="/path/to/juicer_tools.jar"
ABC_PATH="/path/to/ABC-Enhancer-Gene-Prediction/src"

# Input Hi-C files and associated cell types
HIC_FILES=("example1.hic" "example2.hic") # List of Hi-C files
CELL_TYPES=("CellType1" "CellType2")      # Corresponding cell types

# Check if arrays match in length
if [ "${#HIC_FILES[@]}" -ne "${#CELL_TYPES[@]}" ]; then
  echo "Error: HIC_FILES and CELL_TYPES arrays must have the same length."
  exit 1
fi

# Loop through Hi-C files and process them
LEN=${#HIC_FILES[@]}
for ((i = 0 ; i < ${LEN} ; i++)); do
  HIC_FILE=${HIC_FILES[$i]}
  CELL_TYPE=${CELL_TYPES[$i]}

  echo "Processing ${HIC_FILE} for cell type ${CELL_TYPE}..."

  python "${ABC_PATH}/juicebox_dump.py" \
    --hic_file "${HIC_PATH}/${HIC_FILE}" \
    --juicebox "java -jar ${JUICER_PATH}" \
    --outdir "${OUT_DIR}/hic/${CELL_TYPE}/" \
    --chromosomes all
  
  echo "Completed processing for ${HIC_FILE}."
done

echo "All Hi-C files processed successfully."

