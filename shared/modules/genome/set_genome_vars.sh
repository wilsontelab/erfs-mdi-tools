# action:
#     set environment variables for genome file access

# file prefixes for genome+sample-specific output files
export DATA_GENOME_PREFIX=${DATA_FILE_PREFIX}.${GENOME}

# genome bins file (not sample-specific)
export GENOME_BINS_BED=${TASK_DIR}/${GENOME}.bin_${BIN_SIZE}.bed.gz

# get genome size from canonical chromosomes only
# get the list of all placed chromosome sequences, including chrX, chrY, chrM, and chrEBV if present
export GENOME_SIZE=`awk '$1!~/_/ && $1!="chrEBV"{s+=$2}END{print s}' ${GENOME_FASTA}.fai`
export GENOME_CHROMS=`cat ${GENOME_FASTA}.fai | cut -f1 | grep -v _`
