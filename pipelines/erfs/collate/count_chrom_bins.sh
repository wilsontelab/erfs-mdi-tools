# action:
#     count reads in genome bins for a single chrom in a single bam file
# input:
#     $BAM_FILE     as arg1
#     $CHROM        as arg2
#     $FAI_FILE     as arg3
#     $BIN_SIZE
#     $MIN_MAPQ
# outputs:
#     bin counts for all possible chrom bins on STDOUT

# get arguments
export BAM_FILE=$1
export CHROM=$2
export FAI_FILE=$3

# get the maximum possible bin number from the chromosome size
# used by count_chrom_bins.pl
CHROM_SIZE=`awk '$1=="'$CHROM'"' $FAI_FILE | cut -f2`
export MAX_BIN_I0=$((($CHROM_SIZE - 1) / $BIN_SIZE))

# pull all distinct insert 5' endpoints from the bam file
bash $ACTION_DIR/../parse_unique_inserts.sh | 

# count insert endpoints per bin
perl $ACTION_DIR/count_chrom_bins.pl
