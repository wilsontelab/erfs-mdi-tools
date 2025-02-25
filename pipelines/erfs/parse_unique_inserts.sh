# action:
#     prepare bam alignments for bin counting and insert size assessment
# expects:
#     $EXCLUSIONS_BED and $GAPS_FILE
#     $BAM_FILE
#     $CHROM, can be empty to assess all chromosomes
#     $BIN_SIZE
#     $MIN_MAPQ
# outputs:
#     distinct insert endpoints in format "start1\tstrand_flag_1[,strand_flag_2]" on STDOUT

# pull the single read from high-quality alignments
#   --exclude-flags 3844
#       read unmapped (0x4)
#       not primary alignment (0x100)
#       read fails platform/vendor quality checks (0x200)
#       read is PCR or optical duplicate (0x400)
#       supplementary alignment (0x800)
# enforce minimum mapping quality
samtools view --bam --exclude-flags 3844 --min-MQ $MIN_MAPQ $BAM_FILE $CHROM | 

# remove excluded regions as requested
bedtools intersect -v -a - -b ${EXCLUSIONS_BED} |
samtools view |

# parse to start1,strand_flag
awk 'BEGIN{OFS="\t"}{print $4, and($2, 16)}' | 

# ensure than only unique reads pairs are counted, i.e., finish de-duplicating reads
# can have the same read in final output if found in multiple input files
bedtools groupby -g 1 -c 2 -o distinct
