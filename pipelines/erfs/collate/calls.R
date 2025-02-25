# action:
#     assemble external call data

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(parallel)
}))
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'DATA_GENOME_PREFIX',
        'BRCA2_CALLS',
        'RAD51_CALLS',
        'BRCA2_CTL_COV',
        'BRCA2_SHRNA_COV',
        'RAD51_CTL_COV',
        'RAD51_SIRNA_COV'
    )
))
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# loop through all BAM files to determine their coverage in each genome bin
#-------------------------------------------------------------------------------------

message("loading external calls")

BRCA2_CALLS <- fread(env$BRCA2_CALLS)
setnames(BRCA2_CALLS, c("chrom","start0","end1","genes"))

RAD51_CALLS <- fread(env$RAD51_CALLS)
setnames(RAD51_CALLS, c("chrom","start0","end1","n_bins","call_density"))

BRCA2_CTL_COV <- fread(env$BRCA2_CTL_COV)
setnames(BRCA2_CTL_COV, c("chrom","start0","end1","coverage","bp_at_max","bin_size","fraction_at_max"))

BRCA2_SHRNA_COV <- fread(env$BRCA2_SHRNA_COV)
setnames(BRCA2_SHRNA_COV, c("chrom","start0","end1","coverage","bp_at_max","bin_size","fraction_at_max"))

RAD51_CTL_COV <- fread(env$RAD51_CTL_COV)
setnames(RAD51_CTL_COV, c("chrom","start0","end1","coverage","bp_at_max","bin_size","fraction_at_max"))

RAD51_SIRNA_COV <- fread(env$RAD51_SIRNA_COV)
setnames(RAD51_SIRNA_COV, c("chrom","start0","end1","coverage","bp_at_max","bin_size","fraction_at_max"))

message("saving output for app")
obj <- list(
    calls = list(
        brca2 = BRCA2_CALLS[, .(chrom, start0, end1)],
        rad51 = RAD51_CALLS[, .(chrom, start0, end1)]
    ),
    coverage = list(
        brca2_ctl   = BRCA2_CTL_COV[,   .(chrom, start0, end1, coverage)],
        brca2_shrna = BRCA2_SHRNA_COV[, .(chrom, start0, end1, coverage)],
        rad51_ctl   = RAD51_CTL_COV[,   .(chrom, start0, end1, coverage)],
        rad51_sirna = RAD51_SIRNA_COV[, .(chrom, start0, end1, coverage)]
    )
)
saveRDS(
    obj, 
    file = paste(env$DATA_GENOME_PREFIX, "calls.rds", sep = '.')
)
str(obj)
#=====================================================================================
