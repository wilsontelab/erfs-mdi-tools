# action:
#     in parallel over multiple samples:
#         count unique read endpoints in genome bins
# input:
#     ${GENOME_BINS_BED} = bed file of genome bins created by genome/bin pipeline action
#     aligned, sorted, and indexed bam files (not necessarily deduplicated, we will dedup)
# outputs:
#     samples   = list of bam files
#     binCounts = data.table of bin metadata and deduplicated read counts for each sample

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
        'GENOME',
        'GENOME_FASTA',
        'GENOME_BINS_BED',
        'ACTION_DIR',
        'DATA_GENOME_PREFIX',
        'BRCA2_CTL',
        'BRCA2_SHRNA',
        'RAD51_CTL',
        'RAD51_SIRNA'
    ),
    integer = c(
        'MIN_MAPQ',
        'BIN_SIZE',
        'N_CPU'
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

message("loading genome bins")
faiFile <- paste0(env$GENOME_FASTA, '.fai')
bins    <- fread(env$GENOME_BINS_BED) # already restricted by genome/bin to autosomes and chrX/Y
chroms  <- bins[, unique(chrom)]
n_bins  <- nrow(bins)
# bins[, nAlleles := ifelse(chrom %in% c('chrX', 'chrY'), 1L, 2L)]

message("calculating sample bin counts")
countChromBins <- paste('bash', file.path(env$ACTION_DIR, 'count_chrom_bins.sh'))
samples <- list(
    brca2_ctl   = env$BRCA2_CTL, 
    brca2_shrna = env$BRCA2_SHRNA, 
    rad51_ctl   = env$RAD51_CTL, 
    rad51_sirna = env$RAD51_SIRNA
)
nSamples <- length(samples)
binCounts <- do.call(cbind, mclapply(1:nSamples, function(sampleI) {
# binCounts <- do.call(cbind, lapply(1:nSamples, function(sampleI) {
    sampleName <- names(samples)[sampleI]
    message(paste('   ', sampleName))
    bamFile <- samples[[sampleName]]
    unlist(lapply(chroms, function(chrom) { # all bins values over all ordered chroms
        fread(cmd = paste(countChromBins, bamFile, chrom, faiFile)) # one value per bin on chrom
    }))
}, mc.cores = env$N_CPU))
# }))

binCounts <- cbind(bins, binCounts)
setnames(binCounts, c('chrom', 'start0', 'end1', 'excluded', 'pct_gc', names(samples)))

message()
message("saving binCounts for app")
obj <- list(
    bin_size    = env$BIN_SIZE,
    samples     = samples,
    chroms      = chroms,
    binCounts   = binCounts
)
saveRDS(
    obj, 
    file = paste(env$DATA_GENOME_PREFIX, "binCounts.rds", sep = '.')
)
str(obj)
#=====================================================================================
