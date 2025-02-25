# limit values
gcLimits <- c(0.25, 0.7) # optimal for hg38 genome, only used by normalizeGC step (other inherit from collate/scores)

# score types metadata
scoreTypes <- list(
    genome = list(
        gc = list(
            label = "GC",
            unit = "Percent",
            trackHeaderLabel = "GC Percent",
            trackSummaryLabel = "Z Score",
            class = "baseComposition",
            valueLim = gcLimits,
            normValue = "z",
            summaryType = "z"
        )
    ),
    sample = list(
        # gcrz = list(
        #     distUnit = 0.2,
        #     include = c("quantile"),
        #     log10 = FALSE,
        #     label = "GC Residual",
        #     unit = "Z Score",
        #     trackHeaderLabel = "GC Residual Read Count",
        #     trackScoreLabel = "Z Score",
        #     class = "coverage",
        #     valueLim = c(-3.5, 3.5),
        #     deltaLim = c(-3.5, 3.5),
        #     normValue = "z",
        #     summaryType = "quantile"
        # ),
        # gcrz_delta = list(
        #     distUnit = 0.2,
        #     include = c("quantile"),
        #     log10 = FALSE,
        #     label = "GC Residual Delta",
        #     unit = "Z Score",
        #     trackHeaderLabel = "GC Residual Delta",
        #     trackScoreLabel = "Z Score Delta",
        #     class = "coverage",
        #     valueLim = c(-3.5, 3.5),
        #     deltaLim = c(-3.5, 3.5),
        #     normValue = "z",
        #     summaryType = "quantile"
        # ),
        hmmz = list(
            distUnit = 0.2,
            include = c("quantile"),
            label = "Enrichment11",
            unit = "Z Score11",
            trackHeaderLabel = "MiDAS Enrichment",
            trackScoreLabel = "Z Score",
            class = "coverage",
            valueLim = c(-3.5, 3.5),
            deltaLim = c(-3.5, 3.5),
            normValue = "z",
            summaryType = "quantile"
        ),
        hmmz_delta = list(
            distUnit = 0.2,
            include = c("quantile"),
            label = "MiDAS Enrichment Delta",
            unit = "Z Score",
            trackHeaderLabel = "MiDAS Enrichment Delta",
            trackScoreLabel = "Delta Z",
            class = "coverage",
            valueLim = c(-3.5, 3.5),
            deltaLim = c(-3.5, 3.5),
            normValue = "z",
            summaryType = "quantile"
        )
        # ,
        # nrll = list(
        #     distUnit = 0.1,
        #     label = "Protamine Transition Enrichment",
        #     unit = "NRLL",
        #     trackHeaderLabel = "Protamine Transition Enrichment",
        #     trackScoreLabel = "NRLL",
        #     class = "insertSize",
        #     valueLim = c(-1.5, 1.5),
        #     deltaLim = c(-1.5, 0),
        #     normValue = "z",
        #     summaryType = "quantile"
        # )
    )
)

geneTargets <- list(
    brca2 = c("brca2_ctl","brca2_shrna"), 
    rad51 = c("rad51_ctl","rad51_sirna")
)
cellLines <- list(
    brca2 = "H1299",
    rad51 = "U2OS"
)
revGeneTargets <- list(
    brca2_ctl   = "brca2",
    brca2_shrna = "brca2",
    rad51_ctl   = "rad51",
    rad51_sirna = "rad51"
)
getTrainingRegions <- function(sample_name){
    geneTarget <- revGeneTargets[[sample_name]]
    fread(file.path( # empirically determined regions representing CN states in control sample
        serverEnv$MDI_DIR, 
        "suites/definitive/erfs-mdi-tools/shiny/apps/erfs_vis/modules/appSteps/erfs_copyNumber/cn_training_ctl.csv"
    ))[gene_target == geneTarget]
}
getTrainingRegions_all <- function(sample_name){
    geneTarget <- revGeneTargets[[sample_name]]
    fread(file.path( # empirically determined regions representing CN states regardless of stress, i.e., in all samples
        serverEnv$MDI_DIR, 
        "suites/definitive/erfs-mdi-tools/shiny/apps/erfs_vis/modules/appSteps/erfs_copyNumber/cn_training_all.csv"
    ))[gene_target == geneTarget]
}
getCnOptimizedFits <- function(sample_name_){
    fread(file.path( # empirically determined regions representing CN states regardless of stress, i.e., in all samples
        serverEnv$MDI_DIR, 
        "suites/definitive/erfs-mdi-tools/shiny/apps/erfs_vis/modules/appSteps/erfs_copyNumber/cn_count_optim.csv"
    ))[sample_name == sample_name_]
}
traceColors <- list(
    sample = c(
        brca2_ctl   = CONSTANTS$plotlyColors$green,
        brca2_shrna = CONSTANTS$plotlyColors$blue,
        rad51_ctl   = CONSTANTS$plotlyColors$orange,
        rad51_sirna = CONSTANTS$plotlyColors$purple
    ),
    gene_target = c(
        brca2       = CONSTANTS$plotlyColors$blue,
        rad51       = CONSTANTS$plotlyColors$purple
    )
)
cnColors = c(
    CONSTANTS$plotlyColors$green,
    CONSTANTS$plotlyColors$blue,
    CONSTANTS$plotlyColors$orange,
    CONSTANTS$plotlyColors$purple,
    CONSTANTS$plotlyColors$red,
    "grey20"
)

# standardized color palettes
hmColors <- list(
    # BROWN  = rgb(0.2, 0,   0),   # extreme gain    
    # RED    = rgb(0.9, 0,   0),   # 3, CN3, full gain
    # YELLOW = rgb(0.8, 0.8, 0),   # 2.5, mosaic state
    # GREEN  = rgb(0,   0.8, 0),   # 2, CN2, no CNV
    # CYAN   = rgb(0,   0.6, 1),   # 1.5, mosaic state
    # BLUE   = rgb(0,   0,   1),   # 1, CN1, full loss
    # PURPLE = rgb(0.4, 0,   0.4), # extreme loss
    # BLACK  = rgb(0,   0,   0)    # 0, nothing    
    RED     = rgb(0.9, 0,   0),
    GREY    = rgb(0.75, 0.75, 0.75),
    BLUE    = rgb(0,   0,   1)
)
sampleColors <- c(
    CONSTANTS$plotlyColors$black,
    CONSTANTS$plotlyColors$blue,
    CONSTANTS$plotlyColors$orange,
    CONSTANTS$plotlyColors$green,
    CONSTANTS$plotlyColors$purple
)
nTrackMapColorsPerSide <- 30
trackMapColors <- list(
    low  = colorRampPalette(c(hmColors$GREY, hmColors$BLUE))(nTrackMapColorsPerSide + 1), # blue color is cold/depleted,
    high = colorRampPalette(c(hmColors$GREY, hmColors$RED))( nTrackMapColorsPerSide + 1)  # red  color is hot/ enriched
)
