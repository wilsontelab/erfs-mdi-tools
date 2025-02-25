# establish the dynamic color palettes for heat map and other visualizations

# functions that convert different types of normalized scores to color ranges
z_score_color <- function(zScore, maxZScore){
    minZScore <- -maxZScore 
    z <- pmax(minZScore, pmin(maxZScore, zScore))
    I <- floor(nTrackMapColorsPerSide * abs(z) / maxZScore) + 1L
    ifelse(z < 0, trackMapColors$low[I], trackMapColors$high[I])
}
quantile_score_color <- function(quantile, maxQuantile){
    minQuantile <- 1 - maxQuantile
    quantile <- pmax(minQuantile, pmin(maxQuantile, quantile))
    I <- floor(nTrackMapColorsPerSide * abs(quantile - 0.5) / (maxQuantile - 0.5)) + 1L
    ifelse(quantile < 0.5, trackMapColors$low[I], trackMapColors$high[I])
}

# functions to get distribution trace colors
getSampleColorsByStage <- function(allSamples, samples){
    stages <- allSamples[, unique(stage)]
    colors <- sapply(allSamples$stage, function(x) stageColors[which(stages == x)])
    names(colors) <- allSamples$sample_name
    colors[samples$sample_name]
}
getStageColors <- function(allSamples, samples){
    stages <- allSamples[, unique(stage)]
    colors <- stageColors[1:length(stages)]
    names(colors) <- stages
    colors[samples[, unique(stage)]]
}
getStageTypeColors <- function(sourceId, allSamples, samples){
    allStageTypes <- unique(getStageTypesByStage(sourceId, allSamples[, unique(stage)]))
    stageTypes    <- unique(getStageTypesByStage(sourceId,    samples[, unique(stage)]))
    colors <- stageTypeColors[1:length(allStageTypes)]
    names(colors) <- allStageTypes
    colors[stageTypes]
}
getSeriesSummaryColors <- function(scoreTypeName, scoreValues, config){ # used to color the top group-level summary row of every track heatmap group
    switch(
        scoreTypeName,
        gc   = z_score_color(       scoreValues, config$Max_Z_Score),
        txn  = cpm_score_color(     scoreValues, config$Min_Txn_Log10_CPM, config$Max_Txn_Log10_CPM),
        hmmz = quantile_score_color(scoreValues, config$Max_Quantile)
    )
}
getSeriesSampleColors <- function(scoreTypeName, scoreValues, config){ # used to color the subsequent sample-level rows of every track heatmap group
    switch(
        scoreTypeName,
        hmmz = z_score_color(scoreValues,           config$Max_Z_Score),
        hmmz_delta = z_score_color(scoreValues,           config$Max_Z_Score)
    )
}
getSeriesQuantileColors <- function(scoreTypeName, scoreValues, config){ # like above, now coloring the quantile rows for intra-sample/group relative scores
    quantile_score_color(scoreValues, config$Max_Quantile)
}
