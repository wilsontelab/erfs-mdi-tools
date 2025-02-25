# functions to load, parse, and filter bin data

# handle bin exclusions, including at gc extremes
isExcludedBin <- function(excluded, pct_gc){
    excluded == 1 | pct_gc < gcLimits[1] | pct_gc > gcLimits[2]
}
isIncludedAutosomeBin <- function(excluded, pct_gc, nAlleles){
    !isExcludedBin(excluded, pct_gc) & nAlleles == 2
}
getIncudedAutosomeBins <- function(bins){
    isIncludedAutosomeBin(bins$excluded, bins$pct_gc, 2) # , bins$nAlleles
}

# load and format genome bins and read counts (from atat/collate action)
erfsBinData <- function(sourceId){
    startSpinner(session, message = "loading collated bin data")
    filePath <- loadPersistentFile(
        sourceId = sourceId, 
        contentFileType = "binCounts", 
        ttl = CONSTANTS$ttl$month, 
        postProcess = function(bd){
            bd$binCounts$binI <- 1:nrow(bd$binCounts)
            bd
        }
    )
    stopSpinner(session)
    persistentCache[[filePath]]$data
}

# load and format pipeline-external data
erfsCallData <- function(sourceId){
    startSpinner(session, message = "loading external call data")
    filePath <- loadPersistentFile(
        sourceId = sourceId, 
        contentFileType = "externalCalls", 
        ttl = CONSTANTS$ttl$month,
        postProcess = function(cd){
            for(sample in names(cd$coverage)){
                cd$coverage[[sample]]$binI <- 1:nrow(cd$coverage[[sample]])
            }
            # override the pipeline calls, they are from hg19...
            cd$calls$brca2 <- fread(file.path( # empirically determined regions representing CN states regardless of stress, i.e., in all samples
                serverEnv$MDI_DIR, 
                "suites/definitive/erfs-mdi-tools/shiny/apps/erfs_vis/classes/browserTracks/erfsCalls/Tarsounas_peaks_hg38.bed"
            ))[, 1:3]
            setnames(cd$calls$brca2, c("chrom","start0","end1"))
            cd
        }
    )
    stopSpinner(session)
    persistentCache[[filePath]]$data
}

# Z-Score analysis, depends on pipeline bin data and user-fitted models
#----------------------------------------------------------------------
# analyze and aggregate distributions of different bin scores
# all scores are expected to be one a comparable scale between samples, including
replaceNaN <- function(x){
    x[is.nan(x)] <- NA
    x
}
analyzeScoreDist <- function(bins, scores, scoreType){
    scores_wrk <- scores[getIncudedAutosomeBins(bins)] # only use good autosomal bins to analyze score distributions
    dist <- data.table(x = as.integer(round(scores_wrk / scoreType$distUnit)) * scoreType$distUnit)[, .(y = .N), keyby = .(x)]
    dist[, y := y / replaceNaN(sum(y, na.rm = TRUE))]
    mu <- replaceNaN(mean(scores_wrk, na.rm = TRUE))
    sd <- replaceNaN(sd(  scores_wrk, na.rm = TRUE))
    median <- replaceNaN(median(scores_wrk, na.rm = TRUE))
    list(
        score    = scores, # score are dropped downstream when not included in scoreType
        dist     = dist,
        median   = median,
        mean     = mu,
        sd       = sd,
        z        = if("z" %in% scoreType$include) (scores - median) / sd else NULL,         # for normal/parametric scores
        quantile = if("quantile" %in% scoreType$include) ecdf(scores_wrk)(scores) else NULL # for non-parametric scores
    )
}
analyzeSampleScores <- function(bd, scoreFn, scoreType, ...){
    x <- lapply(names(bd$samples), function(sample_name){
        startSpinner(session, message = paste("loading hmmz", sample_name))
        analyzeScoreDist(bd$binCounts, scoreFn(bd, sample_name, ...), scoreType)
    })
    names(x) <- names(bd$samples)
    x
}
aggregateAndAnalyzeScores <- function(bins, sampleScores, sample_names, scoreType){
    x <- as.data.table(sapply(sample_names, function(sample_name){
        replaceNaN(sampleScores[[sample_name]]$score)
    }, simplify = FALSE, USE.NAMES = TRUE))
    analyzeScoreDist(bins, x[[2]] - x[[1]], scoreType)
}
aggregateSampleScores <- function(bd, sampleScores, scoreType){
    x <- lapply(names(geneTargets), function(geneTarget_){
        startSpinner(session, message = paste("loading hmmz delta", geneTarget_))
        sample_names <- geneTargets[[geneTarget_]]
        aggregateAndAnalyzeScores(bd$binCounts, sampleScores, sample_names, scoreType)
    })
    names(x) <- names(geneTargets)
    x
}
hmmZScores <- function(sourceId, gcBiasModels = NULL){
    if(is.null(gcBiasModels)) gcBiasModels <- app$normalizeGC$getGcBiasModels_externalCall(sourceId)
    hmmzCache$get(
        "gcResiduals",
        keyObject = list(
            sourceId = sourceId,
            gcBiasModels = gcBiasModels
        ),
        permanent = TRUE,
        from = "ram",
        create = "asNeeded", #"asNeeded",
        createFn = function(...){
            startSpinner(session, message = "clearing hmmz")
            hmmzCache$clear()
            startSpinner(session, message = "loading hmmz")
            bd <- erfsBinData(sourceId)
            scoreTypeName <- "hmmz"
            scoreType <- scoreTypes$sample[[scoreTypeName]]
            sampleScores <- analyzeSampleScores(bd, function(bd, sample_name){
                optim <- getCnOptimizedFits(sample_name)
                binCounts <- bd$binCounts[[sample_name]]
                ref_sample_name <- sub("_shrna", "_ctl", sample_name)
                ref_sample_name <- sub("_sirna", "_ctl", ref_sample_name)
                cnModel <- gcBiasModels[[ref_sample_name]]$hmm$cn
                fit <- merge(
                    data.table(copy_number = cnModel),
                    optim[, .(copy_number, mu, theta)],
                    by = "copy_number",
                    all.x = TRUE,
                    sort = FALSE
                )
                x <- data.table(
                    a = pnbinom(binCounts + 1L, size = fit$theta, mu = fit$mu),
                    b = pnbinom(binCounts,      size = fit$theta, mu = fit$mu)
                )
                x <- qnorm(x[, rowMeans(.SD)])
                x[abs(x) == Inf] <- 100
                x
            }, scoreType)
            aggregateScores <- aggregateSampleScores(bd, sampleScores, scoreType)
            x <- list(
                sampleScores    = sampleScores,
                aggregateScores = aggregateScores
            )
            stopSpinner(session)
            x
        }
    )$value
}

# functions to load, parse, and filter bin data
getScoreLevel <- function(scoreTypeName){
         if(scoreTypeName %in% names(scoreTypes$genome)) 'genome'
    else if(scoreTypeName %in% names(scoreTypes$sample)) 'sample'
    else "NA"
}
getScoreType <- function(sourceId, scoreTypeName){
    scoreLevel <- getScoreLevel(scoreTypeName)
    scoreTypes[[scoreLevel]][[scoreTypeName]]
}
# getTypedStages <- function(sourceId) unlist(erfsScores(sourceId)$stageTypes)
# getStageTypesByStage <- function(sourceId, stages) unlist(erfsScores(sourceId)$reverseStageTypes[stages])

#----------------------------------------------------------------------
# sample-level score structure summary and associated score object retrieval
#----------------------------------------------------------------------
getGenomeScores <- function(sourceId, scoreTypeName){ # returns a single genome-level score object
    x <- list(erfsScores(sourceId)$scores$genome[[scoreTypeName]])
    names(x) <- scoreTypeName
    x
}
getSampleScoresList <- function(sourceId, scoreTypeName){ # returns a list of sample-level score objects based on GC normalization
    # if(scoreTypeName == "hmmz"){
        # gcResiduals(sourceId)
        hmmZScores(sourceId)
    # } else {
    #     erfsScores(sourceId)$scores$sample[[scoreTypeName]]
    # }
}
getSampleScores <- function(sourceId, scoreTypeName){ # returns a list of sample-level score objects
    getSampleScoresList(sourceId, scoreTypeName)$sampleScores  # [samples]
}
getSampleDeltaScores <- function(sourceId, scoreTypeName){ # returns a list of stage-level score objects matching a list of samples
    getSampleScoresList(sourceId, scoreTypeName)$aggregateScores # $by_stage[unique(samples$stage)]
}
# getStageTypeScores <- function(sourceId, scoreTypeName, samples){ # returns a list of stageType-level score objects matching a list of samples
#     stageTypes <- getStageTypesByStage(sourceId, samples$stage)
#     getSampleScoresList(sourceId, scoreTypeName)$aggregateScores$by_stageType[unique(stageTypes)]
# }
# getStageTypeDeltaScores <- function(sourceId, scoreTypeName, cleanDist = FALSE){ # returns a single score object for the delta between stage types (or a genome-level score object)
#     list(
#         stageType_delta = if(getScoreLevel(scoreTypeName) == "sample") {
#             getSampleScoresList(sourceId, scoreTypeName)$aggregateScores$stageType_delta
#         } else {
#             x <- erfsScores(sourceId)$scores$genome[[scoreTypeName]]
#             if(scoreTypeName == "txn" && cleanDist) x$dist <- x$dist[x %between% scoreTypes$genome$txn$valueLim]
#             x
#         }
#     )
# }
getSeriesAggScores <- function(sourceId, scoreTypeName){
    if(endsWith(scoreTypeName, "_delta")){
        getSampleDeltaScores(sourceId, scoreTypeName)
    } else {
        getSampleScores(sourceId, scoreTypeName)
    }
    # switch(
    #     config$Aggregate_By,
    #     sample      = getSampleScores(sourceId, scoreTypeName, samples),
    #     stage       = getStageScores(sourceId, scoreTypeName, samples),
    #     stage_type  = getStageTypeScores(sourceId, scoreTypeName, samples)
    # )
}

getWindowBinData <- function(sourceId, bd, sample_name, binI){
    ref_sample_name <- sub("_shrna", "_ctl", sample_name)
    ref_sample_name <- sub("_sirna", "_ctl", ref_sample_name)
    gcBiasModels <- app$normalizeGC$getGcBiasModels_externalCall(sourceId)
    cnModel <- gcBiasModels[[ref_sample_name]]$hmm$cn[binI]
    optim <- getCnOptimizedFits(sample_name)
    fit <- merge(
        data.table(
            start0      = bd$binCounts[binI, start0],
            end1        = bd$binCounts[binI, end1],
            bin_count   = bd$binCounts[binI][[sample_name]],
            copy_number = cnModel
        ),
        optim[, .(
            copy_number, 
            mu, 
            theta
        )],
        by = "copy_number",
        all.x = TRUE,
    )
}

