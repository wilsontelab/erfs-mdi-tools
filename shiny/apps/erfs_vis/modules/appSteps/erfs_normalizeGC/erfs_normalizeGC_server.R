#----------------------------------------------------------------------
# server components for the erfs_normalizeGC appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
erfs_normalizeGCServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'erfs_normalizeGC'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# data package sources and source-level data objects derived from pipeline
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("source", selection = "single")
samples <- reactive({ # vector of the names of all co-analyzed samples
    sourceId <- sourceId()
    req(sourceId)
    data.table(sample_name = names(erfsBinData(sourceId)$samples))
})

#----------------------------------------------------------------------
# samples table
#----------------------------------------------------------------------
sampleTable <- bufferedTableServer(
    "sample",
    id,
    input,
    tableData = samples,
    selection = 'single',
    options = list(
        paging = FALSE,
        searching = FALSE  
    )
)
sample <- reactive({
    row <- sampleTable$rows_selected()
    req(row)
    samples()[row]
})

#----------------------------------------------------------------------
# appStep outcomes, saved to disk since these are ~one-time analysis steps
# includes GC bias fit, chromosome-level data and junction fits, but not HMM
#----------------------------------------------------------------------
gcBiasFileName <- "gcBiasModels.rds"
gcBiasFile <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    expandSourceFilePath(sourceId, gcBiasFileName)
})
invalidateGcBiasModels <- reactiveVal(1)
getGcBiasModels <- function(gcBiasFile = NULL, sourceId = NULL){
    if(is.null(gcBiasFile)) gcBiasFile <- expandSourceFilePath(sourceId, gcBiasFileName)
    if(file.exists(gcBiasFile)) readRDS(gcBiasFile) else list()
}
gcBiasModels <- reactive({ # gc bias model from negative binomial are calculated synchronously
    invalidateGcBiasModels()
    getGcBiasModels( gcBiasFile = gcBiasFile() )
})
gcBiasModel <- reactive({
    gcBiasModels <- gcBiasModels()
    sample <- sample()
    gcBiasModels[[sample$sample_name]]
})
cnModel <- reactive({
    sample_name <- sample()$sample_name
    sample_name <- sub("_shrna", "_ctl", sample_name)
    sample_name <- sub("_sirna", "_ctl", sample_name)
    gcBiasModels <- gcBiasModels()
    gcBiasModels[[sample_name]]$hmm$cn
})

#----------------------------------------------------------------------
# interactive GC bias plots, selection cascades to solving negative binomial
#----------------------------------------------------------------------
gcBiasPlotData <- function(){
    sourceId <- sourceId()
    sample <- sample()
    req(sourceId, sample)
    startSpinner(session, message = paste("plotting", sample$sample_name))
    bd <- erfsBinData(sourceId)

    # # incorrectly assumes CN = 2 for all bins, use to help find empirical training regions
    # I <- bd$binCounts[, excluded == 0 & chrom != "chrX" & chrom != "chrY"]
    # data.table(
    #     x = bd$binCounts[I, pct_gc], # same as fractionGC
    #     y = bd$binCounts[I][[sample$sample_name]] / 2, # rpba = reads per bin per allele
    #     nAlleles = 2, # same as binCN
    #     color = CONSTANTS$plotlyColors$blue
    # )[sample.int(.N, min(.N, settings$get("GC_Bias","N_Plotted_Bins")))]

    trainingRegions <- getTrainingRegions(sample$sample_name)
    do.call(rbind, lapply(1:length(cnColors), function(CN){
        trs <- trainingRegions[copy_number == CN]
        do.call(rbind, lapply(1:nrow(trs),function(i){ 
            I <- bd$binCounts[, 
                excluded == 0 & 
                chrom  == trs[i, chrom] &
                start0 >= trs[i, start0] &
                end1   <= trs[i, end1]
            ]
            data.table(
                x = bd$binCounts[I, pct_gc], # same as fractionGC
                y = bd$binCounts[I][[sample$sample_name]] / CN, # rpba = reads per bin per allele
                nAlleles = CN, # same as binCN
                color = cnColors[CN]
            )
        }))
    }))[sample.int(.N, min(.N, settings$get("GC_Bias","N_Plotted_Bins")))]
}
gcOverplotData <- function(){
    gcBiasModel <- gcBiasModel()
    if(!isTruthy(gcBiasModel)) return(NULL)
    sample <- sample()
    startSpinner(session, message = paste("overplotting", sample$sample_name))
    nb <- gcBiasModel$fit
    gc <- nb$model$fractionGC
    data.table(
        x = gc,
        y = predict(nb, gc, type = 'mu') # rpba = reads per bin per allele * nAlleles = rpb
    )
}
gcBiasPlot <- interactiveScatterplotServer(
    "gcBiasPlot",
    plotData = reactive({ 
        x <- gcBiasPlotData() 
        stopSpinner(session)
        x
    }),
    accelerate = TRUE,
    overplot = reactive({
        x <- gcOverplotData()
        stopSpinner(session)
        x
    }),
    overplotMode = "lines",
    overplotColor = CONSTANTS$plotlyColors$red,
    xtitle = "Fraction GC",
    xrange = gcLimits,
    ytitle = "Reads Per Bin Per Allele",
    yrange = function(...) range_pos(..., foldIQR = 5),
    selectable = "lasso"
)
observeEvent(gcBiasPlot$selected(), {
    showUserDialog(
        "Use this GC Bias Fit?", 
        tags$p(
            "If you are happy with your GC bias selection, click OK to run the next, slow actions."
        ), 
        tags$p(
            "If you are not happy, click Cancel and repeat the GC selection."
        ), 
        callback = function(parentInput) {
            removeModal()
            fitGCBiasFromSelected(gcBiasPlot$selected())
        },
        type = 'okCancel', 
        easyClose = FALSE, 
        fade = 250
    )
})

#----------------------------------------------------------------------
# fit negative binomial distribution to GC bias
#----------------------------------------------------------------------
fitGCBiasFromSelected <- function(selected){
    req(selected, nrow(selected) > 10)
    sample <- sample()
    startSpinner(session, message = paste("fitting", sample$sample_name))
    selected <- as.data.table(selected)

    # fit the negative binomial distribution to the selected points
    gcBiasModels <- gcBiasModels()
    d <- gcBiasPlotData() 
    CN <- d[selected$pointNumber + 1, nAlleles]
    fit <- new_nbinomCountsGC(
        binCounts  = selected$y * CN,
        fractionGC = selected$x,
        binCN      = CN,
        method = 'cubic'
    )

    # run the HMM on all bins using the fitted model
    sourceId <- sourceId()
    bd <- erfsBinData(sourceId)
    sample <- sample()
    hmm <- viterbi(
        fit, 
        bd$binCounts[[sample$sample_name]],
        bd$binCounts$pct_gc,
        maxCN = 6, 
        transProb = 1e-7, # options for the HMM
        chroms = bd$binCounts$chrom, # if a vector, use keyedViterbi by chromosome
        forceCNs = NULL, # if a vector, force to HMM output for bins that are not NA
        asRle = FALSE
    )
    message("hmm done")

    # save for future use
    gcBiasModels[[sample$sample_name]] <- list(
        sample = sample,
        fit    = fit,
        hmm    = hmm
    )
    saveRDS(gcBiasModels, file = gcBiasFile())
    stopSpinner(session)
    invalidateGcBiasModels( invalidateGcBiasModels() + 1 ) 
}
# List of 3
#  $ hmm  :List of 3
#   ..$ emissProbs: num [1:308837, 1:7] -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0
# .1 -0.1 ...
#   ..$ transProbs: num [1:7, 1:7] -6.00e-07 -1.61e+01 -1.61e+01 -1.61e+01 -1.61e+
# 01 ...
#   ..$ keys      :List of 2
#   .. ..$ lengths: int [1:24] 24896 24220 19830 19022 18154 17081 15935 14514 138
# 40 13380 ...
#   .. ..$ values : chr [1:24] "chr1" "chr2" "chr3" "chr4" ...
#   .. ..- attr(*, "class")= chr "rle"
#   ..- attr(*, "class")= chr "hmmEPTable"
#  $ cn   : int [1:308837] 0 0 0 0 0 0 0 0 0 0 ...
#  $ maxCN: num 6

#----------------------------------------------------------------------
# interactive GC bias plots, selection cascades to solving negative binomial
#----------------------------------------------------------------------
trainingBinCounts <- reactive({
    sourceId <- sourceId()
    sample <- sample()
    req(sourceId, sample)
    bd <- erfsBinData(sourceId)
    cnModel <- cnModel()
    trainingRegions <- getTrainingRegions_all(sample$sample_name)
    do.call(rbind, lapply(1:nrow(trainingRegions),function(i){ 
        I <- bd$binCounts[,
            excluded == 0 & 
            chrom  == trainingRegions[i, chrom] &
            start0 >= trainingRegions[i, start0] &
            end1   <= trainingRegions[i, end1]
        ]
        data.table(
            binCount = bd$binCounts[I][[sample$sample_name]],
            binCnHmm = cnModel[I]
        )
    }))
})
binCountsPlotData <- function(){
    tr_bins <- trainingBinCounts()
    do.call(rbind, lapply(1:length(cnColors), function(CN_hmm){
        agg <- tr_bins[binCnHmm == CN_hmm, .N, by = .(binCount)]
        data.table(
            x = agg$binCount,
            y = agg$N / sum(agg$N, na.rm = TRUE),
            color = cnColors[CN_hmm]
        )  
    }))
}
binCountsPlot <- staticPlotBoxServer(
    "binCountsPlot",
    maxHeight = "400px",
    points   = TRUE,
    legend  = TRUE,
    margins = FALSE,
    create = function() {
        d <- binCountsPlotData()
        par(mar = c(4.1, 4.1, 0.1, 0.1))
        maxX <- 85
        binCountsPlot$initializeFrame(
            xlim = c(0, maxX),
            ylim = c(0, max(d$y) * 1.05),
            xlab = "Bin Count",
            ylab = "Frequency"
        )
        x <- 0:maxX
        for (CN in 1:5){
            mu <- input[[paste0("mu", CN)]]
            theta <- input[[paste0("theta", CN)]]
            abline(v = mu, col = cnColors[CN])
            x <- 0:maxX
            lines(x, dnbinom(x, size = theta, mu = mu), col = cnColors[CN])
        }
        binCountsPlot$addPoints(
            x = d$x,
            y = d$y,
            col = d$color
        )
    }
)

#----------------------------------------------------------------------
# external utilities for browser support, etc.
#----------------------------------------------------------------------
getGcBiasModels_externalCall <- function(sourceId){
    gcSourceId <- tryCatch(sourceId(), error = function(e) NULL)
    if(isTruthy(gcSourceId) && gcSourceId == sourceId) gcBiasModels()
    else getGcBiasModels(sourceId = sourceId)
}
getBinZScore <- function(sourceId, sampleName, binCounts, fractionGC, binCN){
    gcBiasModels <- getGcBiasModels_externalCall(sourceId)
    gcBiasModel <- gcBiasModels[[sampleName]]
    if(!isTruthy(gcBiasModel)) return(NULL)
    zScore(gcBiasModel$fit, binCounts, fractionGC, binCN)
}

# #----------------------------------------------------------------------
# # download GC residuals
# #----------------------------------------------------------------------
# output$downloadGcResiduals <- downloadHandler(
#     filename = function(){
#         paste(module, "gcrz.rds", sep = ".")
#     },
#     content  = function(tmpFile){
#         sourceId <- sourceId()
#         req(sourceId)
#         startSpinner(session, message = "downloading GC residuals")
#         x <- getSampleScoresList(sourceId, "gcrz")
#         dstr(x)
#         saveRDS(x, tmpFile)
#         stopSpinner(session)
#     }
# )

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    outcomes = list(),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    getGcBiasModels_externalCall = getGcBiasModels_externalCall,
    getBinZScore = getBinZScore,
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
