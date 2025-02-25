#----------------------------------------------------------------------
# server components for the erfs_copyNumber appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
erfs_copyNumberServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'erfs_copyNumber'
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
# interactive GC bias plots, selection cascades to solving negative binomial
#----------------------------------------------------------------------
binCountsPlotData <- function(){
    sourceId <- sourceId()
    sample <- sample()
    req(sourceId, sample)
    startSpinner(session, message = paste("plotting", sample$sample_name))
    bd <- erfsBinData(sourceId)
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
            agg <- data.table(binCount = bd$binCounts[I][[sample$sample_name]])[,
                .N,
                by = .(binCount)
            ]
            data.table(
                x = agg$binCount,
                y = agg$N / sum(agg$N, na.rm = TRUE),
                color = cnColors[CN]
            )   
        }))
    }))
}
binCountsPlot <- interactiveScatterplotServer(
    "binCountsPlot",
    plotData = reactive({ 
        x <- binCountsPlotData() 
        stopSpinner(session)
        x
    }),
    pointSize = 6,
    accelerate = TRUE,
    xtitle = "Bin Count",
    xrange = c(0,60),
    ytitle = "Frequency",
)

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
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
