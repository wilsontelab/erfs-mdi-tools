#----------------------------------------------------------------------
# UI components for the erfs_copyNumber appStep module
#----------------------------------------------------------------------

# module ui function
erfs_copyNumberUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$erfs_copyNumber)

    # UI functions
    plotBox_ <- function(title, ui){
        box(
            title = title,
            width = 6,
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE,
            ui
        )
    }

    # return the UI contents
    standardSequentialTabItem(

        # page header text
        options$longLabel,
        options$leaderText,

        # page header links, uncomment as needed
        id = id,
        # documentation = TRUE,
        # terminal = TRUE,
        console = serverEnv$IS_DEVELOPER,
        code = serverEnv$IS_DEVELOPER,
        settings = TRUE,

        # data source selectors
        fluidRow(
            dataSourceTableUI(
                ns("source"), 
                "Data Source", 
                width = 6, 
                collapsible = FALSE,
                inFluidRow = FALSE
            ),
            bufferedTableUI(
                ns("sample"),
                "Sample",
                width = 6,
                solidHeader = TRUE,
                status = "primary",
                collapsible = FALSE
            )
        ),

        # GC bias plots
        fluidRow(
            plotBox_(
                "Bin Counts",
                interactiveScatterplotUI(ns("binCountsPlot"), height = '400px')
            )
            
        ),
        # fluidRow(
        #     staticPlotBoxUI(
        #         ns("gcBiasFitComposite"), 
        #         "Composite of GC Bias Fits",
        #         width = 6,
        #         status = "primary",
        #         collapsible = TRUE,
        #         solidHeader = TRUE,
        #         collapsed = FALSE
        #     ),
        #     column(
        #         width = 6,
        #         downloadLink(ns("downloadGcResiduals"), label = "download GC residual Z-scores")
        #     )
        # ),
        NULL
    )
}
