#----------------------------------------------------------------------
# UI components for the erfs_scoreSummary appStep module
#----------------------------------------------------------------------

# module ui function
erfs_scoreSummaryUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$erfs_scoreSummary)

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
        settings = FALSE,

        # appStep UI elements, populate as needed
        fluidRow(
            dataSourceTableUI(
                ns("dataSourceTable"),
                "Data Source", 
                width = 6, 
                collapsible = FALSE,
                inFluidRow = FALSE
            )
        ),
        fluidRow(
            box(
                title = "Score Type Switchboard",
                status = "primary",
                solidHeader = TRUE,
                collapsible = FALSE,
                width = 12,
                radioButtons(
                    ns("scoreType"),
                    "Score Type",
                    choiceNames  = c(
                        # "Fraction GC",  
                        "GC Residual Z Score"
                        # "GC Residual Bias" 
                    ),
                    choiceValues = c(
                        # "gc", 
                        "gcrz"
                        # "gcrz_bias"
                    ),
                    selected = "gcrz",
                    inline = TRUE,
                    width = "100%"
                )
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("sampleDistributionPlot"), 
                "Sample Distributions",
                width = 6,
                status = "primary",
                collapsible = TRUE,
                solidHeader = TRUE,
                collapsed = FALSE
            ),
            staticPlotBoxUI(
                ns("sampleDeltaDistributionPlot"), 
                "Gene Target Delta Distributions",
                width = 6,
                status = "primary",
                collapsible = TRUE,
                solidHeader = TRUE,
                collapsed = FALSE
            )
        ),
        NULL
    )
}
