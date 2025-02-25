#----------------------------------------------------------------------
# UI components for the erfs_normalizeGC appStep module
#----------------------------------------------------------------------

# module ui function
erfs_normalizeGCUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$erfs_normalizeGC)

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
                "GC Bias",
                interactiveScatterplotUI(ns("gcBiasPlot"), height = '400px')
            ),
            staticPlotBoxUI(
                ns("binCountsPlot"), 
                title = "Bin Count Redux",
                width = 6,
                collapsible = TRUE,
                collapsed = FALSE,
                solidHeader = TRUE,
                status = "primary"
            )
        ),
        fluidRow(
            column(
                width = 6,
            ),
            column(
                width = 6,
                fluidRow(
                    column(
                        width = 6,
                        numericInput(
                            ns("mu1"),
                            label = "Mu CN 1",
                            value = 10, min = 0, max = 100, step = 1
                        )
                    ),
                    column(
                        width = 6,
                        numericInput(
                            ns("theta1"),
                            label = "Theta CN 1",
                            value = 50, min = 10, max = 1000, step = 10
                        )
                    ),
                ),
                fluidRow(
                    column(
                        width = 6,
                        numericInput(
                            ns("mu2"),
                            label = "Mu CN 2",
                            value = 20, min = 0, max = 100, step = 1
                        )
                    ),
                    column(
                        width = 6,
                        numericInput(
                            ns("theta2"),
                            label = "Theta CN 2",
                            value = 50, min = 10, max = 1000, step = 10
                        )
                    ),
                ),
                fluidRow(
                    column(
                        width = 6,
                        numericInput(
                            ns("mu3"),
                            label = "Mu CN 3",
                            value = 30, min = 0, max = 100, step = 1
                        )
                    ),
                    column(
                        width = 6,
                        numericInput(
                            ns("theta3"),
                            label = "Theta CN 3",
                            value = 50, min = 10, max = 1000, step = 10
                        )
                    ),
                ),  
                fluidRow(
                    column(
                        width = 6,
                        numericInput(
                            ns("mu4"),
                            label = "Mu CN 4",
                            value = 40, min = 0, max = 100, step = 1
                        )
                    ),
                    column(
                        width = 6,
                        numericInput(
                            ns("theta4"),
                            label = "Theta CN 4",
                            value = 50, min = 10, max = 1000, step = 10
                        )
                    ),
                ),
                fluidRow(
                    column(
                        width = 6,
                        numericInput(
                            ns("mu5"),
                            label = "Mu CN 5",
                            value = 50, min = 0, max = 100, step = 1
                        )
                    ),
                    column(
                        width = 6,
                        numericInput(
                            ns("theta5"),
                            label = "Theta CN 5",
                            value = 50, min = 10, max = 1000, step = 10
                        )
                    ),
                ),
            )
        ),
        NULL
    )
}
