#----------------------------------------------------------------------
# erfsCalls trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_erfsCallsTrack <- function(...) {
    list( 
        click = FALSE,
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        expand = FALSE,
        navigation = TRUE, # whether the track offers a custom, additional row of within-track navigation inputs
        NULL
    )
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.erfsCallsTrack <- showTrackSourcesDialog

# build method for the S3 class; REQUIRED
build.erfsCallsTrack <- function(track, reference, coord, layout){
    req(coord, coord$chromosome)

    startSpinner(session, message = "getting bins")
    sourceId <- track$settings$items()[[1]]$Source_ID
    req(sourceId)
    cd <- erfsCallData(sourceId)
    req(cd)

    # use generic methods and any other custom code to determine the track's (dynamic) Y span
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    ylim <- c(0, 2)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = "Calls", # yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"
        y <- 1
        for(gene_target in names(cd$calls)){
            df <- cd$calls[[gene_target]]
            I <- df[, chrom == coord$chromosome & start0 < coord$end & end1 >= coord$start]
            if(any(I)) rect(
                df[I, start0 + 1], y, 
                df[I, end1],       y + 1, 
                col = traceColors$gene_target[[gene_target]], 
                border = NA,
            )
            y <- y - 1
        }
        trackLegend(
            track, coord, ylim, bty = "n", 
            legend = toupper(names(traceColors$gene_target)),
            col = unlist(traceColors$gene_target),
            pch = 16,
            cex = 1.1
        )
    })

    # return the track's magick image and associated metadata
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
# only one navigation set is shown per track, your navigation should decide how to handle multiple regions
navigation.erfsCallsTrack <- function(track, session, id, browser){
    navTable <- initTrackNav(track, session, "navTable")
    trackNavData <- reactive({
        sourceId <- track$settings$items()[[1]]$Source_ID
        req(sourceId)
        cd <- erfsCallData(sourceId)
        req(cd)
        rbind(
            cd$calls$brca2[, .(chrom, start0, end1, size = end1 - start0, target = "BRCA2")],
            cd$calls$rad51[, .(chrom, start0, end1, size = end1 - start0, target = "RAD51")]
        )[order(chrom, start0, end1)]
    })
    tagList(
        trackNavTable(
            track, 
            session, 
            browser$id,
            navTable, # the name as provided by initTrackNav
            tableData = trackNavData, # populate a table based on track settings, etc.
            actionFn = function(selectedRow){
                req(selectedRow)
                d <- trackNavData()[selectedRow]
                handleTrackNavTableClick(1, track, d$chrom, d$start0, d$end1)
            }
            # add other argument to pass to bufferedTableServer, but omit "selection" and "width"
        )
        # etc.
    )
}
