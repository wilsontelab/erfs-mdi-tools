#----------------------------------------------------------------------
# erfsReadDepth trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_erfsReadDepthTrack <- function(...) {
    list( 
        click = FALSE,
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        expand = FALSE,
        NULL
    )
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.erfsReadDepthTrack <- showTrackSourcesDialog

# build method for the S3 class; REQUIRED
build.erfsReadDepthTrack <- function(track, reference, coord, layout){
    req(coord, coord$chromosome)

    startSpinner(session, message = "getting bins")
    sourceId <- track$settings$items()[[1]]$Source_ID
    req(sourceId)
    bd <- erfsBinData(sourceId)
    req(bd)
    binI <- bd$binCounts[chrom == coord$chromosome & start0 < coord$end & end1 >= coord$start, binI]
    req(any(binI))
    # gcBiasModels <- app$normalizeGC$getGcBiasModels_externalCall(sourceId)
    samples <- track$settings$get("ERFS_ReadDepth","Samples")

    # use generic methods and any other custom code to determine the track's (dynamic) Y span
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    y_max <- trimws(track$settings$get("Track","Y_Limit"))
    y_max <- if(y_max == "auto" || y_max == "") 100 else as.integer(y_max)
    muAdjusted <- track$settings$get("ERFS_ReadDepth","Mu_Adjusted")
    ylim <- if(muAdjusted) c(-10, y_max * 1.1) else c(0, y_max * 1.1)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  
            ylab = if(muAdjusted) "Read Count Delta" else "Read Count", 
            # yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"
        for(sample_name in samples){
            d <- getWindowBinData(sourceId, bd, sample_name, binI)
            x <- d[, start0 + (end1 - start0) / 2]
            y <- d$bin_count - (if(muAdjusted) d$mu else 0)
            points(
                x, 
                pmin(y, y_max), 
                pch = 19,
                cex = 1.25,
                col = ifelse(
                    d$copy_number >= 1 & d$copy_number <= 4,
                    traceColors$sample[[sample_name]], # addAlphaToColors(, alpha = 0.5)  
                    "grey40"
                )
            )
            abline(h=0, col = "grey40")
        }
        I <- names(traceColors$sample) %in% samples
        trackLegend(
            track, coord, ylim, bty = "n", 
            legend = toupper(names(traceColors$sample[I])),
            col = unlist(traceColors$sample[I]),
            pch = 16,
            cex = 1.1,
            pt.cex = 1.5
        )
    })

    # return the track's magick image and associated metadata
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}
