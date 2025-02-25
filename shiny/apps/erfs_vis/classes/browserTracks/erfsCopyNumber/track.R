#----------------------------------------------------------------------
# erfsCopyNumber trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_erfsCopyNumberTrack <- function(...) {
    list( 
        click = FALSE,
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        expand = FALSE,
        navigation = FALSE, # whether the track offers a custom, additional row of within-track navigation inputs
        NULL
    )
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.erfsCopyNumberTrack <- showTrackSourcesDialog

# build method for the S3 class; REQUIRED
build.erfsCopyNumberTrack <- function(track, reference, coord, layout){
    req(coord, coord$chromosome)

    startSpinner(session, message = "getting bins")
    sourceId <- track$settings$items()[[1]]$Source_ID
    req(sourceId)
    bd <- erfsBinData(sourceId)
    req(bd)
    binI <- bd$binCounts[chrom == coord$chromosome & start0 < coord$end & end1 >= coord$start, binI]
    req(any(binI))
    geneTarget <- track$settings$get("ERFS_Copy_Number","Gene_Target")

    # use generic methods and any other custom code to determine the track's (dynamic) Y span
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    ylim <- c(-0.5, 6.5)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  
            ylab = paste(cellLines[[geneTarget]],"CN"), 
            # yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

            abline(h = 0:6, col = "grey", lwd = 0.5)

            sample_name <- paste(geneTarget, "_ctl", sep = "")
            d <- getWindowBinData(sourceId, bd, sample_name, binI)
            points(
                d$start0 + (d$end1 - d$start0) / 2, 
                d$copy_number, 
                pch = 19,
                cex = 0.75
            )
    })

    # return the track's magick image and associated metadata
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}
