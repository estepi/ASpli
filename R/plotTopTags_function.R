plotTopTags<-function(auxdf, 
                      genome, 
                      targetsPlot, 
                      output.dir)
{
  alns <- apply(targetsPlot,1,
              function(x)
              {
                AlignmentsTrack(
                  x[1], 
                  isPaired = FALSE, 
                  name = x[2],
                  col=x[3],
                  fill.coverage=x[3])
              }
  )
  auxdf$gene_coordinates <- as.character(auxdf$gene_coordinates)
  for (i in 1:nrow(auxdf))
  {
    title = paste(row.names(auxdf[i,]), auxdf$event[i])
    message("Plot: ", title)
    split1 <- matrix(unlist(strsplit(auxdf$gene_coordinates[i], "[:]") ), byrow=TRUE, ncol=2) 
    split2 <-  matrix(unlist(strsplit(split1[2],  "[-]") ), byrow=TRUE, ncol=2) 
    gen.start <- as.numeric(split2[1])
    gen.end <- as.numeric(split2[2])
    chrom <- split1[1]
    binStart <- auxdf$start[i] -20
    binEnd <- auxdf$end[i] +20
    ###########################################fix
    txTrAT <- GeneRegionTrack(genome, 
                              chromosome = chrom,
                              start = gen.start, 
                              end = gen.end,
                              options(ucscChromosomeNames=FALSE),
                              name = "Gene Model", 
                              transcriptAnnotation = "symbol",
                              size=3, 
                              fill="grey")
    ht1 <- HighlightTrack(
      trackList = c(alns, txTrAT), 
      start = binStart, 
      end = binEnd,
      chromosome = chrom,
      col="black",
      fill=NA,
      lwd=2,
      inBackground=FALSE
    )
    ###############################################################
    currentDir <- getwd()
    outputDir <- paste(currentDir, output.dir, sep = "/")       
    if(!file.exists(output.dir))
    {
      dir.create(outputDir)
    }
    file <- paste(row.names(auxdf[i,]), "png", sep=".")
    file <- sub(":", "_", file, perl=TRUE)      
    png(filename = paste(outputDir, file, sep = "/"))
    plotTracks(ht1, 
               from = gen.start, 
               to = gen.end,
               chromosome = chrom,
               min.height = 0, 
               minCoverageHeight = 0, 
               type="coverage",
               #type=c("pileup","coverage"),
               min.width = 3, 
               min.distance = 5,
               #ackground.title = "darkblue",
               main=title)
    dev.off()
    
  }
  message("End of plots")
}
