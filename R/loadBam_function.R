################### loadBAM###################################
loadBAM <-
  function(targets, cores=NULL) 
  {
    #1.-read files and get a list
    files <- as.character(targets$bam)
    #2.-read alignments
    if (is.null(cores) )
    {
      datac <- lapply(
        
        files,
        function(x) 
        { 
          aln <- readGAlignments(x)
          return(aln) 
        }
      )
    }
    
    else
    {
      
      datac <- mclapply( 
        files, mc.cores = cores, 
        function(x) 
        { 
          aln <- readGAlignments(x)
          return(aln) 
        } 
      )
    }
    
    names(datac) = rownames(targets)
    return(datac)
    
  }
##############################################################
