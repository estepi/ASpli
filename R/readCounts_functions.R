#funciones para read counts
.counterGenes <-
  function(reads, feature, cores=NULL)
  {
    if (is.null(cores) )
      #it wont paralelize, it wouldn require multicore library
    {
      hits <- lapply(reads, function(x){countOverlaps(feature, x,  
                                                    ignore.strand = TRUE) })
    }
    else
    {
      hits <- mclapply(reads, mc.cores=cores, function(x){countOverlaps(feature, x,   
                                                                      ignore.strand = TRUE)})  
    }
    hits.ul <- do.call(cbind.data.frame, hits)
    wi <- width(feature)
    wit <- sum(wi)
    geneStarts <- sapply(start(feature),min)
    geneEnds <- sapply(end(feature),max)
    geneChr <- sapply(seqnames(feature),unique)
    strand <- min(strand(feature))
    strand[strand==1] <- "+"
    strand[strand==2] <- "-"
    genes <- GRanges(seqnames=geneChr,
                   strand=strand,
                   ranges=IRanges(geneStarts,geneEnds), 
                   effective_length=wit) 
    names(genes) <- names(feature)
    gene_coordinates <- feature@elementMetadata$gene_coordinates
    mcols(genes) <- append(mcols(genes), DataFrame(gene_coordinates=gene_coordinates))
    locus_overlap <- feature@elementMetadata$locus_overlap
    mcols(genes) <- append(mcols(genes), DataFrame(locus_overlap=locus_overlap))
    symbol <- feature@elementMetadata$symbol #symbol
    mcols(genes) <- append(mcols(genes), DataFrame(symbol=symbol))  
    aa <- data.frame(as.data.frame(genes@elementMetadata$symbol), 
                   as.data.frame(genes@elementMetadata$locus_overlap), 
                   as.data.frame(genes@elementMetadata$gene_coordinates), 
                   as.data.frame(genes@ranges), 
                   effective_length=as.data.frame(genes@elementMetadata$effective_length), 
                   hits)
    colnames(aa)[1:8] <- c("symbol","locus_overlap","gene_coordinates",
                          "start","end", "length","name","effective_length" ) 
    aa$name <- NULL
    return(aa)
  }
#########################################################
.counterBin <-
  function(reads, feature, genes, cores=NULL) 
  { 
    if (is.null(cores) )#no aplica paralelizacion, no necesita library de paralelizacion
    {
      hits <- lapply(reads,function(x){countOverlaps(feature, x,  
                                                   ignore.strand = TRUE)})#OK
    }
    else
    {
      hits <- mclapply(reads,mc.cores=cores,function(x){countOverlaps(feature, x,  
                                                                    ignore.strand = TRUE)})#OK  
    }
    hits.ul <- do.call(cbind.data.frame, hits)
    te <- match(feature@elementMetadata$locus, rownames(genes))
    gene_coordinates <- genes$gene_coordinates[te]
    aa <- data.frame(as.data.frame(feature@elementMetadata$feature),
                   as.data.frame(feature@elementMetadata$event),   
                   as.data.frame(feature@elementMetadata$locus), 
                   as.data.frame(feature@elementMetadata$locus_overlap), 
                   as.data.frame(feature@elementMetadata$symbol), 
                   gene_coordinates, 
                   as.data.frame(feature@ranges),  hits)
    aa$names <-  NULL
    colnames(aa)[1:9] <- c("feature","event","locus","locus_overlap","symbol",
                          "gene_coordinates","start","end","length")
    return(aa)
  }
#########################################################
.counterJbin <-
  function(reads, feature, genes, cores=NULL, l) 
  {
    ungapped <- lapply(reads, function(x) {x[njunc(x)==0,]}) #extraigo los que no tienen GAPS
    if (is.null(cores) )#no aplica paralelizacion, no necesita library de paralelizacion
    {
      hits <- lapply(ungapped,function(x){countOverlaps(feature, x, 
                                                      ignore.strand = TRUE, minoverlap = l)})#OK
    }
    else
    {
      hits <- mclapply(ungapped,mc.cores=cores,function(x){countOverlaps(feature, x,  
                                                                       ignore.strand = TRUE,  minoverlap = l)})#OK  
    }
    hits.ul <- do.call(cbind.data.frame, hits)
    te <- match(feature@elementMetadata$locus, rownames(genes))
    gene_coordinates <- genes$gene_coordinates[te]
    aa <- data.frame(as.data.frame(feature@elementMetadata$event),   
                   as.data.frame(feature@elementMetadata$locus), 
                   as.data.frame(feature@elementMetadata$locus_overlap), 
                   as.data.frame(feature@elementMetadata$symbol), 
                   gene_coordinates, 
                   as.data.frame(feature@ranges),  
                   hits)
    
    aa$names <- NULL
    colnames(aa)[1:8] <- c("event","locus","locus_overlap","symbol",
                          "gene_coordinates","start","end","length")
    return(aa)
  }
###############################################################
.ovBinJunction<-function(features, jranges)
{
  annJunctions <- featuresj(features)
  jname <- rep("-", length(jranges))
  hitBin <- rep("-", length(jranges))
  hitGen <- rep("-", length(jranges))
  hitGenStrand <- rep("*", length(jranges))
  gene_coordinates <- rep("-", length(jranges))
  ambiguos <- rep("-", length(jranges))
  j_within_bin <- rep("-", length(jranges))
  ##############nuevos#####################
  feature <- featuresg(features)
  geneStarts <- sapply(start(feature),min)
  geneEnds <- sapply(end(feature),max)
  geneChr<-sapply(seqnames(feature),unique)
  strand <- min(strand(feature))
  strand[strand==1] <- "+"
  strand[strand==2] <- "-"
  genes <- GRanges(seqnames=geneChr,
                 strand=strand,
                 ranges=IRanges(geneStarts,geneEnds), 
                 gene_coordinates=feature@elementMetadata$gene_coordinates,
                 symbol=feature@elementMetadata$symbol)
  names(genes) <- names(feature)
  overGene <- findOverlaps(jranges, genes, type="within")
  overGeneDF <- as.data.frame(overGene)#get a df
  posJrange <- overGeneDF$queryHits 
  #replace index numbers by names
  posGene <- overGeneDF$subjectHits
  #replace index numbers by names; posGene[1:10]
  overGeneDF$queryHits <- names(jranges)[as.numeric(overGeneDF$queryHits)]
  overGeneDF$subjectHits <- names(genes)[as.numeric(overGeneDF$subjectHits)]
  table <- table(overGeneDF$queryHits)
  ttG <- data.frame(aggregate(subjectHits ~ queryHits, data = overGeneDF, paste, collapse=";")) 
  dd0 <- match(ttG$queryHits,names(jranges))
  hitGen[dd0] <- ttG$subjectHits
  dd <- match(ttG$queryHits,names(table))
  ttG$undef <- table[dd]
  ttG$tag <- rep("-",nrow(ttG))
  ttG$tag[ttG$undef>1] <- "yes"
  ambiguos[dd0] <- ttG$tag
  ################################################
  hitGen[posJrange] <- names(genes[posGene])
  #remplazo usando el indice de la juntura, el nombre del gen
  hitGen[-posJrange] <- "noHit"
  hitGenStrand[posJrange] <- as.character(strand(genes)[posGene])
  gene_coordinates[posJrange] <- genes@elementMetadata$gene_coordinates[posGene] 
  #short coord viene del objeto genes
  #######################################################################
  overJ <- findOverlaps(jranges, annJunctions, type="equal") #identify annotated junctions
  overJDF <- as.data.frame(overJ) #get a df
  namesJ <- as.numeric(overJDF[,1])  #get index of jrangs that hit against annJunctions
  namesAnnJ <- as.numeric(overJDF[,2]) #get index of annJunctions thta hit against jranges
  jname[namesJ] <- names(jranges[namesJ])
  hitBin[namesJ] <- names(annJunctions[namesAnnJ]) #ok, metadata vector
  hitBin[-namesJ] <- "noHit" #ok,  metadata vector. 
  #name of the annotated junction  in the positon of the experimental one
  #Identify which gene contain the junction. 
  #In case I have not hit against annotated junction, this is a very useful
  #information
  ##########spanning exons bins
  #3 useful information about which bins the junctions span. 
  #any krahit is considered # aca esta el problema
  exonsBins <- featuresb(features)[featuresb(features)@elementMetadata$feature=="E",]
  over <- findOverlaps(jranges, exonsBins)
  overDF <- as.data.frame(over)
  namesJ <- as.numeric(overDF[,1])
  overDF[,1] <- names(jranges[namesJ])
  namesBins <- as.numeric(overDF[,2])
  overDF[,2] <- names(exonsBins[namesBins])
  tt <- data.frame(aggregate(subjectHits ~ queryHits, data = overDF, paste, collapse=";")) 
  span <- rep("-", length(jranges))
  te <- match(names(jranges), tt$queryHits) #ok
  span <- tt$subjectHits[te]
  #####################################################################
  overJunctionWithinBins <- findOverlaps(jranges, exonsBins, type="within")
  overJunctionWithinBinsDF <- as.data.frame(overJunctionWithinBins)
  namesJ <- as.numeric(overJunctionWithinBinsDF[,1])
  namesB <- as.numeric(overJunctionWithinBinsDF[,2])
  overJunctionWithinBinsDF[,1] <- names(jranges[namesJ])
  overJunctionWithinBinsDF[,2] <- names(exonsBins[namesB])
  agtt <- data.frame(aggregate(subjectHits ~ queryHits,
                            data = overJunctionWithinBinsDF, paste, collapse=";")) 
  tw <- match(names(jranges), agtt$queryHits) #ok;
  j_within_bin <- agtt$subjectHits[tw]
  symbol <- rep("-", length(jranges))    
  symbol[posJrange] <- as.character(genes@elementMetadata$symbol[posGene])
  mcols(jranges) <- append(mcols(jranges), DataFrame(hitBin=hitBin, 
                                                    hitGen=hitGen, 
                                                    hitGenStrand=hitGenStrand,
                                                    gene_coordinates=gene_coordinates,
                                                    undef=ambiguos,
                                                    bin_spanned=span,
                                                    j_within_bin=j_within_bin,
                                                    symbol=symbol))
  return(jranges)
}  
##########################################################
.counterJunctions <-
  function(features, bam, cores, maxISize)
  {
    if (is.null(cores) )#no aplica paralelizacion, no necesita library de paralelizacion
    {
      
      ujunctions <- lapply (bam, function(x)    {  
        junctions <- unlist(junctions(x) )
        strand(junctions) <- "*"
        start(junctions) <- start(junctions)-1
        end(junctions) <- end(junctions)+1
        ujunctions <- unique(junctions)
        return(ujunctions)   
      } )
    }
    else
    {
      ujunctions <- mclapply(bam, mc.cores=cores, function(x) {
        junctions <- unlist(junctions(x) )
        strand(junctions) <- "*"
        start(junctions) <- start(junctions)-1
        end(junctions) <- end(junctions)+1
        ujunctions <- unique(junctions)
        return(ujunctions)
        
      })
    }  
    #here I have unique elements of all the junctiosn
    jranges <- unique(unlist(GRangesList(unlist(ujunctions))))
    maxWidth <-  maxISize+2
    jranges <- jranges[width(jranges)<= maxISize]
    #Here I summarize hits agains the element
    fcoord <- paste(seqnames(jranges), 
                  start(jranges),
                  end(jranges) , sep="." )
    jranges@ranges@NAMES <- fcoord
    #########################################
    jcounts<-lapply(bam, function(x)
    {
      junctions <- unlist(junctions(x) )
      strand(junctions)<- "*"
      start(junctions) <- start(junctions)-1
      end(junctions) <- end(junctions)+1
      count <- countMatches(jranges, junctions)    
      jc <- data.frame(row.names=names(jranges), count)
      return(jc)
    })
    df <- do.call("cbind", jcounts); head(df)
    colnames(df) <- names(jcounts)
    #desde aca la bifurcacion:parte critica!
    jranges <- .ovBinJunction(features, jranges)
    
    jrdf <- data.frame(
      as.data.frame(jranges@elementMetadata$hitBin), 
      as.data.frame(jranges@elementMetadata$hitGen), 
      as.data.frame(jranges@elementMetadata$hitGenStrand),
      as.data.frame(jranges@elementMetadata$undef),
      as.data.frame(jranges@elementMetadata$symbol), 
      as.data.frame(jranges@elementMetadata$gene_coordinates),
      as.data.frame(jranges@elementMetadata$bin_spanned),
      as.data.frame(jranges@elementMetadata$j_within_bin),
      row.names=names(jranges)  )
    colnames(jrdf) <- c("junction", 
                     "gene", 
                     "strand",
                     "multipleHit",
                     "symbol",
                     "gene_coordinates",
                     "bin_spanned",
                     "j_within_bin")
    
    rownames(jrdf) <- names(jranges)
    aa <- merge(jrdf, df, by.x="row.names", by.y="row.names", sort=FALSE)
    rnames <- paste(start(jranges)-1,end(jranges)+1, sep="-" )
    rownames(aa) <- fcoord
    aa$Row.names <- NULL
    return(aa)
  }
################################################################
