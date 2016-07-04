##########################################################################
setClass(
    Class = "ASpliFeatures",
    representation = representation(
      genes = "GRangesList",
      bins = "GRanges",
      junctions = "GRanges"))
##########################################################################
setClass(
  Class = "ASpliCounts",
  representation = representation(
    gene.counts = "data.frame", 
    exon.intron.counts = "data.frame",
    junction.counts = "data.frame",
    e1i.counts = "data.frame", 
    ie2.counts = "data.frame",
    gene.rd = "data.frame",
    bin.rd = "data.frame"))
##########################################################################
setClass(
  Class="ASpliAS",
  representation = representation(
    irPIR = "data.frame",
    altPSI = "data.frame",
    esPSI = "data.frame",
    junctionsPIR = "data.frame",
    junctionsPSI = "data.frame",
    join = "data.frame")
  )
##########################################################################
  setClass(
    Class = "ASpliDU",
    representation = representation(
      genes = "data.frame",
      bins = "data.frame",
      junctions = "data.frame"))
##########################################################################  
setGeneric (
    name = "binGenome",
    def = function(genome, md = NULL)
      standardGeneric("binGenome"))
##########################################################################
setMethod(
  f = "binGenome",
  signature = "TxDb",
  definition = function (genome,md = NULL){
    features <- new(Class = "ASpliFeatures")
    if (is.null(md))
    {
      md <- data.frame(names(transcriptsBy(genome)), stringsAsFactors = FALSE)
      row.names(md) <- names(transcriptsBy(genome))
      colnames(md) <- "symbol" 
    }
    sink("ASpli_binFeatures.log")
    genes.by.exons <- .createGRangesGenes(genome, md) 
    lg <- length(genes.by.exons)
    message("* Number of extracted Genes =",lg,"\n")
    cat("* Number of extracted Genes =")
    cat(lg)
    cat("\n")
##########################################################################
    exon.bins <- .createGRangesExons(genome, md)
    #add locus_overlap
    index <- match(exon.bins@elementMetadata$locus, names(genes.by.exons))
    locus_overlap <- rep("-", length(exon.bins))
    locus_overlap <- genes.by.exons@elementMetadata$locus_overlap[index]
    mcols(exon.bins) <- append(mcols(exon.bins), DataFrame(locus_overlap=locus_overlap))
    le <- length(exon.bins)
    message("* Number of extracted Exon Bins =",le)
    cat("* Number of extracted Exon Bins =")
    cat(le)
    cat("\n")
##########################################################################
    intron.tot <- .createGRangesIntrons(genome, md)
    #add locus_overlap
    index <- match(intron.tot@elementMetadata$locus, names(genes.by.exons))
    locus_overlap <- rep("-", length(intron.tot))
    locus_overlap <- genes.by.exons@elementMetadata$locus_overlap[index]
    mcols(intron.tot) <- append(mcols(intron.tot), 
                                DataFrame(locus_overlap=locus_overlap))
    li <- length(intron.tot)
    message("* Number of extracted intron bins =",li)
    cat("* Number of extracted intron bins =")
    cat(li)
    cat("\n")
##########################################################################
    transcripts <- .createGRangesTranscripts(genome)
    lt <- length(unlist(transcripts))
    message("* Number of extracted trascripts =",lt)
    cat("* Number of extracted trascripts =")
    cat(lt)
    cat("\n")
##########################################################################
    junctions <- .createGRangesJunctions(genome) 
    #add locus_overlap
    index <- match(junctions@elementMetadata$locus, names(genes.by.exons))
    locus_overlap <- rep("-", length(junctions))
    locus_overlap <- genes.by.exons@elementMetadata$locus_overlap[index]
    mcols(junctions) <- append(mcols(junctions), 
                               DataFrame(locus_overlap=locus_overlap))
    lj<-length(junctions)
    message("* Number of extracted junctions =",lj)
    cat("* Number of extracted junctions =")
    cat(lj)
    cat("\n")
##########################################################################
    intron.bins <- intron.tot[intron.tot@elementMetadata$feature == "I"]
    intron.orig <- intron.tot[intron.tot@elementMetadata$feature == "Io"]
    class <- rep("fullI", length(intron.orig))
    mcols(intron.orig) <- append(mcols(intron.orig), DataFrame(class=class))
    event <- rep("-", length(intron.orig))
    eventJ <- rep("-", length(intron.orig))
    mcols(intron.orig) <- append(mcols(intron.orig), DataFrame(event=event))
    mcols(intron.orig) <- append(mcols(intron.orig), DataFrame(eventJ=eventJ))
    ########aca la modificacion:
    exons.introns <- .findASBin(exon.bins, intron.bins, transcripts, junctions)
    fullT <- c(exons.introns,intron.orig)
    #add locus_overlap #queda agregado por la combinacion
    ############aca habria que armar la lista.
    features@genes <- genes.by.exons
    features@bins <- fullT
    features@junctions <- junctions
    sink()
    return(features)})
##########################################################################
setGeneric (
  name = "rds",
  def = function(counts, targets)
  standardGeneric("rds"))
##########################################################################
setMethod(
  f= "rds",
  signature = "ASpliCounts",
  definition = function(counts, targets)
  {
    geneStart <- ncol(countsg(counts))-nrow(targets)+1
    gene.rd <- cbind( countsg(counts)[,1:geneStart-1], 
                      countsg(counts)[,geneStart:ncol(countsg(counts))]/countsg(counts)$effective_length
                   )
    binStart <- ncol(countsb(counts))-nrow(targets)+1
    bin.rd <- cbind(countsb(counts)[, 1:binStart-1], 
                  countsb(counts)[,binStart:ncol(countsb(counts))]
                  /countsb(counts)$length)
                  
    tb <- match(bin.rd$locus, rownames(gene.rd))
          rdfinalb=cbind(bin.rd, bin.rd[,binStart:ncol(bin.rd)]
                   /gene.rd[tb,geneStart:ncol(countsg(counts))])
                            
    counts@gene.rd <- gene.rd
    counts@bin.rd <- rdfinalb
    return(counts)
  }
)
##############################################################################
setGeneric (
  name = "readCounts",
  def = function(features, bam, cores=NULL, l, maxISize, minAnchor=NULL)
  standardGeneric("readCounts"))
##########################################################################
setMethod(
  f = "readCounts",
  signature = "ASpliFeatures",
  definition = function(features, bam, cores=NULL, l, maxISize, minAnchor=NULL)
  {
      counts <- new(Class="ASpliCounts")
      if(is.null(minAnchor)){minAnchor=10}
      minA <- round(minAnchor*l/100)
      gene.hits <- .counterGenes(bam, featuresg(features), cores)
      message("read summarization by gene completed")
      counts@gene.counts <- gene.hits
      ####################################################################
      #exons
      bins <- featuresb(features)
      totCounts <- .counterBin(bam, bins, gene.hits, cores)
      message("read summarization by bin completed")
      ####################################################################
      #introns
      introns <- c(bins[bins@elementMetadata$feature == "I"], 
                 bins[bins@elementMetadata$feature == "Io" ],
                 bins[bins@elementMetadata$eventJ=="IR"])
      ####################################################################
      e1i <- introns
      start(e1i) <- start(introns)-(l-minAnchor)
      end(e1i) <- start(introns)+(l-minAnchor)
      e1i.hits <- .counterJbin(bam, e1i, gene.hits, cores, l)
      message("read summarization by ei1 region completed")
      #####################################################################
      ie2<-introns
      start(ie2) <- end(introns)-(l-minA)
      end(ie2) <- end(introns)+(l-minA)
      ie2.hits <- .counterJbin(bam, ie2, gene.hits, cores, l)
      message("read summarization by ie2 region completed")
      ####################################################################
       #junctions
      junction.hits = .counterJunctions(features, bam, cores, maxISize)
      message("junction summarization completed")
      #####################################################################
      counts@gene.counts <- gene.hits
      counts@exon.intron.counts <- totCounts
      counts@junction.counts <- junction.hits 
      counts@e1i.counts <- e1i.hits  
      counts@ie2.counts <- ie2.hits
      counts <- rds(counts, targets)
      return(counts)
})
##########################################################################
##########################################################################
setGeneric (
  name= "AsDiscover",
  def = function(counts, 
               targets, 
               features, 
               bam, 
               l, 
               pair,
               threshold=NULL, 
               cores=NULL)
  standardGeneric("AsDiscover")
  )
##########################################################################
setMethod(
    f = "AsDiscover",
  signature = "ASpliCounts",
  definition =function(counts, 
               targets, 
               features, 
               bam, 
               l, 
               pair,
               threshold=NULL, 
               cores=NULL)
    
  {
    as <- new(Class = "ASpliAS")
    if(is.null(threshold)){threshold=5}
    ### Common Part ############################################
    df0 <- countsj(counts)[countsj(counts)$multipleHit=="-",]
    df0 <- df0[df0$gene!="noHit",]
    jcounts <- .filterJunctionBySample(df0=df0, 
                                    targets=targets, 
                                    threshold=threshold)  #mean > one of the condition
##########################################################################
    #Junctions PSI:
    junctionsPSI <- .junctionsPSI_SUM(df0, targets, pair)
    as@junctionsPSI <- junctionsPSI
    message("Junctions PSI completed")
    ############es el mismo dataset#########################
    #JUNTURAS#######################
    #We analyze each junction as a putative intron. 
    #So, we estimate PIR using e1i, ie2 counts as an intron retention
    junctionsPIR <- .junctionsDiscover(df=jcounts, 
                                     bam, 
                                     cores=NULL, 
                                     l, 
                                     targets, 
                                     features, 
                                     pair)
    message("Junctions PIR completed")
    as@junctionsPIR <- junctionsPIR
    ############es el mismo dataset#########################
    split <- data.frame(matrix(unlist(strsplit(rownames(jcounts),  "[.]") ),
                            byrow=TRUE, ncol=3),   
                            row.names=rownames(jcounts),
                            stringsAsFactors=FALSE )  
    colnames(split) <- c("chrom","start","end")
    split$start <- as.numeric(split$start)
    split$end <- as.numeric(split$end)
    split$chrom <- as.character(split$chrom)
    jranges <- GRanges( 
    seqnames <- split$chrom,
    ranges  <- IRanges(start= split$start, 
                      end = split$end, 
                      names = rownames(split)))
    ########################################################
    #parte 1: irPIR intronJPIR_SUM_COND
    ic <- rbind(countsb(counts)[countsb(counts)$feature == "I",], 
              countsb(counts)[countsb(counts)$feature == "Io",], 
              countsb(counts)[countsb(counts)$event == "IR*",],
              countsb(counts)[countsb(counts)$event == "IR",])
    intranges <- featuresb(features)[rownames(ic)]
##########################################################################
    dfe1e2 <- .e1e2JPIR(intranges, jcounts)
    dfe1e2$J3 <- as.character(dfe1e2$J3)
    #junction spanning intron, from e1 to e2. Deafault junction
    #in functions...
##########################################################################
    countsSt <- ncol(counts@e1i.counts)-nrow(targets)+1
    e1i <- counts@e1i.counts[,countsSt:ncol(counts@e1i.counts)]
    ie2 <- counts@ie2.counts[,countsSt:ncol(counts@ie2.counts)] 
##########################################################################
    J1 <- paste(rownames(e1i),"E1I", sep="_")
    J2 <- paste(rownames(ie2),"IE2", sep="_")
    dffinal <- data.frame(J1, e1i, J2, ie2)
    #no tiene junturas, agregamos las columnas J1, J2
    empty <- data.frame(matrix(NA, 
                            nrow =  nrow(dffinal), 
                            ncol =  ncol(dfe1e2)), 
                            stringsAsFactors=FALSE)  
    colnames(empty) <- colnames(dfe1e2)
    int <- match(dfe1e2$jbin, row.names(dffinal))
    empty[int,] <- dfe1e2
    dmerge <- cbind(dffinal, empty); 
    dmerge$jbin <- NULL
    rownames(dmerge) <- names(intranges)
##########################################################################    
    intPIR <- function(x){
            x[is.na(x)] <- 0 # we have to remove NAs
            res <- (x[1]+x[2])/(x[1]+x[2]+2*x[3])
            return(res)
    }
##########################################################################
    df <- dmerge
    df$J1 <- NULL
    df$J2 <- NULL
    df$J3 <- NULL
    ff <- rep(targets$condition,3)
    myfactor <- factor(paste(ff, rep(1:3,each=length(targets$condition)),sep="."), 
                     levels=paste(rep(pair, each=3 ), 1:3, sep="."))
    colnames(df) <- paste(ff, rep(1:3,each=length(targets$condition)), sep=".")
    dfSum <- t(apply(df, 1, function(x){tapply(as.numeric(x), 
                                             INDEX=myfactor, sum  )}))
    colnames(dfSum) <- rep(unique(targets$condition), each=3) 
    myfactor <- factor(rep(unique(targets$condition),each=3), levels=pair)
    intPIRres <- cbind(event=ic$event,
                     dmerge, 
                     t(apply(dfSum, 1, 
                             function(x){tapply(as.numeric(x), 
                             INDEX=myfactor, 
                              intPIR  )})))
    message("Junctions IR PIR completed")
    as@irPIR <- intPIRres
##########################################################################
    #parte 2:altPSI, esPSI exonJPSI_SUM_COND
    ec <- countsb(counts)[countsb(counts)$feature == "E",]
    ec <- ec[ec$event != "IR",]
    ec <- ec[ec$event != "IR*",]
    exranges <- featuresb(features)[rownames(ec)]
    dfstart <- .startJPSI(jranges, exranges, jcounts, targets) 
    #identifies start junctio
    dfend <- .endJPSI(jranges, exranges, jcounts, targets)  #identifies end junction 
    dfwithin <- .withinJPSI(jranges, exranges, jcounts, targets)#identifies within junction
    #add the information at the end of the DF of exonsde
    corr.start <- match(rownames(dfstart), names(exranges))
    corr.end <- match(rownames(dfend), names(exranges))
    corr.w <- match(rownames(dfwithin), names(exranges))
##########################################################################
    dffinal <- data.frame(matrix(NA, 
                              nrow =  length(exranges), 
                              ncol = 3*ncol(dfstart)), 
                              stringsAsFactors=FALSE)
    colnames(dffinal) <- c(colnames(dfstart), colnames(dfend), colnames(dfwithin))
    colN <- rep(c("J", rownames(targets)), 3)
    colnames(dffinal) <- colN
    rownames(dffinal) <- names(exranges)
    startStart <- 1
    startEnd <- ncol(dfstart)
    endStart <- startEnd+ 1
    endEnd <- endStart + ncol(dfend) - 1
    withinStart <- endEnd + 1
    withinEnd <-  withinStart +ncol(dfwithin) - 1
##########################################################################
    dffinal[corr.start, startStart:startEnd] <- dfstart
    dffinal[corr.end, endStart:endEnd] <- dfend
    dffinal[corr.w, withinStart:withinEnd] <- dfwithin
    cN <- colnames(dffinal)
    cN <- cN[cN!="J"]
    dfEvents <- cbind(event=as.data.frame(exranges@elementMetadata$event), dffinal)
    colnames(dfEvents)[1] <- "event"
##########################################################################
    dfAlt <- rbind(dfEvents[dfEvents$event=="Alt3ss",], 
                 dfEvents[dfEvents$event=="Alt5ss",], 
                 dfEvents[dfEvents$event=="Alt3ss*",], 
                 dfEvents[dfEvents$event=="Alt5ss*",])
    dfAltSub <- dfAlt[,colnames(dfAlt)!="J"]
    dfAltSub <- dfAltSub[,colnames(dfAltSub)!="event"]
    colnames(dfAltSub)  <- cN
##########################################################################
    altPSI <- function(x){
      x[is.na(x)] <- 0 # we have to remove NAs
      res <- (x[1]+x[2])/(x[1]+x[2]+x[3])
      return(res)
    }
##########################################################################
    ff <- rep(targets$condition,3)
    colnames(dfAltSub) <- paste(ff, rep(1:3,each=length(targets$condition)))
    myfactor <- factor(paste(ff, rep(1:3,each=length(targets$condition))), 
                     levels=paste(rep(pair, each=3 ), 1:3))
    
    dfAltSum <- t(apply(dfAltSub, 1, 
                function(x){tapply(as.numeric(x), 
                                         INDEX=myfactor,
                                         sum  )}))
    
    colnames(dfAltSum) <- rep(unique(targets$condition), each= 3 )
    myfactor <- factor(rep(unique(targets$condition),each= 3 ), levels=pair)
    altPSIRes <- cbind(dfAlt, t(apply(dfAltSum, 1, 
                                    function(x){tapply(as.numeric(x), 
                                                       INDEX=myfactor, 
                                                       altPSI  )})))
    colnames(altPSIRes) <- colnames(intPIRres)
    message("Junctions AltSS PSI completed")
    as@altPSI <- altPSIRes
##########################################################################
    dfES <- rbind(dfEvents[dfEvents$event=="ES",], dfEvents[dfEvents$event=="-",], dfEvents[dfEvents$event=="ES*",] )
    dfESSub <- dfES[,colnames(dfES)!="J"]
    dfESSub <- dfESSub[,colnames(dfESSub)!="event"]
    colnames(dfESSub)  <-cN
##########################################################################
    EsPSI <- function(x){
      x[is.na(x)] <- 0 # we have to remove NAs
      res <- (x[1]+x[2])/(x[1]+x[2]+x[3]*2)
      return(res)
    }
##########################################################################
    colnames(dfESSub) <- paste(ff, rep(1:3,each=length(targets$condition)))
    myfactor <- factor(paste(ff, rep(1:3,each=length(targets$condition))), 
                 levels=paste(rep(pair, each=3 ), 1:3))
    dfESSum <- t(apply(dfESSub, 1, 
               function(x){tapply(as.numeric(x), 
                                        INDEX=myfactor,
                                        sum  )}))
    colnames(dfESSum) <- rep(unique(targets$condition), each=3)
    myfactor <- factor(rep(unique(targets$condition), each=3), levels=pair)
    EsPSIRes <- cbind(dfES, t(apply(dfESSum, 1, 
                                  function(x){tapply(as.numeric(x),
                                                     INDEX=myfactor, 
                                                     EsPSI  )})))
##########################################################################
    colnames(EsPSIRes) <- colnames(intPIRres)
    message("Junctions ES PSI completed")
    as@esPSI <- EsPSIRes 
    #junctions discovery
    as@join <- rbind(altPSIRes,EsPSIRes, intPIRres)
    return(as)
})
##########################################################################
setGeneric (
  name = "writeAS",
  def = function(as, output.dir="as")
  standardGeneric("writeAS"))
##########################################################################
setMethod(
  f = "writeAS",
  signature = "ASpliAS",
  definition =function(as, output.dir="as")
  {
    currentDir <- getwd()
      outputDir <- paste(currentDir, output.dir, sep = "/")       
      exonsDir <- paste(outputDir, "exons", sep = "/")       
      intronsDir <- paste(outputDir, "introns", sep = "/")
      junctionsDir <- paste(outputDir, "junctions", sep = "/")                     	     	     	     	    
    
    if (!file.exists(output.dir))
    {
      dir.create(outputDir)
      dir.create(exonsDir)
      dir.create(intronsDir)
      dir.create(junctionsDir)
    }
    ################ EXONS ####################################
    file <- paste(exonsDir, "exon.altPSI.tab", sep="/")
    write.table(altPSI(as), file, sep="\t", quote=FALSE, col.names=NA)
    file <- paste(exonsDir, "exon.altES.tab", sep="/")
    write.table(esPSI(as), file, sep="\t", quote=FALSE, col.names=NA)
    ################ INTRONS ##################################
    file <- paste(intronsDir, "intron.irPIR.tab", sep="/")
    write.table(irPIR(as), file, sep="\t", quote=FALSE, col.names=NA)
    ################ Junctions ################################
    file <- paste(junctionsDir, "junction.PIR.tab", sep="/")
    write.table(as@junctionsPIR, file, sep="\t", quote=FALSE, col.names=NA)
    file <- paste(junctionsDir, "junction.PSI.tab", sep="/")
    write.table(junctionsPSI(as), file, sep="\t", quote=FALSE, col.names=NA)
    ################ ALL ######################################
    currentDir=getwd()
    file <- paste(outputDir, "as_discovery.tab", sep="/")
    write.table(joint(as), file, sep="\t", quote=FALSE, col.names=NA)
    }
)
##########################################################################
setGeneric (
  name = "DUreport",
  def = function(counts, targets, pair, group,
               minGenReads=NULL,
               minBinReads=NULL,
               minRds=NULL,
               ignoreExternal=NULL,
               threshold=NULL)
  standardGeneric("DUreport"))
##########################################################################
setGeneric (
  name = "DUreport_DEXSeq",
  def = function(counts, targets, pair, group,
               minGenReads=NULL,
               minBinReads=NULL,
               minRds=NULL,
               threshold=NULL)
  standardGeneric("DUreport_DEXSeq"))
##########################################################################
setMethod(
  f = "DUreport",
  signature = "ASpliCounts",
  definition =function(counts, targets, pair, group, 
                      minGenReads=NULL,
                      minBinReads=NULL,
                      minRds=NULL,
                      ignoreExternal=NULL,
                      threshold=NULL)
  {
    du <- new(Class="ASpliDU")
    #define parameters#
    if(is.null(minGenReads)){minGenReads=10}
    if(is.null(minBinReads)){minBinReads=5}
    if(is.null(minRds)){minRds=0.05}
    if(is.null(threshold)){threshold=5}
    ##################################
    df0 <- countsg(counts)
    dfG0 <- .filterByReads(df0=df0,
                    targets=targets,
                    min=minGenReads,
                    type="any")
    dfGen <- .filterByRdGen(df0=dfG0,
                     targets=targets,
                     min=minRds,
                     type="any")
    genesde <- .genesDE(df=dfGen, 
                 targets=targets, 
                 pair=pair,
                 group=group) 
    message("Genes DE completed")
    du@genes <- genesde
    ##########################################################################
    dfG0 <- .filterByReads(df0=df0,
                    targets=targets,
                    min=minGenReads,
                    type="all")
    dfGen <- .filterByRdGen(df0=dfG0,
                     targets=targets,
                     min=minRds,
                     type="all")
    dfBin <- countsb(counts)[countsb(counts)[,"locus"]%in%row.names(dfGen),]
    df1 <- .filterByReads(df0=dfBin,
                   targets=targets, 
                   min=minBinReads,
                   type="any")
    df2 <- .filterByRdBinRATIO(
    dfBin=df1,
    dfGen=dfGen,
    targets=targets, 
    min=minRds,
    type="any")
#bins con AS en binsN
  if (is.null(ignoreExternal)){ignoreExternal=TRUE}
  binsdu <- .binsDU(df=df2,
                  targets,
                  dfGen,
                  ignoreExternal=ignoreExternal, 
                  pair,
                  group) 
du@bins <- binsdu
message("Bins DU completed")
##########################################################################
df0 <- countsj(counts)[countsj(counts)[,"gene"]%in%rownames(dfGen),]
df <- .filterJunctionBySample(df0=df0, 
                             targets=targets, 
                             threshold=threshold)  #mean > one of the condition
df <- df[df$multipleHit=="-",]
junctionsdeSUM <- .junctionsDU_SUM(df,
                                    targets, 
                                    genesde,
                                    pair, 
                                    group,dfGen)
du@junctions <- junctionsdeSUM
message("Junctions DU completed")
return(du)
}
)
#########################################################################
setMethod(
  f = "DUreport_DEXSeq",
  signature = "ASpliCounts",
  definition = function(counts, targets, pair, group, 
                      minGenReads=NULL,
                      minBinReads=NULL,
                      minRds=NULL,
                      threshold=NULL)
  {
    du <- new(Class="ASpliDU")
    #define parameters#
    if(is.null(minGenReads)){minGenReads=10}
    if(is.null(minBinReads)){minBinReads=5}
    if(is.null(minRds)){minRds=0.05}
    if(is.null(threshold)){threshold=5}
    ###############################################
    df0 <- countsg(counts)
    dfG0 <- .filterByReads(df0=df0,
                         targets=targets,
                         min=minGenReads,
                         type="any")
    dfGen <- .filterByRdGen(df0=dfG0,
                          targets=targets,
                          min=minRds,
                          type="any")
    genesde <- .genesDE_DESeq(df=dfGen, 
                            targets=targets, 
                            pair=pair) 
    du@genes <- genesde #idem ???
    ###############################################################
    dfG0 <- .filterByReads(df0=df0,
                         targets=targets,
                         min=minGenReads,
                         type="all")
    dfGen <- .filterByRdGen(df0=dfG0,
                          targets=targets,
                          min=minRds,
                          type="all")
    
    dfBin <- countsb(counts)[countsb(counts)[,"locus"]%in%row.names(dfGen),]
    df1 <- .filterByReads(df0=dfBin,
                        targets=targets, 
                        min=minBinReads,
                        type="any")
    df2 <- .filterByRdBinRATIO(
      dfBin=df1,
      dfGen=dfGen,
      targets=targets, 
      min=minRds,
      type="any")
    #bins con AS en binsN
    binsdu <- .binsDU_DEXSeq(df=df2,
                           targets=targets,
                           group=group) 
    du@bins <- binsdu
    ########################################################################
    df0 <- countsj(counts)[countsj(counts)[,"gene"]%in%rownames(dfGen),]
    df <- .filterJunctionBySample(df0=df0, 
                                targets=targets, 
                                threshold=threshold)  #mean > one of the condition
    df <- df[df$multipleHit=="-",]
    junctionsdeSUM <- .junctionsDU_SUM_DEXSeq(df,
                                            targets=targets, 
                                            genesde=genesde,
                                            group=group)
    du@junctions <- junctionsdeSUM
    return(du)
  }
)
##########################################################################
setGeneric (
  name = "writeDU",
  def =function(du, output.dir="du")
  standardGeneric("writeDU"))
##########################################################################
setMethod(
  f = "writeDU",
  signature = "ASpliDU",
  definition =function(du, output.dir="du")
  {
    currentDir <- getwd()
    outputDir <- paste(currentDir, output.dir, sep = "/")       
    genesDir <- paste(outputDir, "genes", sep = "/")       
    exonsDir <- paste(outputDir, "exons", sep = "/")       
    intronsDir <- paste(outputDir, "introns", sep = "/")
    junctionsDir <- paste(outputDir, "junctions", sep = "/")                                   	     	    
    if (!file.exists(output.dir))
    {
      dir.create(outputDir)
      dir.create(genesDir)
      dir.create(exonsDir)
      dir.create(intronsDir)
      dir.create(junctionsDir)
      }
      ################ GENES ####################################
    file <- paste(genesDir, "gene.de.tab", sep="/")
    write.table(du@genes, file, sep="\t", quote=FALSE, col.names=NA)
    ################ EXONS ###################################
    exonBins <- binsDU(du)[binsDU(du)$feature == "E",]
    exonBins <- exonBins[exonBins$event !="IR",]
    file <- paste(exonsDir, "exon.du.tab", sep="/")
    write.table(exonBins, file, sep="\t", quote=FALSE, col.names=NA)
    ################ INTRONS ##################################
    intronBins<-rbind(binsDU(du)[binsDU(du)$feature == "I",], 
                      binsDU(du)[binsDU(du)$feature == "Io",],
                      binsDU(du)[binsDU(du)$event == "IR",])
    file <- paste(intronsDir, "intron.du.tab", sep="/")
    write.table(intronBins, file, sep="\t", quote=FALSE, col.names=NA)
    ################ Junctions ##################################
    file <- paste(junctionsDir, "junction.du.tab", sep="/")
    write.table(junctionsDU(du), file, sep="\t", quote=FALSE, col.names=NA)
  }
)
##########################################################################
