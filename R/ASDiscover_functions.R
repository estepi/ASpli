#######################ASdiscovery FUNCTIONS###############################
.filterJunctionBySample <- function(df0, targets, threshold)
{
  
  a <- table(targets$condition)
  start <- ncol(df0)-sum(a)+1
  jm <- df0[,start:ncol(df0)]
  colnames(jm) <- targets$condition
  ffJ <- t(apply(jm, 1, 
          function(x){tapply(as.numeric(x), 
                                 INDEX=colnames(jm), 
                                  min )})>threshold)
  ii <- rowSums(ffJ)>0#at least in one condition
  df <- df0[ii,]
  return (df)
}
################################################################
.e1e2JPIR<-function(intranges, jcounts)
{
  split <- data.frame(matrix(unlist(strsplit(rownames(jcounts),  "[.]") ),
                   byrow=TRUE, ncol=3),   
                   row.names=rownames(jcounts),
                   stringsAsFactors=FALSE )  
  colnames(split) <- c("chrom","start","end")
  split$start <- as.numeric(split$start)  +1
  split$end <- as.numeric(split$end) -1
  split$chrom <- as.character(split$chrom)
  jranges <- GRanges( 
    seqnames = split$chrom,
                    ranges =IRanges(start= split$start, 
                    end = split$end, 
                    names = rownames(split)))
  jranges <- sort(jranges)
  #we have junctions coordenates equal to introns cordenates
  ############################################################
  jInt <- findOverlaps(jranges, intranges, type="equal") 
  #junctions equal introns bins
  jIntDF <- as.data.frame(jInt) 
  jIntDF$queryHits <- names(jranges)[jIntDF$queryHits]
  jIntDF$subjectHits <- names(intranges)[jIntDF$subjectHits]
  matchJ <- match(jIntDF$queryHits, rownames(jcounts))
  jc <- jcounts[matchJ, setdiff(colnames(jcounts),
                              c("junction",
                                "gene",
                                "strand",
                                "multipleHit",
                                "symbol",
                                "gene_coordinates", 
                                "bin_spanned","j_within_bin"))]
  #extraigo los conts de las junturas  que dan match igual
  ###########hago la estadistica ###############################
  dfe1e2 <- data.frame(J3=jIntDF$queryHits, 
                     jbin=jIntDF$subjectHits,
                     jc)
  return(dfe1e2)
}
###################################################
.startJPSI <- function(jranges, exranges, jcounts, targets)
{
  jranges.start <- jranges
  start(jranges.start) <- end(jranges.start)
  bin.start <- findOverlaps(exranges, jranges.start, type="start") 
  #query exon -- subject junctio
  bin.start.df <- as.data.frame(bin.start)
  bin.start.df$names <- names(exranges[bin.start.df$queryHits])
  bin.start.df$J1 <- names(jranges.start[bin.start.df$subjectHits])
  start <- ncol(jcounts) -nrow(targets) +1
  end <- ncol(jcounts)
  df1 <- data.frame(names=bin.start.df$names, 
                 jcounts[bin.start.df$subjectHits,start:end],
                 row.names=NULL) #recover counts
  tt1 <- data.frame(aggregate(. ~ names, data = df1, sum))
  tt2 <- data.frame(aggregate(J1 ~ names, data = bin.start.df, paste, collapse=";"))
  aa <- merge(tt2,tt1, by.x="names", by.y="names")#merge both dataframes
  rownames(aa) <- aa$names
  aa$names <-  NULL
  return(aa)
}
########################################################
.endJPSI <- function(jranges, exranges, jcounts, targets)
{
  jranges.end <- jranges
  end(jranges.end) <- start(jranges.end)
  bin.end <- findOverlaps(exranges, jranges.end, type="end") 
  bin.end.df <- as.data.frame(bin.end)
  bin.end.df$names <- names(exranges[bin.end.df$queryHits])
  bin.end.df$J2 <- names(jranges.end[bin.end.df$subjectHits])
  start <- ncol(jcounts) -nrow(targets) +1
  end <- ncol(jcounts)
  df1 <- data.frame(names=bin.end.df$names, 
                 jcounts[bin.end.df$subjectHits,start:end],
                 row.names=NULL) #recover counts
  tt1 <- data.frame(aggregate(. ~ names, data = df1, sum))   #aggregate counts
  tt2 <- data.frame(aggregate(J2 ~ names, data = bin.end.df, 
                           paste, collapse=";")) #aggregate by name
  aa <- merge(tt2,tt1, by.x="names", by.y="names")#merge both dataframes
  rownames(aa) <-  aa$names
  aa$names <-  NULL
  return(aa)
}
#######################################################
.withinJPSI <- function(jranges, exranges, jcounts, targets)
{
  bin.within <- findOverlaps(exranges, jranges, type="within")
  bin.within.df <- as.data.frame(bin.within)
  bin.within.df$names <- names(exranges[bin.within.df$queryHits])
  bin.within.df$J3 <- names(jranges[bin.within.df$subjectHits])
  start <- ncol(jcounts) -nrow(targets) +1
  end <- ncol(jcounts)
  df1<- data.frame(names=bin.within.df$names, 
                 jcounts[bin.within.df$subjectHits,start:end],
                 row.names=NULL) #recover counts
  tt1 <- data.frame(aggregate(. ~ names, data = df1, sum)) #aggregate counts
  tt2 <- data.frame(aggregate( J3 ~ names, data = bin.within.df, paste, collapse=";")) 
  #aggregate by name
  aa <- merge(tt2,tt1, by.x="names", by.y="names")#merge both dataframes
  rownames(aa) <- aa$names
  aa$names  <-  NULL
  return(aa)
}
###############################################################
.createGRangesExpJunctions <-
  function(jnames) 
  {
    split <- data.frame(matrix(unlist(strsplit(jnames,  "[.]") ),
                     byrow=TRUE, ncol=3),   
                     row.names=jnames,
                     stringsAsFactors=FALSE )  
    colnames(split)  <- c("chrom","start","end")
    split$start <- as.numeric(split$start)  
    split$end <- as.numeric(split$end)
    split$chrom <- as.character(split$chrom)
    jranges <- GRanges( seqnames = split$chrom, 
                        ranges =IRanges(
                        start= split$start, 
                        end = split$end,
                        names = jnames)    )
    jranges <- sort(jranges)
    return(jranges)
  }
###############################################################
.junctionsDiscover <- function(df, bam, cores, l, targets, features, pair) {
  countsSt <- ncol(df)-nrow(targets)+1
  jcounts <- df[,countsSt:ncol(df)] 
  jranges <- .createGRangesExpJunctions(rownames(df))
  e1i <- jranges
  minAnchor <- round(8*l/100)
  start(e1i) <- start(jranges)-(l-minAnchor)-1
  end(e1i) <- start(jranges)+(l-minAnchor)-1
  ie2 <- jranges
  start(ie2) <- end(jranges)-(l-minAnchor)-1
  end(ie2) <- end(jranges)+(l-minAnchor)-1
  #######################################################
  ungapped <- lapply(bam, function(x) {x[njunc(x)==0,]}) #extraigo los que no tienen GAPS
  hits.e1i <- lapply(ungapped,function(x){countOverlaps(e1i, x,  
                                                      ignore.strand = TRUE, 
                                                      minoverlap = l)})
  hits.e1i.ul <- do.call(cbind.data.frame, hits.e1i)
  hits.ie2 <- lapply(ungapped,function(x){countOverlaps(ie2, 
                                                      x,  
                                                      ignore.strand = TRUE, 
                                                      minoverlap = l)})
  hits.ie2.ul <- do.call(cbind.data.frame, hits.ie2)
  #######################################################
  dfPIR <- cbind(hits.e1i.ul, hits.ie2.ul, jcounts)
  cndfPIR <- colnames(dfPIR)
  #only counts, repeat sample names, they are used for factorization
  intPIR <- function(x){
    x[is.na(x)] <- 0 # we have to remove NAs
    res <- (x[1]+x[2])/(x[1]+x[2]+2*x[3])
    return(res)
  }
  ff <- rep(targets$condition,3)
  colnames(dfPIR) <- paste(ff, rep(1:3,each=length(targets$condition)), sep="."); head(dfPIR)
  #here colnames are converted into condition
  myfactor <- factor(paste(ff, rep(1:3,each=length(targets$condition)),sep="."), 
                   levels=paste(rep(pair, each=3 ), 1:3, sep="."))
  dfSum <- t(apply(dfPIR, 1, function(x){tapply(as.numeric(x),
                                              INDEX=myfactor, sum  )}))
  #sumo las junturas por condicion
  colnames(dfSum) <- rep(unique(targets$condition), each=3)
  myfactor <- factor(rep(unique(targets$condition),each=3), levels=pair)
  intPIRres <- t(apply(dfSum, 1, 
                function(x){tapply(as.numeric(x), 
                            INDEX=myfactor, 
                            intPIR  )}))
  #########################################
  hitIntron <- rep("-",  nrow(df))
  hitIntronEvent <- rep("-", nrow(df))
  #########################################
  intronBins <- featuresb(features)
  start(intronBins) <- start(intronBins)-1
  end(intronBins)<-end(intronBins)+1
  overJunctionEqualBins <- findOverlaps(jranges, intronBins,ignore.strand = TRUE, type="equal")
  overJunctionEqualBinsDF <- as.data.frame(overJunctionEqualBins); head(overJunctionEqualBinsDF)
  namesJ <- as.numeric(overJunctionEqualBinsDF[,1])
  namesB <- as.numeric(overJunctionEqualBinsDF[,2])
  overJunctionEqualBinsDF[,1] <- names(jranges[namesJ])
  overJunctionEqualBinsDF[,2] <- names(intronBins[namesB])
  overJunctionEqualBinsDF[,3] <- intronBins@elementMetadata$event[namesB]
  colnames(overJunctionEqualBinsDF)[3] <- "event"
  tw <- match(names(jranges), overJunctionEqualBinsDF$queryHits) #ok;
  hitIntron <- overJunctionEqualBinsDF$subjectHits[tw]
  hitIntronEvent <- overJunctionEqualBinsDF$event[tw]
  #############################################################
  colnames(dfPIR) <- cndfPIR
  pir_large <- data.frame(hitIntron, hitIntronEvent, dfPIR, intPIRres)
  fcoord <- paste(seqnames(jranges), 
                start(jranges),
                end(jranges) , sep=".")
  rownames(pir_large) <- fcoord
  return(pir_large)
}
######################################################################3
.junctionsPSI_SUM <- function(df, 
                            targets,
                            pair)
{
  jratio <- function(x){
    x[is.na(x)] <- 0 # we have to remove NAs
    res <- x[1]/(x[1]+x[2])
    return(res)
  }
  ############those sharing 3', 5', total etc#################
  countsSt <- ncol(df)-nrow(targets)+1
  jranges <- .createGRangesExpJunctions(rownames(df))
  
  if (packageVersion("IRanges")<2.6) 
  {
    j.start <- findOverlaps(jranges, ignoreSelf=TRUE,
                            ignoreRedundant=FALSE,type="start")
  }
  else
  {
    j.start <- findOverlaps(jranges, drop.self=TRUE,
                        drop.redundant=FALSE,type="start")
    
  }
  jjstart <- as.data.frame(j.start)
  jjstart$queryHits <- names(jranges[jjstart$queryHits])
  jjstart$subjectHits <- names(jranges[jjstart$subjectHits])
  shareStart <- data.frame(aggregate(subjectHits ~ queryHits, data = jjstart, paste, collapse=";")) 
  #counts matrix
  start <- ncol(df)-nrow(targets) +1 #ok
  end <- ncol(df)#ok
  #aca no hay cuentas solo recupero los counts usando el indexado de subjectHits
  dfCountsStart <- data.frame(
                            names=jjstart$queryHits, 
                            df[jjstart$subjectHits,start:end],
                            row.names=NULL) #recover counts
  #tiene como names al query y como counts todos los subjects que dan con ese query. 
  #en teoria podria haber mas de 1 sbject, por eso hace el aggregate-
  #si es solo 1 hit deberian coincidir, cruzandose
  #aca hace el agregate y tengo los counts de todas las junturas 
  #que comparten start con esa menos ella
  dfSumStart <- data.frame(aggregate(. ~ names, data = dfCountsStart, sum))
  sumJ <- paste(colnames(dfSumStart), "jsum", sep=".")
  colnames(dfSumStart) <- sumJ
  rownames(dfSumStart) <- dfSumStart$names.jsum
  dfSumStart$names.jsum <- NULL
  #armo un nuevo data frame
  dffStart <- data.frame(matrix(NA, nrow =  nrow(df), ncol = ncol(dfSumStart)) )
  rownames(dffStart) <- rownames(df)
  colnames(dffStart) <- colnames(dfSumStart)
  mSumStart <- match(row.names(dfSumStart), row.names(dffStart)) 
  #reordeno el dffSumStart de acuerdo al df
  dffStart[mSumStart,] <- dfSumStart#OK
  dffStart[is.na(dffStart)] <- 0
  mStartHit <- match(shareStart$queryHits, row.names(dffStart))
  #aca reacomodo el starthit con el indexado de dffStart
  StartHit <- as.character(rep("-", nrow(dffStart)) )
  StartHit[mStartHit] <- shareStart$subjectHits
  ################################################
  ratioStart <- data.frame(df[,countsSt:ncol(df)],
                         dffStart)
  colnames(ratioStart) <- rep(rownames(targets),2)
  #aca hay que armar un df itnermedio con la suma por condicion:
  ff <- factor(rep(targets$condition,2))
  colnames(ratioStart) <- factor(paste(ff, rep(1:2,each=length(targets$condition))), 
                               levels=paste(rep(pair, each=2 ), 1:2))
  myfactor <- factor(paste(ff, rep(1:2,each=length(targets$condition))), 
                   levels=paste(rep(pair, each=2 ), 1:2))
  ####
  dfSum <- t(apply(ratioStart, 1, function(x){tapply(as.numeric(x), 
                                                   INDEX=myfactor, sum  )}))
  myfactor <- factor(rep(unique(targets$condition),each=2), levels=pair)
  #######################################
  #solo este seria el df que hay que recuperar para una tabla resumida:
  jratioStartRes <- t(apply(dfSum, 1, 
                     function(x){tapply(as.numeric(x), 
                                             INDEX=myfactor, 
                                             jratio )}))
  colnames(jratioStartRes) <- paste(colnames(jratioStartRes), "start", sep=".")
  ################################################
  if (packageVersion("IRanges")<2.6) 
  {
    j.end=findOverlaps(jranges, 
                       ignoreSelf=TRUE,
                       ignoreRedundant=FALSE,
                       type="end")
    
  }
    else
  {
    j.end=findOverlaps(jranges, 
                     drop.self=TRUE,
                     drop.redundant=FALSE,
                     type="end")
  }
  jjend <- as.data.frame(j.end)
  jjend$queryHits <- names(jranges[jjend$queryHits])
  jjend$subjectHits <- names(jranges[jjend$subjectHits])
  shareEnd <- data.frame(aggregate(subjectHits ~ queryHits, data = jjend, paste, collapse=";")) 
  dfCountsEnd <- data.frame( names=jjend$queryHits, 
                          df[jjend$subjectHits,start:end],
                          row.names=NULL) #recover counts
  dfSumEnd <- data.frame(aggregate(. ~ names, data = dfCountsEnd, sum))   
  sumJ <- paste(colnames(dfSumEnd), "jsum", sep=".")
  colnames(dfSumEnd) <- sumJ
  rownames(dfSumEnd) <- dfSumEnd$names.jsum
  dfSumEnd$names.jsum <- NULL
  dffEnd <- data.frame(matrix(NA, nrow =  nrow(df), ncol = ncol(dfSumEnd)) )
  rownames(dffEnd) <- rownames(df)
  colnames(dffEnd) <- colnames(dfSumEnd)
  ########################################################################
  mSumEnd <- match(row.names(dfSumEnd), row.names(dffEnd))
  dffEnd[mSumEnd,] <- dfSumEnd
  dffEnd[is.na(dffEnd)] <- 0
  ########################################################################
  mEndHit <- match(shareEnd$queryHits, row.names(dffEnd))
  EndHit <- rep("-", nrow(dffEnd))
  EndHit[mEndHit] <- shareEnd$subjectHits
  ########################################################################
  ratioEnd <- data.frame(df[,countsSt:ncol(df)],dffEnd)
  ff <- rep(targets$condition,2)
  colnames(ratioEnd) <- paste(ff, rep(1:2,each=length(targets$condition)))
  myfactor <- factor(paste(ff, rep(1:2,each=length(targets$condition))), 
                   levels=paste(rep(pair, each=2 ), 1:2))
  dfSum <- t(apply(ratioEnd, 1, function(x){tapply(as.numeric(x), 
                                                 INDEX=myfactor,
                                                 sum  )}))
  colnames(dfSum) <- rep(unique(targets$condition),each=2)
  myfactor <- factor(rep(unique(targets$condition),each=2), levels=pair)
  #########################################################################
  jratioEndRes <- t(apply(dfSum, 1, function(x){tapply(as.numeric(x), 
                                                     INDEX=myfactor, 
                                                     jratio )}))
  colnames(jratioEndRes) <- paste(colnames(jratioEndRes), "end", sep=".")
  #########################################################################
  #eventos putativos:
  pAS <- rep("-", nrow(df))
  pAS[StartHit!="-" & df$strand =="-" | EndHit!="-" & df$strand =="+"] <- "Alt5ss"
  pAS[StartHit!="-" & df$strand =="+" | EndHit!="-" & df$strand =="-"] <- "Alt3ss"
  pAS[StartHit!="-" & EndHit!="-"] <- "ES" #OK
  ##########################################################################
  psi_merge_large <- data.frame(df,                       
                              StartHit,
                              dffStart,
                              jratioStartRes,
                              EndHit,
                              dffEnd,
                              jratioEndRes,
                              pAS)
  return(psi_merge_large)
  }                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
