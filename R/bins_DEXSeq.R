.binsDU_DEXSeq <-
  function(df,targets, group)
  {
    countsSt <- ncol(df)-nrow(targets)+1
    GeneInfo <- data.frame(Gene.Exon=rownames(df), 
                         matrix(unlist(strsplit(as.character(rownames(df)), ":")), 
                                nrow=nrow(df), ncol=2, byrow=TRUE ) )
    colnames(GeneInfo) = c("GeneExon","GeneID","ExonID")
    countData <- df[,countsSt:ncol(df)]
    sampleData <- data.frame(condition=group)
    design1 <- formula( ~ sample + exon + condition:exon )
    groupID <- GeneInfo$GeneID
    featureID <- GeneInfo$GeneExon
    data <- DEXSeqDataSet( countData, sampleData, design1,featureID, groupID )
    dx <- estimateSizeFactors( data )
    dx <- estimateDispersions( dx)
    dx <- testForDEU( dx )
    dx <- DEXSeqResults( dx)
    rownames(dx) <- rownames(df)
    logFC <- dx$exonBaseMean
    pvalue <- dx$pvalue
    bin.fdr <- dx$padj
  
    splicing_full <- data.frame (df[,-(countsSt:ncol(df))],
                                logFC, 
                                pvalue,
                                bin.fdr,
                                stringsAsFactors=FALSE)
    rownames(splicing_full) <- rownames(df)
    return(splicing_full)
  }
##################################################################
.junctionsDU_SUM_DEXSeq <- function(df, 
                           targets, 
                           genesde, group)
  
{
  ######################################################
  jratio<-function(x){
    x[is.na(x)] <- 0 # we have to remove NAs
    res <- x[1]/(x[1]+x[2])
    return(res)  }
  ##############################################################################################################3  
  countsSt <- ncol(df)-nrow(targets)+1
  GeneInfo <- data.frame(Gene.Exon=rownames(df), 
                       matrix(unlist(strsplit(as.character(rownames(df)), ":")), 
                       nrow=nrow(df), ncol=2, byrow=TRUE ) )
  colnames(GeneInfo) <- c("GeneExon","GeneID","ExonID")
  countData <- df[,countsSt:ncol(df)]
  sampleData <- data.frame(condition=group)
  design1 <- formula( ~ sample + exon + condition:exon )
  groupID <- GeneInfo$GeneID
  featureID <- GeneInfo$GeneExon
  data <- DEXSeqDataSet( countData, sampleData, design1,featureID, groupID )
  data <- estimateSizeFactors( data )
  dispersions(data) <- 0.1
  dx <- testForDEU( data )
  dx <- DEXSeqResults( dx)
  rownames(dx) <- rownames(df)
  logFC <- dx$exonBaseMean
  pvalue <- dx$pvalue
  fdr <- dx$padj
  ############those sharing 3', 5', total etc#################
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
  shareStart <- data.frame(aggregate(subjectHits ~ queryHits, 
                                     data = jjstart, paste, collapse=";")) 
  #counts matrix
  start <- ncol(df) - nrow(targets) +1 #ok
  end <- ncol(df)#ok
  #aca no hay cuentas solo recupero los counts usando el indexado de subjectHits
  dfCountsStart <- data.frame( names=jjstart$queryHits, 
                            df[jjstart$subjectHits,start:end],
                            row.names=NULL) #recover counts
  #tiene como names al query y como counts todos los subjects que dan con ese query. 
  #en teoria podria haber mas de 1 sbject, por eso hace el aggregate-
  #si es solo 1 hit deberian coincidir, cruzandose
  #aca hace el agregate y tengo los counts de todas las junturas que comparten start con esa menos ella
  dfSumStart <- data.frame(aggregate(. ~ names, data = dfCountsStart, sum))
  sumJ <- paste(colnames(dfSumStart), "jsum", sep=".")
  colnames(dfSumStart) <- sumJ
  rownames(dfSumStart) <-  dfSumStart$names.jsum
  dfSumStart$names.jsum <- NULL
  #armo un nuevo data frame
  dffStart <- data.frame(matrix(NA, nrow =  nrow(df), 
                              ncol = ncol(dfSumStart)) )
  rownames(dffStart) <- rownames(df)
  colnames(dffStart) <- colnames(dfSumStart)
  mSumStart <- match(row.names(dfSumStart), row.names(dffStart)) 
  #reordeno el dffSumStart de acuerdo al df
  dffStart[mSumStart,] <- dfSumStart#OK
  dffStart[is.na(dffStart)] <- 0
  mbin_start_hit <- match(shareStart$queryHits, row.names(dffStart))
  #aca reacomodo el bin_start_hit con el indexado de dffStart
  bin_start_hit <- as.character(rep("-", nrow(dffStart)) )
  bin_start_hit[mbin_start_hit] <- shareStart$subjectHits
  ################################################
  ratioStart <- data.frame(df[,countsSt:ncol(df)],
                         dffStart)
  colnames(ratioStart) <- rep(rownames(targets),2)
  #aca hay que armar un df itnermedio con la suma por condicion:
  ff <- rep(targets$condition,2)
  colnames(ratioStart) <- paste(ff, rep(1:2,each=length(targets$condition)))
  dfSum <- t(apply(ratioStart, 1, function(x){tapply(as.numeric(x), INDEX=colnames(ratioStart), sum  )}))
  colnames(dfSum) <- rep(unique(targets$condition),each=2)
  jratioStartRes <- t(apply(dfSum, 1, function(x){tapply(as.numeric(x), INDEX=colnames(dfSum), jratio )}))
  ################################################
  if (packageVersion("IRanges")<2.6) 
  {
    j.end <- findOverlaps(jranges, 
                       ignoreSelf=TRUE,
                       ignoreRedundant=FALSE,
                       type="end")
  }
  else
  {
  j.end <- findOverlaps(jranges, 
                     drop.self=TRUE,
                     drop.redundant=FALSE,
                     type="end")
  }
  jjend <- as.data.frame(j.end)
  jjend$queryHits <- names(jranges[jjend$queryHits])
  jjend$subjectHits <- names(jranges[jjend$subjectHits])
  shareEnd <- data.frame(aggregate(subjectHits ~ queryHits,
                                data = jjend, paste, collapse=";")) 
  dfCountsEnd <- data.frame( names=jjend$queryHits, 
                          df[jjend$subjectHits,start:end],
                          row.names=NULL) #recover counts
  dfSumEnd <- data.frame(aggregate(. ~ names, data = dfCountsEnd, sum))   
  sumJ <- paste(colnames(dfSumEnd), "jsum", sep=".")
  colnames(dfSumEnd) <- sumJ
  rownames(dfSumEnd) <- dfSumEnd$names.jsum
  dfSumEnd$names.jsum <- NULL
  dffEnd  <- data.frame(matrix(NA, nrow =  nrow(df), ncol = ncol(dfSumEnd)) )
  rownames(dffEnd) <- rownames(df)
  colnames(dffEnd) <- colnames(dfSumEnd)
  ########################################################################
  mSumEnd <- match(row.names(dfSumEnd), row.names(dffEnd))
  dffEnd[mSumEnd,] <-  dfSumEnd
  dffEnd[is.na(dffEnd)] <- 0
  ########################################################################
  mbin_end_hit <- match(shareEnd$queryHits, row.names(dffEnd))
  bin_end_hit <- rep("-", nrow(dffEnd))
  bin_end_hit[mbin_end_hit] <- shareEnd$subjectHits
  ########################################################################
  ratioEnd <- data.frame(df[,countsSt:ncol(df)],dffEnd)
  ff <- rep(targets$condition,2)
  colnames(ratioEnd) <- paste(ff, rep(1:2,each=length(targets$condition)))
  dfSum <- t(apply(ratioEnd, 1, function(x){tapply(as.numeric(x), 
                                                 INDEX=colnames(ratioEnd), sum  )}))
  colnames(dfSum) <- rep(levels(targets$condition),each=2)
  jratioEndRes <- t(apply(dfSum, 1, function(x){tapply(as.numeric(x), 
                                                     INDEX=colnames(dfSum), jratio )}))
  #########################################################################
  et_merge <- data.frame(df,                       
                       logFC,
                       pvalue, 
                       fdr,
                       bin_start_hit,
                       dffStart,
                       jratioStartRes,
                       bin_end_hit,
                       dffEnd,
                       jratioEndRes)
  return(et_merge)
}                        
##################################################
.genesDE_DESeq <-
  function(df, targets, pair)
  { 
    countsSt <- ncol(df)-nrow(targets)+1
    colData <- data.frame(condition=targets[,2], row.names=row.names(targets))
    #########################################################
    dds <- DESeqDataSetFromMatrix(countData = df[,countsSt:ncol(df)],
                                  colData = colData,
                                  design = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds, contrast=c("condition",pair))
    fdr.gen <- res$padj
    #########################################################
    genes_full <- data.frame(df[,-(countsSt:ncol(df))],
                           logFC=as.numeric(res$log2FoldChange), 
                           pvalue=as.numeric(res$pvalue), 
                           gen.fdr=as.numeric(fdr.gen), 
                           stringsAsFactors=FALSE)
    rownames(genes_full) <- rownames(df)
    return(genes_full)
  }
##################################################
