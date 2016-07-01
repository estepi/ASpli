.filterByReads <-
    function(df0, targets, min, type) 
    {
      a <- table(targets$condition) 
      start <- ncol(df0)-sum(a)+1
      cn <- colnames(df0)[start:ncol(df0)]
      colnames(df0)[start:ncol(df0)] <- as.character(targets$condition)
      #independiente del numero de condiciones
      list <- matrix(unlist(
        lapply(
          unique(targets$condition), 
          function(x) rowMeans(df0[,colnames(df0) == x]) >= min  )),  
        nrow=nrow(df0), 
        byrow=FALSE)
      colnames(df0)[start:ncol(df0)] <- cn
      if (type == "all"){
        df=df0[rowSums(list)==ncol(list),]
      }else {
        df=df0[rowSums(list)>0,]
      }
      return (df)
    }
#########################################################################
.filterByRdGen <-
  function(df0, targets, min, type) 
  {
    a <- table(targets$condition) 
    start <- ncol(df0)-sum(a)+1
    frd <- df0[,start:ncol(df0)]/df0$effective_length
    colnames(frd) <- targets$condition
    ######################  
    list <- matrix(unlist(
      lapply(unique(colnames(frd)), 
              function(x) rowMeans(frd[,colnames(frd) == x])  >= min )),  
      nrow = nrow(frd), 
      byrow = FALSE)
    #######################
    #keeps those genes which ave rd > min in any condition
    if (type=="any")  { 
      ii <- rowSums(list)>0  
      df <- df0[ii,]
    }else  { 
      ii <- rowSums(list)==ncol(list)
      df <- df0[ii,]
    }
    return (df)
  }
#####################################################################
.genesDE <-
  function(df, targets, pair, group)
  { 
    countsSt <- ncol(df) - nrow(targets) + 1
    #########################################################
    er <- DGEList(counts=df[,countsSt:ncol(df)],group=group)
    er <- calcNormFactors(er) #methods TMM
    er <- estimateCommonDisp(er)
    er <- estimateTagwiseDisp(er)
    et <- exactTest(er, pair=pair)
    fdr.gen <- p.adjust(et$table$PValue, method="BH")
    et_merge <- cbind(et$table, fdr.gen) #ww should have read density by sample
    #########################################################
    genes_full <- data.frame(df[,-(countsSt:ncol(df))],
                           logFC=as.numeric(et$table$logFC), 
                           pvalue=as.numeric(et$table$PValue), 
                           gen.fdr=as.numeric(fdr.gen), 
                           stringsAsFactors=FALSE)
    rownames(genes_full) <- rownames(df)
    return(genes_full)
  }
################################################################
.filterByRdBinRATIO <-
  function(
    dfGen,
    dfBin,
    targets, 
    min, 
    type) 
  {
    ########################################
    a <- table(targets$condition) 
    startG <- ncol(dfGen)-sum(a)+1
    startB <- ncol(dfBin)-sum(a)+1
    genes.rd <- dfGen[,startG:ncol(dfGen)]/dfGen$effective_length #OK
    bins.rd <- dfBin[,startB:ncol(dfBin)]/dfBin$length #OK
    ###############################
    colnames(bins.rd) <- targets$condition
    avRdBin <- matrix(unlist(
      lapply(unique(colnames(bins.rd)), 
              function(x) 
              rowMeans(bins.rd[,colnames(bins.rd) == x]))),  
      nrow = nrow(bins.rd), 
      byrow = FALSE)#OK
    #################################################
    colnames(genes.rd) <- targets$condition
    avRdGen <- matrix(unlist(
      lapply(unique(colnames(genes.rd)), 
             function(x) rowMeans(genes.rd[,colnames(genes.rd) == x]))),  
      nrow = nrow(genes.rd), 
      byrow = FALSE)#Ok
    #####################################################
    te <- match(dfBin$locus, rownames(dfGen)) 
    gen.rdb <- avRdGen[te,]
    bin.gen.rd <-  avRdBin/gen.rdb
    bin.gen.rd[is.na(bin.gen.rd)] <- 0 
    
    if (type=="any")  { 
      ii  <- rowSums(bin.gen.rd >= min)>0
      dfBin <- dfBin[ii,]
    } else{ 
      ii <- rowSums(bin.gen.rd >= min )==ncol(bin.gen.rd)
      dfBin <- dfBin[ii,]
    }
    return (dfBin)
  }
##############################################################
.normalizeByGenFeature <-
  function(feature, gene, targets)
  {
    f.index <- match(feature$locus, rownames(gene)) #identify counts of each gene
    counts.genes.f <- gene[f.index,] 
    #have a matrix of gene counts, same dim as exons counts
    #no matter about the size of the matrix
    startg <-  ncol(counts.genes.f) - nrow(targets) +1
    endg <- ncol(counts.genes.f)
    
    startf <- ncol(feature) - nrow(targets) +1
    endf <- ncol(feature)
    gen.mean.f <- rowMeans(counts.genes.f[,startg:endg]) #get the mean of gene counts  
    #OK, independent of number of samples  
    #round the division  
    counts.f.N <- round(as.matrix(feature[,startf:endf]) / 
                        as.matrix(counts.genes.f[,startg:endg]) *gen.mean.f) #round the division
    counts.f.N[is.na(counts.f.N)] <- 0 # we have to remove NAs
    counts.f.N[is.infinite(counts.f.N)] <- 0 #we also have to remove Inf values
    #correct
    feature.counts.n<-cbind(feature[,1:startf -1], counts.f.N)
    #have a df equal to original df counts (without norm)
    return(feature.counts.n)
  }
#################################################################
.binsDU <-
  function(df,
           targets, 
           dfGen,
           ignoreExternal=NULL, 
           pair,
           group)
  {
    df <- .normalizeByGenFeature(feature=df, gene=dfGen, targets) #OK  
    if (ignoreExternal==TRUE)
    {
      df=df[df$event!="external",] 
    }
    countsSt <- ncol(df)-nrow(targets)+1
    er <- DGEList(counts=df[,countsSt:ncol(df)], 
                  group=group)
    er <- calcNormFactors(er)
    er <- estimateCommonDisp(er)
    er <- estimateTagwiseDisp(er)
    et <- exactTest(er, pair=pair)
    fdr.bin <- p.adjust(et$table$PValue, method="BH")
    
    logFC <- et$table$logFC
    pvalue <- et$table$PValue
    bin.fdr <- fdr.bin
    splicing_full <- data.frame (df[,-(countsSt:ncol(df))],
                                logFC, 
                                pvalue,
                                bin.fdr,
                                stringsAsFactors=FALSE)
    rownames(splicing_full) <- rownames(df)
    return(splicing_full)
  }
##################################################################
.junctionsDU_SUM <- function(df, 
                          targets, 
                          genesde, 
                          pair, 
                          group,
                          dfGen)  
{
  ######################################################
    jratio<-function(x){
    x[is.na(x)] <- 0 # we have to remove NAs
    res<-x[1]/(x[1]+x[2])
    return(res)  }
  ##############################################################################################################3  
  #normalize:
  colnames(df)[2] <- "locus"
  df<-.normalizeByGenFeature(feature=df, gene=dfGen, targets) #OK  
  colnames(df)[2] <- "gene"
  countsSt <- ncol(df) - nrow(targets) + 1
  er <- DGEList(counts=df[,countsSt:ncol(df)],
                group=group)
  er <- calcNormFactors(er)
  er <- estimateCommonDisp(er)
  er <- estimateTagwiseDisp(er)
  et <- exactTest(er, pair=pair)
  fdr <- p.adjust(et$table$PValue, method="BH")
  logFC <- et$table$logFC
  pvalue <- et$table$PValue
  ############those sharing 3', 5', total etc#################
  jranges <- .createGRangesExpJunctions(rownames(df))
  if (packageVersion("IRanges")<2.6) 
  {
    j.start <- findOverlaps(jranges, ignoreSelf=TRUE, ignore.redundant=FALSE,
                            type="start")  
  }
  else
  {
    
  j.start <- findOverlaps(jranges, drop.self=TRUE, drop.redundant=FALSE,
                          type="start")
  }
  jjstart <- as.data.frame(j.start)
  jjstart$queryHits <- names(jranges[jjstart$queryHits])
  jjstart$subjectHits <- names(jranges[jjstart$subjectHits])
  shareStart <- data.frame(aggregate(subjectHits ~ queryHits, 
                                     data = jjstart, paste, collapse=";")) 
  #counts matrix
  start <- ncol(df) - nrow(targets) + 1 #ok
  end <- ncol(df)#ok
  #aca no hay cuentas solo recupero los counts usando el indexado de subjectHits
  dfCountsStart <- data.frame( names=jjstart$queryHits, 
                            df[jjstart$subjectHits,start:end],
                            row.names=NULL) #recover counts
  #tiene como names al query y como counts todos los subjects que dan con ese query. 
  #en teoria podria haber mas de 1 sbject, por eso hace el aggregate-
  #si es solo 1 hit deberian coincidir, cruzandose
  #aca hace el agregate y tengo los counts de todas las junturas 
  #que comparten start con esa menos ella
  dfSumStart <- data.frame(aggregate(. ~ names, data = dfCountsStart, sum))
  sumJ <- paste(colnames(dfSumStart), "jsum", sep=".")
  colnames(dfSumStart) <-  sumJ
  rownames(dfSumStart) <- dfSumStart$names.jsum
  dfSumStart$names.jsum <- NULL
  #armo un nuevo data frame
  dffStart <- data.frame(matrix(NA, nrow =  nrow(df), ncol = ncol(dfSumStart)) )
  rownames(dffStart) <- rownames(df)
  colnames(dffStart) <- colnames(dfSumStart)
  mSumStart <- match(row.names(dfSumStart), row.names(dffStart)) 
  #reordeno el dffSumStart de acuerdo al df
  dffStart[mSumStart,] <-  dfSumStart#OK
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
  dfSum <- t(apply(ratioStart, 1, function(x){tapply(as.numeric(x), 
                                                   INDEX=colnames(ratioStart), 
                                                   sum)}))
  colnames(dfSum) <- rep(unique(targets$condition), each=2)
  
  jratioStartRes <- t(apply(dfSum, 1, function(x){tapply(as.numeric(x), 
                                                       INDEX=colnames(dfSum), 
                                                       jratio )}))
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
  dffEnd =data.frame(matrix(NA, nrow =  nrow(df), 
                            ncol = ncol(dfSumEnd)) )
  rownames(dffEnd) <- rownames(df)
  colnames(dffEnd) <- colnames(dfSumEnd)
  ########################################################################
  mSumEnd <- match(row.names(dfSumEnd), row.names(dffEnd))
  dffEnd[mSumEnd,] <- dfSumEnd
  dffEnd[is.na(dffEnd)] <- 0
  ########################################################################
  mbin_end_hit <- match(shareEnd$queryHits, row.names(dffEnd))
  bin_end_hit <- rep("-", nrow(dffEnd))
  bin_end_hit[mbin_end_hit] <- shareEnd$subjectHits
  ########################################################################
  ratioEnd <- data.frame(df[,countsSt:ncol(df)],dffEnd)
  ff <- rep(targets$condition,2)
  colnames(ratioEnd) <- paste(ff, rep(1:2,length(targets$condition)))
  dfSum <- t(apply(ratioEnd, 1, function(x){tapply(as.numeric(x), 
                                                 INDEX=colnames(ratioEnd), sum  )}))
  colnames(dfSum) <- rep(unique(targets$condition),each=2)
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
                       jratioEndRes )
return(et_merge)
}
##################################################
