##########################################################################
setClass(
    Class="ASpli-features",
    representation=representation(
      genes = "GRangesList",
      bins="GRanges",
      junctions="GRanges"))
##########################################################################
setClass(
  Class="ASpli-counts",
  representation=representation(
    gene.counts= "data.frame", 
    exon.intron.counts= "data.frame",
    junction.counts= "data.frame",
    e1i.counts= "data.frame", 
    ie2.counts= "data.frame",
    gene.rd="data.frame",
    bin.rd="data.frame"))
##########################################################################
setClass(
  Class="ASpli-AS",
  representation=representation(
    irPIR = "data.frame",
    altPSI="data.frame",
    esPSI="data.frame",
    junctionsPIR="data.frame",
    junctionsPSI="data.frame",
    join="data.frame")
  )
##########################################################################
  setClass(
    Class="ASpli-DU",
    representation=representation(
      genes = "data.frame",
      bins="data.frame",
      junctions="data.frame"))
##########################################################################  
setGeneric (
    name= "binGenome",
    def=function(genome, md=NULL)
      {standardGeneric("binGenome")})
##########################################################################
setMethod(
  f= "binGenome",
  signature= "TxDb",
  definition=function (genome,md=NULL){
    
    features<-new(Class="ASpli-features")
    
    if (is.null(md))
    {
      md<-data.frame(names(transcriptsBy(genome)), 
      row.names=names(transcriptsBy(genome)), 
      stringsAsFactors=FALSE )
      colnames(md)="symbol" 
    }
    
    genes.by.exons<-.createGRangesGenes(genome, md) 
    ll<-length(genes.by.exons)
    message("* Number of extracted Genes =",ll,"\n")
##########################################################################
    exon.bins<-.createGRangesExons(genome, md)
    #add locus_overlap
    index<-match(exon.bins@elementMetadata$locus, names(genes.by.exons))
    locus_overlap<-rep("-", length(exon.bins))
    locus_overlap<-genes.by.exons@elementMetadata$locus_overlap[index]
    mcols(exon.bins) <- append(mcols(exon.bins), DataFrame(locus_overlap=locus_overlap))
    ll<-length(exon.bins)
    message("* Number of extracted Exon Bins =",ll)
##########################################################################
    intron.tot<-.createGRangesIntrons(genome, md)
    #add locus_overlap
    index<-match(intron.tot@elementMetadata$locus, names(genes.by.exons))
    locus_overlap<-rep("-", length(intron.tot))
    locus_overlap<-genes.by.exons@elementMetadata$locus_overlap[index]
    mcols(intron.tot) <- append(mcols(intron.tot), 
                                DataFrame(locus_overlap=locus_overlap))
    ll<-length(intron.tot)
    message("* Number of extracted intron bins =",ll)
##########################################################################
    transcripts<-.createGRangesTranscripts(genome)
    ll<-length(unlist(transcripts))
    message("* Number of extracted trascripts =",ll)
##########################################################################
    junctions<-.createGRangesJunctions(genome) 
    #add locus_overlap
    index<-match(junctions@elementMetadata$locus, names(genes.by.exons))
    locus_overlap<-rep("-", length(junctions))
    locus_overlap<-genes.by.exons@elementMetadata$locus_overlap[index]
    mcols(junctions) <- append(mcols(junctions), 
                               DataFrame(locus_overlap=locus_overlap))
    ll<-length(junctions)
    message("* Number of extracted junctions =",ll)
##########################################################################
    intron.bins<- intron.tot[intron.tot@elementMetadata$feature == "I"]
    intron.orig <-intron.tot[intron.tot@elementMetadata$feature == "Io"]
    class<-rep("fullI", length(intron.orig))
    mcols(intron.orig) <- append(mcols(intron.orig), DataFrame(class=class))
    event<-rep("-", length(intron.orig))
    eventJ<-rep("-", length(intron.orig))
    mcols(intron.orig) <- append(mcols(intron.orig), DataFrame(event=event))
    mcols(intron.orig) <- append(mcols(intron.orig), DataFrame(eventJ=eventJ))
    ########aca la modificacion:
    exons.introns<-.findASBin(exon.bins, intron.bins, transcripts, junctions)
    fullT<-c(exons.introns,intron.orig)
    #add locus_overlap #queda agregado por la combinacion
    ############aca habria que armar la lista.
    features@genes<-genes.by.exons
    features@bins<-fullT
    features@junctions<-junctions
    return(features)})
##########################################################################
setGeneric (
  name= "rds",
  def=function(counts, targets)
  {standardGeneric("rds")})
##########################################################################
setMethod(
  f= "rds",
  signature= "ASpli-counts",
  definition=function(counts, targets)
  {
    geneStart<-ncol(countsg(counts))-nrow(targets)+1
    gene.rd<-cbind(countsg(counts)[,1:geneStart-1], 
                   countsg(counts)[,geneStart:ncol(countsg(counts))]/countsg(counts)$effective_length
                   )
    binStart<-ncol(countsb(counts))-nrow(targets)+1
    bin.rd<-cbind(countsb(counts)[, 1:binStart-1], 
                  countsb(counts)[,binStart:ncol(countsb(counts))]
                  /countsb(counts)$length)
                  
    tb<-match(bin.rd$locus, rownames(gene.rd))
    rdfinalb=cbind(bin.rd, bin.rd[,binStart:ncol(bin.rd)]
                   /gene.rd[tb,geneStart:ncol(countsg(counts))])
                            
    counts@gene.rd=gene.rd
    counts@bin.rd=rdfinalb
    return(counts)
  }
)
##############################################################################
setGeneric (
  name="readCounts",
  def=function(features, bam, cores=NULL, l, maxISize, minAnchor=NULL)
  {standardGeneric("readCounts")})
##########################################################################
setMethod(
  f= "readCounts",
  signature= "ASpli-features",
  definition=function(features, bam, cores=NULL, l, maxISize, minAnchor=NULL)
  {
      counts<-new(Class="ASpli-counts")
      if(is.null(minAnchor)){minAnchor=10}
      minA=round(minAnchor*l/100)
      gene.hits=.counterGenes(bam, featuresg(features), cores)
      message("read summarization by gene completed")
      counts@gene.counts=gene.hits
      ####################################################################
      #exons
      bins<-featuresb(features)
      totCounts<-.counterBin(bam, bins, gene.hits, cores)
      message("read summarization by bin completed")
      ####################################################################
      #introns
      introns<-c(bins[bins@elementMetadata$feature == "I"], 
                 bins[bins@elementMetadata$feature == "Io" ],
                 bins[bins@elementMetadata$eventJ=="IR"])
      ####################################################################
      e1i<-introns
      start(e1i)<-start(introns)-(l-minAnchor)
      end(e1i)<-start(introns)+(l-minAnchor)
      e1i.hits<-.counterJbin(bam, e1i, gene.hits, cores, l)
      message("read summarization by ei1 region completed")
      #####################################################################
      ie2<-introns
      start(ie2)<-end(introns)-(l-minA)
      end(ie2)<-end(introns)+(l-minA)
      ie2.hits<-.counterJbin(bam, ie2, gene.hits, cores, l)
      message("read summarization by ie2 region completed")
      ####################################################################
       #junctions
      junction.hits = .counterJunctions(features, bam, cores, maxISize)
      message("junction summarization completed")
      #####################################################################
      counts@gene.counts=gene.hits
      counts@exon.intron.counts=totCounts
      counts@junction.counts=junction.hits 
      counts@e1i.counts=e1i.hits  
      counts@ie2.counts=ie2.hits
      counts<-rds(counts, targets)
      return(counts)
})
##########################################################################
setGeneric (
  name= "writeCounts",
  def=function(counts, output.dir="counts")
  {standardGeneric("writeCounts")})
##########################################################################
setMethod(
  f= "writeCounts",
  signature= "ASpli-counts",
  definition=function(counts, output.dir="counts")
  {
   currentDir<-getwd()
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
  
  file<-paste(genesDir, "gene.counts.tab", sep="/")
  write.table(countsg(counts), file, sep="\t", quote=FALSE, col.names=NA)
##########################################################################
  ec<-countsb(counts)[countsb(counts)$feature=="E",]
  ec<-ec[ec$event != "IR",]
  ec<-ec[ec$event != "IR*",]
  file<-paste(exonsDir, "exon.counts.tab", sep="/")
  write.table(ec, file, sep="\t", quote=FALSE, col.names=NA)
################ INTRONS ##################################
   file<-paste(intronsDir, "intron.counts.tab", sep="/")
   ic<-rbind(countsb(counts)[countsb(counts)$feature == "I",], 
             countsb(counts)[countsb(counts)$feature == "Io",], 
             countsb(counts)[countsb(counts)$event == "IR",],
             countsb(counts)[countsb(counts)$event == "IR*",])
    write.table(ic, file, sep="\t", quote=FALSE, col.names=NA)
    file<-paste(intronsDir, "e1i.counts.tab", sep="/")
    write.table(countse1i(counts), file, sep="\t", quote=FALSE, col.names=NA)
    file<-paste(intronsDir, "ie2.counts.tab", sep="/")
    write.table(countsie2(counts), file, sep="\t", quote=FALSE, col.names=NA)
    file<-paste(junctionsDir, "junction.counts.tab", sep="/" )
    write.table(countsj(counts), file, sep="\t", quote=FALSE,  col.names=NA)
  }
)
##########################################################################
setGeneric (
  name= "writeRds",
  def=function(counts, output.dir="rds")
  {standardGeneric("writeRds")})
##########################################################################
setMethod(
  f= "writeRds",
  signature= "ASpli-counts",
  definition=function(counts, output.dir="rds")
  {
    currentDir<-getwd()
    if (file.exists(output.dir))
    {
      outputDir <- paste(currentDir, output.dir, sep = "/")       
      genesDir <- paste(outputDir, "genes", sep = "/")       
      exonsDir <- paste(outputDir, "exons", sep = "/")       
      intronsDir <- paste(outputDir, "introns", sep = "/")
    }
    else 
    {
      outputDir <- paste(currentDir, output.dir, sep = "/")       
      genesDir <- paste(outputDir, "genes", sep = "/")       
      exonsDir <- paste(outputDir, "exons", sep = "/")	     
      intronsDir <- paste(outputDir, "introns", sep = "/")
      dir.create(outputDir)
      dir.create(genesDir)
      dir.create(exonsDir)
      dir.create(intronsDir)
    }
    file<-paste(genesDir, "gene.rd.tab", sep="/")
    write.table(rdsg(counts), file, sep="\t", quote=FALSE, col.names=NA)
    ################ EXONS ####################################
    erd<-rdsb(counts)[rdsb(counts)$feature == "E",]
    erd<-erd[erd$event != "IR",]
    erd<-erd[erd$event != "IR*",]
    
    file<-paste(exonsDir, "exon.rd.tab", sep="/")
    write.table(erd, file, sep="\t", quote=FALSE, col.names=NA)
    ################ INTRONS ##################################
    file<-paste(intronsDir, "intron.rd.tab", sep="/")
    ird<-rbind(rdsb(counts)[rdsb(counts)$feature == "I",], 
              rdsb(counts)[rdsb(counts)$feature == "Io",], 
              rdsb(counts)[rdsb(counts)$eventJ == "IR",])
    write.table(ird, file, sep="\t", quote=FALSE, col.names=NA)
    }
)
##########################################################################
setGeneric (
  name= "AsDiscover",
  def=function(counts, 
               targets, 
               features, 
               bam, 
               l, 
               pair,
               threshold=NULL, 
               cores=NULL)
  {
    standardGeneric("AsDiscover")
  }
  )
##########################################################################
setMethod(
    f= "AsDiscover",
  signature= "ASpli-counts",
  definition=function(counts, 
               targets, 
               features, 
               bam, 
               l, 
               pair,
               threshold=NULL, 
               cores=NULL)
    
  {
    as<-new(Class="ASpli-AS")
    if(is.null(threshold)){threshold=5}
    ### Common Part ############################################
    df0<-countsj(counts)[countsj(counts)$multipleHit=="-",]
    df0<-df0[df0$gene!="noHit",]
    jcounts<-.filterJunctionBySample(df0=df0, 
                                    targets=targets, 
                                    threshold=threshold)  #mean > one of the condition
##########################################################################
    #Junctions PSI:
    junctionsPSI<-.junctionsPSI_SUM(df0, targets, pair)
    as@junctionsPSI<-junctionsPSI
    message("Junctions PSI completed")
    ############es el mismo dataset#########################
    #JUNTURAS#######################
    #We analyze each junction as a putative intron. 
    #So, we estimate PIR using e1i, ie2 counts as an intron retention
    junctionsPIR<-.junctionsDiscover(df=jcounts, 
                                     bam, 
                                     cores=NULL, 
                                     l, 
                                     targets, 
                                     features, 
                                     pair)
    message("Junctions PIR completed")
    as@junctionsPIR<-junctionsPIR
    ############es el mismo dataset#########################
    split=data.frame(matrix(unlist(strsplit(rownames(jcounts),  "[.]") ),
                            byrow=TRUE, ncol=3),   
                            row.names=rownames(jcounts),
                            stringsAsFactors=FALSE )  
    colnames(split)=c("chrom","start","end")
    split$start<-as.numeric(split$start)
    split$end<-as.numeric(split$end)
    split$chrom<-as.character(split$chrom)
    jranges <- GRanges( 
    seqnames = split$chrom,
    ranges =IRanges(start= split$start, 
                      end = split$end, 
                      names = rownames(split)))
    ########################################################
    #parte 1: irPIR intronJPIR_SUM_COND
    ic<-rbind(countsb(counts)[countsb(counts)$feature == "I",], 
              countsb(counts)[countsb(counts)$feature == "Io",], 
              countsb(counts)[countsb(counts)$event == "IR*",],
              countsb(counts)[countsb(counts)$event == "IR",])
    intranges<-featuresb(features)[rownames(ic)]
##########################################################################
    dfe1e2<-.e1e2JPIR(intranges, jcounts)
    dfe1e2$J3<-as.character(dfe1e2$J3)
    #junction spanning intron, from e1 to e2. Deafault junction
    #in functions...
##########################################################################
    countsSt<-ncol(counts@e1i.counts)-nrow(targets)+1
    e1i=counts@e1i.counts[,countsSt:ncol(counts@e1i.counts)]; 
    ie2=counts@ie2.counts[,countsSt:ncol(counts@ie2.counts)] 
##########################################################################
    J1<-paste(rownames(e1i),"E1I", sep="_")
    J2<-paste(rownames(ie2),"IE2", sep="_")
    dffinal<-data.frame(J1, e1i, J2, ie2)
    #no tiene junturas, agregamos las columnas J1, J2
    empty=data.frame(matrix(NA, 
                            nrow =  nrow(dffinal), 
                            ncol =  ncol(dfe1e2)), 
                            stringsAsFactors=FALSE)  
    colnames(empty)=colnames(dfe1e2)
    int<-match(dfe1e2$jbin, row.names(dffinal))
    empty[int,]<-dfe1e2
    dmerge<-cbind(dffinal, empty); 
    dmerge$jbin=NULL
    rownames(dmerge)=names(intranges)
##########################################################################    
    intPIR<-function(x){
      x[is.na(x)] <- 0 # we have to remove NAs
      res<-(x[1]+x[2])/(x[1]+x[2]+2*x[3])
      return(res)
    }
##########################################################################
    df<-dmerge
    df$J1=NULL
    df$J2=NULL
    df$J3=NULL
    ff<-rep(targets$condition,3)
    myfactor<-factor(paste(ff, rep(1:3,each=length(targets$condition)),sep="."), 
                     levels=paste(rep(pair, each=3 ), 1:3, sep="."))
    colnames(df)<-paste(ff, rep(1:3,each=length(targets$condition)), sep=".")
    dfSum<-t(apply(df, 1, function(x){tapply(as.numeric(x), 
                                             INDEX=myfactor, sum  )}))
    colnames(dfSum)<-rep(unique(targets$condition), each=3) 
    myfactor<-factor(rep(unique(targets$condition),each=3), levels=pair)
    intPIRres<-cbind(event=ic$event,
                     dmerge, 
                     t(apply(dfSum, 1, 
                             function(x){tapply(as.numeric(x), 
                             INDEX=myfactor, 
                              intPIR  )})))
    message("Junctions IR PIR completed")
    as@irPIR<- intPIRres
##########################################################################
    #parte 2:altPSI, esPSI exonJPSI_SUM_COND
    ec<-countsb(counts)[countsb(counts)$feature == "E",]
    ec<-ec[ec$event != "IR",]
    ec<-ec[ec$event != "IR*",]
    exranges<-featuresb(features)[rownames(ec)]; dim(ec)
    dfstart<-.startJPSI(jranges, exranges, jcounts, targets) 
    #identifies start junctio
    dfend<-.endJPSI(jranges, exranges, jcounts, targets)  #identifies end junction 
    dfwithin<-.withinJPSI(jranges, exranges, jcounts, targets)#identifies within junction
    #add the information at the end of the DF of exonsde
    corr.start=match(rownames(dfstart), names(exranges))
    corr.end=match(rownames(dfend), names(exranges))
    corr.w=match(rownames(dfwithin), names(exranges))
##########################################################################
    dffinal=data.frame(matrix(NA, 
                              nrow =  length(exranges), 
                              ncol = 3*ncol(dfstart)), 
                       stringsAsFactors=FALSE)
    colnames(dffinal)=c(colnames(dfstart), colnames(dfend), colnames(dfwithin))
    colN<-rep(c("J", rownames(targets)), 3)
    colnames(dffinal)<-colN
    rownames(dffinal)=names(exranges)
    startStart<-1
    startEnd<-ncol(dfstart)
    endStart<-startEnd +1
    endEnd<-endStart+ ncol(dfend)-1
    withinStart<-endEnd+1
    withinEnd<-  withinStart +ncol(dfwithin)-1
##########################################################################
    dffinal[corr.start, startStart:startEnd]= dfstart
    dffinal[corr.end, endStart:endEnd]= dfend
    dffinal[corr.w, withinStart:withinEnd]= dfwithin
    cN<-colnames(dffinal)
    cN<-cN[cN!="J"]
    dfEvents<-cbind(event=as.data.frame(exranges@elementMetadata$event), dffinal)
    dim(dfEvents)#6885
    colnames(dfEvents)[1]<-"event"
##########################################################################
    dfAlt<-rbind(dfEvents[dfEvents$event=="Alt3ss",], 
                 dfEvents[dfEvents$event=="Alt5ss",], 
                 dfEvents[dfEvents$event=="Alt3ss*",], 
                 dfEvents[dfEvents$event=="Alt5ss*",])
    dfAltSub<-dfAlt[,colnames(dfAlt)!="J"]
    dfAltSub<-dfAltSub[,colnames(dfAltSub)!="event"]
    colnames(dfAltSub)  <-cN
##########################################################################
    altPSI<-function(x){
      x[is.na(x)] <- 0 # we have to remove NAs
      res<-(x[1]+x[2])/(x[1]+x[2]+x[3])
      return(res)
    }
##########################################################################
    ff<-rep(targets$condition,3)
    colnames(dfAltSub)<-paste(ff, rep(1:3,each=length(targets$condition)))
    myfactor<-factor(paste(ff, rep(1:3,each=length(targets$condition))), 
                     levels=paste(rep(pair, each=3 ), 1:3))
    
    dfAltSum<-t(apply(dfAltSub, 1, 
                      function(x){tapply(as.numeric(x), 
                                         INDEX=myfactor,
                                         sum  )}))
    
    colnames(dfAltSum)<-rep(unique(targets$condition), each=3)
    myfactor<-factor(rep(unique(targets$condition),each=3), levels=pair)
    altPSIRes<-cbind(dfAlt, t(apply(dfAltSum, 1, 
                                    function(x){tapply(as.numeric(x), 
                                                       INDEX=myfactor, 
                                                       altPSI  )})))
    colnames(altPSIRes)<-colnames(intPIRres)
    message("Junctions Alt ss PSI completed")
    as@altPSI<-altPSIRes
##########################################################################
    dfES<-rbind(dfEvents[dfEvents$event=="ES",], dfEvents[dfEvents$event=="-",], dfEvents[dfEvents$event=="ES*",] )
    dfESSub<-dfES[,colnames(dfES)!="J"]
    dfESSub<-dfESSub[,colnames(dfESSub)!="event"]
    colnames(dfESSub)  <-cN
##########################################################################
    EsPSI<-function(x){
      x[is.na(x)] <- 0 # we have to remove NAs
      res<-(x[1]+x[2])/(x[1]+x[2]+x[3]*2)
      return(res)
    }
##########################################################################
    colnames(dfESSub)<-paste(ff, rep(1:3,each=length(targets$condition)))
    myfactor<-factor(paste(ff, rep(1:3,each=length(targets$condition))), 
                     levels=paste(rep(pair, each=3 ), 1:3))
    dfESSum<-t(apply(dfESSub, 1, 
                     function(x){tapply(as.numeric(x), 
                                        INDEX=myfactor,
                                        sum  )}))
    colnames(dfESSum)<-rep(unique(targets$condition), each=3)
    myfactor<-factor(rep(unique(targets$condition),each=3), levels=pair)
    EsPSIRes<-cbind(dfES, t(apply(dfESSum, 1, 
                                  function(x){tapply(as.numeric(x),
                                                     INDEX=myfactor, 
                                                     EsPSI  )})))
##########################################################################
    colnames(EsPSIRes)<-colnames(intPIRres)
    message("Junctions ES ss PSI completed")
    as@esPSI<- EsPSIRes 
    #junctions discovery
    as@join<-rbind(altPSIRes,EsPSIRes, intPIRres)
    return(as)
})
##########################################################################
setGeneric (
  name= "writeAS",
  def=function(as, output.dir="as")
  {standardGeneric("writeAS")})
##########################################################################
setMethod(
  f= "writeAS",
  signature= "ASpli-AS",
  definition=function(as, output.dir="as")
  {
    currentDir<-getwd()
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
    file<-paste(exonsDir, "exon.altPSI.tab", sep="/")
    write.table(altPSI(as), file, sep="\t", quote=FALSE, col.names=NA)
    file<-paste(exonsDir, "exon.altES.tab", sep="/")
    write.table(esPSI(as), file, sep="\t", quote=FALSE, col.names=NA)
    ################ INTRONS ##################################
    file<-paste(intronsDir, "intron.irPIR.tab", sep="/")
    write.table(irPIR(as), file, sep="\t", quote=FALSE, col.names=NA)
    ################ Junctions ################################
    file<-paste(junctionsDir, "junction.PIR.tab", sep="/")
    write.table(as@junctionsPIR, file, sep="\t", quote=FALSE, col.names=NA)
    file<-paste(junctionsDir, "junction.PSI.tab", sep="/")
    write.table(junctionsPSI(as), file, sep="\t", quote=FALSE, col.names=NA)
    ################ ALL ######################################
    currentDir=getwd()
    file<-paste(outputDir, "as_discovery.tab", sep="/")
    write.table(joint(as), file, sep="\t", quote=FALSE, col.names=NA)
    }
)
##########################################################################
setGeneric (
  name= "DUreport",
  def=function(counts, targets, pair, group,
               minGenReads=NULL,
               minBinReads=NULL,
               minRds=NULL,
               ignoreExternal=NULL,
               threshold=NULL)
  {standardGeneric("DUreport")})
##########################################################################
setMethod(
  f= "DUreport",
  signature= "ASpli-counts",
  definition=function(counts, targets, pair, group, 
                      minGenReads=NULL,
                      minBinReads=NULL,
                      minRds=NULL,
                      ignoreExternal=NULL,
                      threshold=NULL)
  {
    du<-new(Class="ASpli-DU")
    #define parameters#
    if(is.null(minGenReads)){minGenReads=10}
    if(is.null(minBinReads)){minBinReads=5}
    if(is.null(minRds)){minRds=0.05}
    if(is.null(threshold)){threshold=5}
    ##################################
    df0<-countsg(counts)
    dfG0<-.filterByReads(df0=df0,
                    targets=targets,
                    min=minGenReads,
                    type="any")
    dfGen<-.filterByRdGen(df0=dfG0,
                     targets=targets,
                     min=minRds,
                     type="any")
    genesde<-.genesDE(df=dfGen, 
                 targets=targets, 
                 pair=pair,
                 group=group) 
    message("Genes DE completed")
    du@genes<-genesde
    ##########################################################################
    dfG0<-.filterByReads(df0=df0,
                    targets=targets,
                    min=minGenReads,
                    type="all")
    dfGen<-.filterByRdGen(df0=dfG0,
                     targets=targets,
                     min=minRds,
                     type="all")
    dfBin<-countsb(counts)[countsb(counts)[,"locus"]%in%row.names(dfGen),]
    df1<-.filterByReads(df0=dfBin,
                   targets=targets, 
                   min=minBinReads,
                   type="any")
    df2<-.filterByRdBinRATIO(
    dfBin=df1,
    dfGen=dfGen,
    targets=targets, 
    min=minRds,
    type="any")
#bins con AS en binsN
  if (is.null(ignoreExternal)){ignoreExternal=TRUE}
  binsdu<-.binsDU(df=df2,
                  targets,
                  dfGen,
                  ignoreExternal=ignoreExternal, 
                  pair,
                  group) 
du@bins<-binsdu
message("Bins DU completed")
##########################################################################
df0<-countsj(counts)[countsj(counts)[,"gene"]%in%rownames(dfGen),]
df<-.filterJunctionBySample(df0=df0, 
                             targets=targets, 
                             threshold=threshold)  #mean > one of the condition
df<-df[df$multipleHit=="-",]
junctionsdeSUM<-.junctionsDU_SUM(df,
                                    targets, 
                                    genesde,
                                    pair, 
                                    group,dfGen)
du@junctions<-junctionsdeSUM
message("Junctions DU completed")
return(du)
}
)
##########################################################################
setGeneric (
  name= "writeDU",
  def=function(du, output.dir="du")
  {standardGeneric("writeDU")})
##########################################################################
setMethod(
  f= "writeDU",
  signature= "ASpli-DU",
  definition=function(du, output.dir="du")
  {
    currentDir<-getwd()
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
    file<-paste(genesDir, "gene.de.tab", sep="/")
    write.table(du@genes, file, sep="\t", quote=FALSE, col.names=NA)
    ################ EXONS ###################################
    exonBins<-binsDU(du)[binsDU(du)$feature == "E",]
    exonBins<-exonBins[exonBins$event !="IR",]
    file<-paste(exonsDir, "exon.du.tab", sep="/")
    write.table(exonBins, file, sep="\t", quote=FALSE, col.names=NA)
    ################ INTRONS ##################################
    intronBins<-rbind(binsDU(du)[binsDU(du)$feature == "I",], 
                      binsDU(du)[binsDU(du)$feature == "Io",],
                      binsDU(du)[binsDU(du)$event == "IR",])
    file<-paste(intronsDir, "intron.du.tab", sep="/")
    write.table(intronBins, file, sep="\t", quote=FALSE, col.names=NA)
    ################ Junctions ##################################
    file<-paste(junctionsDir, "junction.du.tab", sep="/")
    write.table(junctionsDU(du), file, sep="\t", quote=FALSE, col.names=NA)
  }
)
##########################################################################
setGeneric (
  name= "writeAll",
  def=function(counts, du, as, output.dir="output")
  {standardGeneric("writeAll")})
##########################################################################
setMethod(
  f= "writeAll",
  definition=function(counts, du, as, output.dir="output")
  {
    writeCounts(counts, output.dir)
    writeRds(counts, output.dir)
    writeAS(as, output.dir)
    writeDU(du, output.dir)
    #armo el df
    currentDir=getwd()
    colnames(as@irPIR)<-colnames(as@altPSI)
    conP<-rbind(altPSI(as), 
                esPSI(as),
                irPIR(as))
    
    ii<-match(rownames(binsDU(du)), row.names(conP))
    bins.join<-data.frame(binsDU(du), conP[ii,])
    bins.join$feature= NULL
    bins.join$event.1= NULL
    summary<-bins.join[,c(1:11) ]
    summary<-cbind(summary,bins.join[,colnames(bins.join)==levels(group)])
    currentDir<-getwd()
    outputDir <- paste(currentDir, output.dir, sep = "/")       
    file<-paste(outputDir, "bins_du_psi_pir.tab", sep="/")
    write.table(bins.join, file, sep="\t", quote=FALSE, col.names=NA)  
    file<-paste(outputDir, "summary.tab", sep="/")
    write.table(summary, file, sep="\t", quote=FALSE, col.names=NA)  
    })
#############Accesors###############################
setGeneric (
  name= "featuresg",
  def=function(x){standardGeneric("featuresg")})
setMethod (f="featuresg",
           signature= "ASpli-features",
           definition=function(x){
             x@genes
           })
#################################################
setGeneric (
  name= "featuresj",
  def=function(x){standardGeneric("featuresj")})
setMethod (f="featuresj",
           signature= "ASpli-features",
           definition=function(x){
             x@junctions
           })
#################################################
setGeneric (
  name= "featuresb",
  def=function(x){standardGeneric("featuresb")})
setMethod (f="featuresb",
           signature= "ASpli-features",
           definition=function(x){
             x@bins
           })
#################################################
#Counts
setGeneric (
  name= "countsg",
  def=function(x){standardGeneric("countsg")})
setMethod (f="countsg",
           signature= "ASpli-counts",
           definition=function(x){
             x@gene.counts
           })
#################################################
setGeneric (
  name= "countsj",
  def=function(x){standardGeneric("countsj")})
setMethod (f="countsj",
           signature= "ASpli-counts",
           definition=function(x){
             x@junction.counts
           })
#################################################
setGeneric (
  name= "countsb",
  def=function(x){standardGeneric("countsb")})
setMethod (f="countsb",
           signature= "ASpli-counts",
           definition=function(x){
             x@exon.intron.counts
           })
#################################################
setGeneric (
  name= "countse1i",
  def=function(x){standardGeneric("countse1i")})
setMethod (f="countse1i",
           signature= "ASpli-counts",
           definition=function(x){
             x@e1i.counts
           })
#################################################
setGeneric (
  name= "countsie2",
  def=function(x){standardGeneric("countsie2")})
setMethod (f="countsie2",
           signature= "ASpli-counts",
           definition=function(x){
             x@ie2.counts
           })
#################################################
setGeneric (
  name= "rdsg",
  def=function(x){standardGeneric("rdsg")})
setMethod (f="rdsg",
           signature= "ASpli-counts",
           definition=function(x){
             x@gene.rd
           })
#################################################
setGeneric (
name= "rdsb",
  def=function(x){standardGeneric("rdsb")})
setMethod (f="rdsb",
           signature= "ASpli-counts",
           definition=function(x){
             x@bin.rd
           })
#################################################
#Accesores AS
setGeneric (
  name= "irPIR",
  def=function(x){standardGeneric("irPIR")})
setMethod (f="irPIR",
           signature= "ASpli-AS",
           definition=function(x){
             x@irPIR
           })
#################################################
setGeneric (
  name= "altPSI",
  def=function(x){standardGeneric("altPSI")})
setMethod (f="altPSI",
           signature= "ASpli-AS",
           definition=function(x){
             x@altPSI
           })
#################################################
setGeneric (
  name= "esPSI",
  def=function(x){standardGeneric("esPSI")})
setMethod (f="esPSI",
           signature= "ASpli-AS",
           definition=function(x){
             x@esPSI
           })
#################################################
setGeneric (
  name= "junctionsPIR",
  def=function(x){standardGeneric("junctionsPIR")})
setMethod (f="junctionsPIR",
           signature= "ASpli-AS",
           definition=function(x){
             x@junctionsPIR
           })
#################################################
setGeneric (
  name= "junctionsPSI",
  def=function(x){standardGeneric("junctionsPSI")})
setMethod (f="junctionsPSI",
           signature= "ASpli-AS",
           definition=function(x){
             x@junctionsPSI
           })
#################################################
setGeneric (
  name= "joint",
  def=function(x){standardGeneric("joint")})
  setMethod (f="joint",
           signature= "ASpli-AS",
           definition=function(x){
             x@join
           })
#################################################
####DU DE
setGeneric (
  name= "genesDE",
  def=function(x){standardGeneric("genesDE")})
setMethod (f="genesDE",
           signature= "ASpli-DU",
           definition=function(x){
             x@genes
           })
#################################################
setGeneric (
  name= "binsDU",
  def=function(x){standardGeneric("binsDU")})
setMethod (f="binsDU",
           signature= "ASpli-DU",
           definition=function(x){
             x@bins
           })
#################################################
setGeneric (
  name= "junctionsDU",
  def=function(x){standardGeneric("junctionsDU")})
setMethod (f="junctionsDU",
           signature= "ASpli-DU",
           definition=function(x){
             x@junctions
           })
#################################################




