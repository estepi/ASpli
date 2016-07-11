#############Accesors###############################
setGeneric (
  name = "featuresg",
  def = function(x)standardGeneric("featuresg"))
#####################################################
setMethod (f = "featuresg",
           signature = "ASpliFeatures",
           definition = function(x){  x@genes })
#################################################
setGeneric (
  name = "featuresj",
  def = function(x){standardGeneric("featuresj")})
setMethod (f = "featuresj",
           signature = "ASpliFeatures",
           definition = function(x){
             x@junctions
           })
#################################################
setGeneric (
  name = "featuresb",
  def = function(x)
    standardGeneric("featuresb"))
#################################################
setMethod (f = "featuresb",
           signature = "ASpliFeatures",
           definition = function(x){x@bins })
#################################################
#Counts
setGeneric (
  name= "countsg",
  def=function(x) standardGeneric("countsg"))
#################################################
setMethod (f = "countsg",
           signature = "ASpliCounts",
           definition = function(x){ x@gene.counts})
#################################################
setGeneric (
  name = "countsj",
  def = function(x)
    standardGeneric("countsj"))
#################################################
setMethod (f = "countsj",
           signature = "ASpliCounts",
           definition = function(x){x@junction.counts })
#################################################
setGeneric (
  name = "countsb",
  def = function(x)
    standardGeneric("countsb"))
#################################################
setMethod (f = "countsb",
           signature = "ASpliCounts",
           definition = function(x){ x@exon.intron.counts })
#################################################
setGeneric (
  name = "countse1i",
  def = function(x) 
    standardGeneric("countse1i"))
#################################################
setMethod (f = "countse1i",
           signature = "ASpliCounts",
           definition =function(x){ x@e1i.counts })
#################################################
setGeneric (
  name = "countsie2",
  def = function(x)
    standardGeneric("countsie2"))
#################################################
setMethod (f = "countsie2",
           signature = "ASpliCounts",
           definition = function(x){ x@ie2.counts })
#################################################
setGeneric (
  name = "rdsg",
  def = function(x)
    standardGeneric("rdsg"))
#################################################
setMethod (f = "rdsg",
           signature = "ASpliCounts",
           definition = function(x){ x@gene.rd })
#################################################
setGeneric (
  name = "rdsb",
  def = function(x)
    standardGeneric("rdsb"))
#################################################
setMethod(f = "rdsb",
          signature = "ASpliCounts",
          definition = function(x){ x@bin.rd })
#################################################
#Accesores AS
setGeneric (
  name =  "irPIR",
  def = function(x)
    standardGeneric("irPIR"))
#################################################
setMethod (f = "irPIR",
           signature =  "ASpliAS",
           definition = function(x){ x@irPIR })
#################################################
setGeneric (
  name =  "altPSI",
  def = function(x)
    standardGeneric("altPSI"))
#################################################
setMethod (f = "altPSI",
           signature = "ASpliAS",
           definition = function(x) {x@altPSI })
#################################################
setGeneric (
  name = "esPSI",
  def = function(x)
    standardGeneric("esPSI"))
#################################################
setMethod (f = "esPSI",
           signature = "ASpliAS",
           definition = function(x) {x@esPSI })
#################################################
setGeneric (
  name = "junctionsPIR",
  def = function(x)
    standardGeneric("junctionsPIR"))
#################################################
setMethod (f ="junctionsPIR",
           signature = "ASpliAS",
           definition = function(x){ x@junctionsPIR })
#################################################
setGeneric (
  name = "junctionsPSI",
  def = function(x)
    standardGeneric("junctionsPSI"))
#################################################
setMethod (f = "junctionsPSI",
           signature = "ASpliAS",
           definition = function(x){x@junctionsPSI })
#################################################
setGeneric (
  name = "joint",
  def = function(x)
    standardGeneric("joint"))
#################################################
setMethod (f = "joint",
           signature = "ASpliAS",
           definition = function(x){ x@join })
#################################################
####DU DE
setGeneric (
  name = "genesDE",
  def = function(x)
    standardGeneric("genesDE"))
#################################################
setMethod(f = "genesDE",
          signature = "ASpliDU",
          definition = function(x){ x@genes })
#################################################
setGeneric(
  name = "binsDU",
  def = function(x)
    standardGeneric("binsDU"))
#################################################
setMethod (f = "binsDU",
           signature = "ASpliDU",
           definition = function(x){ x@bins })
#################################################
setGeneric (
  name = "junctionsDU",
  def = function(x)
    standardGeneric("junctionsDU"))
#################################################
setMethod (f = "junctionsDU",
           signature = "ASpliDU",
           definition = function(x){ x@junctions})
#################################################
#show
setMethod('show', 'ASpliFeatures', 
          function(object)
          {
            cat("Object of class", class(object),"\n")
            cat("Genes: GRangesList of length ", 
                length(object@genes),
               "Access using featuresg(object)", "\n")
            cat("Bins: GRanges of length", 
                length(object@bins),
                "Access using featuresb(object)", "\n")
            cat("Junctions: GRanges of length", 
                length(object@junctions),
                "Access using featuresj(object)")
          })
#################################################
setMethod('show', 'ASpliCounts', 
          function(object)
          {
            cat("Object of class", class(object),"\n")
            cat("Gene counts:", 
                dim(object@gene.counts)[1], "genes analysed.",
                "Access using countsg(object)", "\n")
            cat("Gene RD:", 
                dim(object@gene.rd)[1], "genes analysed.",
                "Access using rdsg(object)", "\n")
            cat("Bin counts:", 
                dim(object@exon.intron.counts)[1], "bins analysed.",
                "Access using countsb(object)", "\n")
            cat("Bin RD:", 
                dim(object@bin.rd)[1],"bins analysed.",
                "Access using rdsb(object)", "\n")
            cat("Junction counts:", 
                dim(object@junction.counts)[1], "junctions analysed.",
              "Access using countsj(object)", "\n")
            })
################################################################
setMethod('show', 'ASpliAS', 
          function(object)
          {
            cat("Object of class", class(object),"\n")
            cat("IR PIR: ", 
                dim(object@irPIR)[1], "intron bins analysed.",
                "Access using irPIR(object)", "\n")
            cat("ES PSI:", 
                dim(object@esPSI)[1], "exon bins analysed.",
                " Access using esPSI(object)", "\n")
            cat("AltSS PSI:", 
                dim(object@altPSI)[1], "exon bins analysed.",
                " Access using altPSI(object)", "\n")
            cat("Junctions PIR:", 
                dim(object@junctionsPIR)[1], "junctions analysed.",
                "Access using junctionsPIR(object)", "\n")
            cat("Junctions PSI:", 
                dim(object@junctionsPSI)[1], "junctions analysed.",
                "Access using junctionsPSI(object)")
            })
################################################################
setMethod('show', 'ASpliDU', 
          function(object)
          {
            cat("Object of class", class(object),"\n")
            cat("Gene DE:", 
                dim(object@genes)[1], "genes analysed.",
                "Access using genesDE(object)", "\n")
            
            cat("Bins DU:", 
                dim(object@bins)[1], "bins analysed.",
                "Access using binsDU(object)", "\n")
            
            cat("Junctions DU:", 
                dim(object@junctions)[1],"junctions analysed.",
                "Access using junctionsDU(object)")
          })
################################################################
#Write methods
setGeneric (
  name =  "writeCounts",
  def = function(counts, output.dir="counts")
    standardGeneric("writeCounts"))
##########################################################################
setMethod(
  f = "writeCounts",
  signature = "ASpliCounts",
  definition =function(counts, output.dir="counts")
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
    
    file <- paste(genesDir, "gene.counts.tab", sep="/")
    write.table(countsg(counts), file, sep="\t", quote=FALSE, col.names=NA)
    ##########################################################################
    ec <- countsb(counts)[countsb(counts)$feature=="E",]
    ec <- ec[ec$event != "IR",]
    ec <- ec[ec$event != "IR*",]
    file <- paste(exonsDir, "exon.counts.tab", sep="/")
    write.table(ec, file, sep="\t", quote=FALSE, col.names=NA)
    ################ INTRONS ##################################
    file <- paste(intronsDir, "intron.counts.tab", sep="/")
    ic <- rbind(countsb(counts)[countsb(counts)$feature == "I",], 
                countsb(counts)[countsb(counts)$feature == "Io",], 
                countsb(counts)[countsb(counts)$event == "IR",],
                countsb(counts)[countsb(counts)$event == "IR*",])
    write.table(ic, file, sep="\t", quote=FALSE, col.names=NA)
    file <- paste(intronsDir, "e1i.counts.tab", sep="/")
    write.table(countse1i(counts), file, sep="\t", quote=FALSE, col.names=NA)
    file <- paste(intronsDir, "ie2.counts.tab", sep="/")
    write.table(countsie2(counts), file, sep="\t", quote=FALSE, col.names=NA)
    file <- paste(junctionsDir, "junction.counts.tab", sep="/" )
    write.table(countsj(counts), file, sep="\t", quote=FALSE,  col.names=NA)
  }
)
##########################################################################
setGeneric (
  name = "writeRds",
  def = function(counts, output.dir="rds")
    standardGeneric("writeRds"))
##########################################################################
setMethod(
  f = "writeRds",
  signature = "ASpliCounts",
  definition =function(counts, output.dir="rds")
  {
    currentDir <- getwd()
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
    file <- paste(genesDir, "gene.rd.tab", sep="/")
    write.table(rdsg(counts), file, sep="\t", quote=FALSE, col.names=NA)
    ################ EXONS ####################################
    erd <- rdsb(counts)[rdsb(counts)$feature == "E",]
    erd <- erd[erd$event != "IR",]
    erd <- erd[erd$event != "IR*",]
    file <- paste(exonsDir, "exon.rd.tab", sep="/")
    write.table(erd, file, sep="\t", quote=FALSE, col.names=NA)
    ################ INTRONS ##################################
    file <- paste(intronsDir, "intron.rd.tab", sep="/")
    ird <- rbind(rdsb(counts)[rdsb(counts)$feature == "I",], 
                 rdsb(counts)[rdsb(counts)$feature == "Io",], 
                 rdsb(counts)[rdsb(counts)$eventJ == "IR",])
    write.table(ird, file, sep="\t", quote=FALSE, col.names=NA)
  }
)
##########################################################################
setGeneric (
  name = "writeAll",
  def = function(counts, du, as, output.dir="output")
    standardGeneric("writeAll"))
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
    currentDir <- getwd()
    colnames(as@irPIR) <- colnames(as@altPSI)
    conP <- rbind(altPSI(as), 
                  esPSI(as),
                  irPIR(as))
    ii <- match(rownames(binsDU(du)), row.names(conP))
    bins.join <- data.frame(binsDU(du), conP[ii,])
    bins.join$feature <- NULL
    bins.join$event.1 <- NULL
    summary <- bins.join[,c(1:11) ]
    summary <- cbind(summary,bins.join[,colnames(bins.join)==levels(group)])
    currentDir <- getwd()
    outputDir <- paste(currentDir, output.dir, sep = "/")       
    file <- paste(outputDir, "bins_du_psi_pir.tab", sep="/")
    write.table(bins.join, file, sep="\t", quote=FALSE, col.names=NA)  
    file<-paste(outputDir, "summary.tab", sep="/")
    write.table(summary, file, sep="\t", quote=FALSE, col.names=NA)  
  })
