.createGRangesGenes <-
  function(genome, md)
  {
    exons <- exonsBy(genome, by="gen") #extract exons by gene
    exons.by.gene.disjoint <- disjoin(exons) 
    geneChr <- sapply(seqnames(exons.by.gene.disjoint),unique)
    geneChr <- as.character(geneChr)
    geneStarts <- sapply(start(exons.by.gene.disjoint),min)
    geneEnds <- sapply(end(exons.by.gene.disjoint),max)
    browser <- paste(geneChr,geneStarts, sep=":") 
    #coordinates for genome browser
    gene_coordinates <- paste(browser,geneEnds, sep="-")
    #agrego el shareLocus
    if (packageVersion("IRanges")<2.6)
       {
        locus <- findOverlaps(exons.by.gene.disjoint, ignoreSelf=TRUE) 
        #locus by locus
        }
    else {
        locus <- findOverlaps(exons.by.gene.disjoint, drop.self=TRUE) 
        #locus by locus   
        }  
    locus.df <- as.data.frame(locus)
    locus_overlap <- rep("-", length(exons.by.gene.disjoint))  
    #add a conditional 
    if (nrow(locus.df)>0) {
      locus.df$names <- names(exons.by.gene.disjoint[locus.df$subjectHits])
      tt1 <- data.frame(aggregate(names ~ queryHits, 
                                  data = locus.df, 
                                  paste, collapse = ';'))
      locus_overlap[tt1$queryHits] <- tt1$names
    }
    ##############
        mcols(exons.by.gene.disjoint) <- append(
            mcols(exons.by.gene.disjoint), 
            DataFrame(gene_coordinates = gene_coordinates))
    mcols(exons.by.gene.disjoint) <- append(
            mcols(exons.by.gene.disjoint), 
            DataFrame(locus_overlap=locus_overlap))
    symbol <- md[names(exons.by.gene.disjoint),]
    mcols(exons.by.gene.disjoint) <- append(mcols(exons.by.gene.disjoint), 
                                     DataFrame(symbol=symbol))
    return(exons.by.gene.disjoint)
  }
#######################################################
.createGRangesExons <-
  function(genome, md) {
    exons <- exonsBy(genome, by="gen") #extract exons by gen 
    exons.partitioning <- as.data.frame(exons@partitioning)
    genes.multiexon.index <- exons.partitioning$width>1 
    names(genes.multiexon.index) <- exons.partitioning$names
    exonsM <- exons[genes.multiexon.index] #this exons are extracted by gen 
    exon.by.gene.disjoint <- disjoin(exonsM)# gene name is kept
    exon.by.gene.disjoint.unlist <- unlist(exon.by.gene.disjoint ) 
    exon.bins <- exon.by.gene.disjoint.unlist[!duplicated(exon.by.gene.disjoint.unlist)] 
    #exon.gene.names <- table(exon.bins@ranges@NAMES)
    #change
    exon.gene.names<-table(factor(exon.bins@ranges@NAMES, 
                                  levels=unique(exon.bins@ranges@NAMES)))
########################################################
    exon.bins.num <- lapply(exon.gene.names,function(x){seq(1:x)})
    exon.bins.num <- unlist(exon.bins.num)
    exon.bins.id <- sprintf('E%03d', exon.bins.num)
    feature.exon <- rep("E", length(exon.bins.id))
    exon.gene.names.repeated <- rep(names(exon.gene.names), exon.gene.names) 
    exon.gene.names.repeated <- rep(names(exon.gene.names), exon.gene.names) 
    exon.gene.names.repeated <- as.character(exon.gene.names.repeated)
    #add metadata
    #add locus
    mcols(exon.bins) <- append(mcols(exon.bins), 
                               DataFrame(locus=exon.gene.names.repeated) )
    #add gene coord
    mcols(exon.bins) <- append(mcols(exon.bins), DataFrame(bin=exon.bins.id))
    mcols(exon.bins) <- append(mcols(exon.bins), DataFrame(feature=feature.exon))
    exon.bins.names <- paste(exon.gene.names.repeated, exon.bins.id ,sep=":")
    exon.bins@ranges@NAMES <- exon.bins.names
    ###add symbol ####
    tt <- match(exon.gene.names.repeated, rownames(md))  
    symbol <- md[tt, ]
    mcols(exon.bins) <- append(mcols(exon.bins), 
                               DataFrame(symbol=symbol))
    return (exon.bins)
  }
########################################################
.createGRangesIntrons <-
  function(genome, md)
  {
    introns <- intronsByTranscript(genome, use.names=TRUE) 
    #extract introns by transcript, from genome
    introns.ulst <- unlist(introns)     #unlist element
    introns.no.dups <- introns.ulst[!duplicated(introns.ulst)] #remove dulicates 
    introns.tx.names <- names(introns.no.dups) 
    #extract names from transcripts keep number of transcript
    introns.map <- select(genome, keys=introns.tx.names, keytype='TXNAME', columns='GENEID') 
    #select genes names. It is a data.frame with introns.map$TXNAME and introns.map$GENEID  cols
    introns.idx <- introns.map$GENEID[!is.na(introns.map$GENEID)] #character vector  
    introns.by.gene <- split(introns.no.dups[!is.na(introns.map$GENEID)], introns.idx) 
    #split original granges in introns by gene and not by tx as original 
    #-> it should be the same number as multiexonic genes.
    introns.by.gene@unlistData@ranges@NAMES <- NULL #delete old names
    intron.o <- unlist(introns.by.gene) #disjoint in bins
    
    #intron.o.gene.names <- table(intron.o@ranges@NAMES) #get vector with gene names
    #change
    intron.o.gene.names<-table(factor(intron.o@ranges@NAMES,
                                      levels=unique(intron.o@ranges@NAMES)))    
    
    intron.o.gene.names.repeated <- rep(names(intron.o.gene.names), intron.o.gene.names) 
    #repeat character vector by numeric vector 
    intron.o.gene.names.repeated <- as.character(intron.o.gene.names.repeated) #as character
    mcols(intron.o) <- append(mcols(intron.o), DataFrame(locus=intron.o.gene.names.repeated)) 
    #add metadata
    intron.o.num <- lapply(intron.o.gene.names,function(x){seq(1:x)})
    intron.o.num <- unlist(intron.o.num)
    intron.o.num <- sprintf('Io%03d', intron.o.num)
    feature.o.intron <- rep("Io", length(intron.o.num))
    mcols(intron.o) <- append(mcols(intron.o), DataFrame(bin=intron.o.num))  
    #add metadata
    mcols(intron.o) <- append(mcols(intron.o), DataFrame(feature=feature.o.intron)) 
    #add metadata
    intron.o.names<-paste(intron.o.gene.names.repeated, intron.o.num,sep=":") 
    #combine gene names and bin number 
    intron.o@ranges@NAMES <- intron.o.names
    ##########################################################################
    introns.by.gene.disjoint <- disjoin(introns.by.gene) 
    #disjoin an get INTRONS BINS 
    introns.by.gene@unlistData@ranges@NAMES <- NULL 
    #delete old names
    intron.bins <- unlist(introns.by.gene.disjoint) 
    #disjoint in bins
    intron.gene.names <- table(intron.bins@ranges@NAMES) 
    #get vector with gene names
    intron.gene.names.repeated <- rep(names(intron.gene.names), intron.gene.names) 
    #repeat character vector by numeric vector 
    intron.gene.names.repeated <- as.character(intron.gene.names.repeated) #as character
    mcols(intron.bins) <- append(mcols(intron.bins), 
                                 DataFrame(locus=intron.gene.names.repeated)) 
    #add metadata
    intron.bins.num <- lapply(intron.gene.names,function(x){seq(1:x)})
    intron.bins.num <- unlist(intron.bins.num)
    intron.bins.num <- sprintf('I%03d',  intron.bins.num)
    feature.intron <- rep("I", length(intron.bins.num) )
    mcols(intron.bins) <- append(mcols(intron.bins), DataFrame(bin=intron.bins.num))  #add metadata
    mcols(intron.bins) <- append(mcols(intron.bins), DataFrame(feature=feature.intron))  #add metadata
    intron.bin.names <- paste(intron.gene.names.repeated, intron.bins.num,sep=":") #combine gene names and bin number 
    intron.bins@ranges@NAMES <- intron.bin.names
    ################################################
    intrones.totales <- c(intron.bins, intron.o)
    intron.tot.u <- sort(unique(intrones.totales))
    ###add symbol ####
    tt <- match(intron.tot.u@elementMetadata$locus, rownames(md))#sort
    symbol <- md[tt, ]
    mcols(intron.tot.u) <-append(mcols(intron.tot.u), DataFrame(symbol=symbol))
    return (intron.tot.u)
  }
########################################################
.createGRangesTranscripts <-
  function(genome) 
  {
    transcripts <- transcriptsBy(genome) #extract transcripts coordinates
    
    return(transcripts)
  }
########################################################
  .createGRangesJunctions <-
    function(genome) {
      introns <- intronsByTranscript(genome, use.names=TRUE) 
      introns.ulst <- unlist(introns)          
      introns.no.dups <- introns.ulst[!duplicated(introns.ulst)] 
      introns.tx.names <- names(introns.no.dups) 
      introns.map <- select(genome, keys=introns.tx.names, keytype='TXNAME', columns='GENEID') 
      introns.idx <- introns.map$GENEID[!is.na(introns.map$GENEID)] #character vector  
      introns.by.gene <- split(introns.no.dups[!is.na(introns.map$GENEID)], introns.idx)
      introns.by.gene@unlistData@ranges@NAMES <- NULL #delete old names
      introns.ul <- unlist(introns.by.gene) 
      introns.gene.names <- table(introns.ul@ranges@NAMES) #get vector with gene names
      introns.gene.names.repeated <- rep(names(introns.gene.names), introns.gene.names) #repeat char
      introns.gene.names.repeated <- as.character(introns.gene.names.repeated) #as character
      mcols(introns.ul) <- append(mcols(introns.ul), 
                                  DataFrame(locus=introns.gene.names.repeated)) #add metadata
      introns.num <- lapply(introns.gene.names,function(x){seq(1:x)})
      introns.num <- unlist(introns.num)
      introns.num <- sprintf('J%03d',  introns.num)
      junction.names <- paste(introns.gene.names.repeated, introns.num,sep=":") 
      #combine gene names and bin number 
      junctions <- introns.ul
      start(junctions) <- start(introns.ul)-1
      end(junctions) <- end(introns.ul)+1
      junctions@ranges@NAMES <- junction.names
      junctions <- sort(junctions)
      return (junctions)
    }
#########################################################  
  .findASBin <-
    function(exon.bins, intron.bins, transcripts, junctions)
    {
      
      as.bins <- findOverlaps(exon.bins, intron.bins, type=c("equal")) 
      exon.as <- rep("-", as.bins@nLnode) #vector for exons metadata
      intron.as <- rep("-", as.bins@nRnode) #vector for introns metadata
      #Asign AS first
      exon.as[as.bins@from] <- "as" #introns exons as;
      intron.as[as.bins@to] <- "as" #identify introns as
      transcripts.unlist  <- unlist(transcripts) #41671  
      find.external.start <- findOverlaps(exon.bins, transcripts.unlist, type=c("start")) #
      find.external.end   <- findOverlaps(exon.bins, transcripts.unlist, type=c("end"))
      exon.as[find.external.start@from] <- "external" #identify start exons bins
      exon.as[find.external.end@from] <- "external" #identify end exons bins
      #Asign external second
      #add metadata 
      #if a bin is AS and external, AS tag will be replaced by external.
      mcols(intron.bins) <- append(mcols(intron.bins), DataFrame(class=intron.as)) 
      ###add introns.bins metadata
      mcols(exon.bins) <- append(mcols(exon.bins), DataFrame(class=exon.as)) 
      ###add introns.bins metadata
      exons.introns <- append(exon.bins, intron.bins) 
      exons.introns.sort <- sort(exons.introns) ##sort and merge
      #deberiamos no dejar los unique sino todos!
      exons.introns.unique <- unique(exons.introns.sort) ##delete duplicates
      #build dataframe for the loop
      auxdf <- data.frame(feature=as.character(exons.introns.unique@elementMetadata$feature), 
                        class=as.character(exons.introns.unique@elementMetadata$class),
                        strand=as.character(exons.introns.unique@strand))
    ###########classifying AS bins ################################
      #events <- rep("-",nrow(auxdf))
      events<-as.character(auxdf$class)
      #eventsJ <- rep("-",nrow(auxdf))
      eventsJ<-as.character(auxdf$class)
      #counters, just for   check
      IR <- 0
      ES <- 0
      Alt5ss <- 0 
      Alt3ss <- 0 
      totAS <- 0
      mult <- 0
      multPos <- c()
      multIR <- 0
      multES <- 0
      multAlt5ss <- 0 
      multAlt3ss <- 0 
      AsNotExternal <- sum(intron.as=="as") 
      for (i in 1:nrow(auxdf)){  
        if  (auxdf$class[i] =="as")  {
          totAS = totAS +1  #OK
          if ( auxdf$class[i+1] == "as"| auxdf$class[i-1] == "as") 
          {#undef
            mult=mult + 1; multPos=append(multPos, c(i, i+1))
            #primero confirmo si no existe una juntura igual, 
            # todos los as bins son exones por default 
            ji <- junctions
            start(ji) <- start(junctions)+1
            end(ji) <- end(junctions)-1
            if (length(findOverlaps(exons.introns.unique[i,], ji, type="equal"))!= 0 ){ 
              events[i] <- "IR*"
              eventsJ[i] <- "IR"
              multIR = multIR+1
            } else  {
              jranges.start <- junctions
              start(jranges.start) <- end(jranges.start)
              st <- length(findOverlaps(type = "start", exons.introns.unique[i,], jranges.start))
              jranges.end <- junctions
              end(jranges.end) <- start(jranges.end)
              end <- length(findOverlaps(type="end", exons.introns.unique[i,], jranges.end))
              
              if ( (st > 0) && (end > 0) ) { multES = multES+1;  
                                             events[i] <- "ES*"; 
                                             eventsJ[i] <- "ES"}
              
              else if ((st > 0) && (auxdf$strand[i] =="+")) { multAlt3ss = multAlt3ss +1;
                                                               events[i] <- "Alt3ss*";  
                                                               eventsJ[i]<-"Alt3ss" }     
              else if ((st > 0) && (auxdf$strand[i]=="-" )) { multAlt3ss = multAlt3ss +1; 
                                                               events[i] <- "Alt3ss*";  
                                                               eventsJ[i]<-"Alt3ss"} 
              else if ((end > 0) & (auxdf$strand[i]=="+" )) { multAlt5ss = multAlt5ss + 1; 
                                                               events[i] <- "Alt5ss*"; 
                                                               eventsJ[i]<-"Alt5ss"}
              else if ((end > 0) && (auxdf$strand[i]=="-")) { multAlt5ss = multAlt5ss + 1; 
                                                              events[i]  <-  "Alt5ss*"; 
                                                              eventsJ[i] <- "Alt5ss"}           
            } 
          }
          #not undef
          else{      
            if ( (auxdf$feature[i-1] == "E") &&  (auxdf$feature[i+1] == "E") ) { 
              IR <- IR + 1; 
              events[i] <- "IR"; 
              eventsJ[i] <- "IR" }
            
            else if ((auxdf$feature[i-1] == "I") && (auxdf$feature[i+1] == "I") ) { 
              ES <- ES + 1; 
              events[i]  <- "ES"; 
              eventsJ[i] <- "ES"}
            else if  ((auxdf$feature[i-1] == "E") &&  (auxdf$feature[i+1] == "I") 
                      &&  (auxdf$strand[i] == "+"))  {
              Alt5ss <- Alt5ss + 1; 
              events[i] <- "Alt5ss"; 
              eventsJ[i] <- "Alt5ss"}
            else if  ((auxdf$feature[i-1] == "E") &&  (auxdf$feature[i+1] == "I") 
                      &&  (auxdf$strand[i] == "-"))  {
              Alt3ss <-  Alt3ss + 1; 
              events[i] <- "Alt3ss"; 
              eventsJ[i] <- "Alt3ss"}
            else if	((auxdf$feature[i-1] == "I") &&  (auxdf$feature[i+1] == "E") 
                     &&  (auxdf$strand[i] == "+"))  {
              Alt3ss <- Alt3ss + 1; 
              events[i] <- "Alt3ss"; 
              eventsJ[i] <- "Alt3ss"}
            else if	((auxdf$feature[i-1] == "I") &&  (auxdf$feature[i+1] == "E") 
                     &&  (auxdf$strand[i] == "-"))  {
              Alt5ss <-  Alt5ss + 1; 
              events[i] <- "Alt5ss"; 
              eventsJ[i] <- "Alt5ss"}
          }
        }
      }
      message("* Number of AS bins (not include external) =", totAS)
      message("* Number of AS bins (include external) =", AsNotExternal)
      message("* Classified as:", "\n", 
              "\t", "ES bins = ", ES, "\t","(",round(ES/totAS*100), "%)" , "\n", 
              "\t","IR bins = " , IR,"\t","(",round(IR/totAS*100), "%)" , "\n",
              "\t","Alt5'ss bins = ", Alt5ss, "\t","(",round(Alt5ss/totAS*100), "%)" , "\n", 
              "\t","Alt3'ss bins = ", Alt3ss,"\t","(",round(Alt3ss/totAS*100), "%)" , "\n",
              "\t", "multiple AS bins = ", mult, "\t","(",round(mult/totAS*100), "%)" , "\n", 
              "classified as:", "\n",
              "\t\t\t", "ES bins = ", multES, "\t","(",round(multES/mult*100), "%)" , "\n", 
              "\t\t\t","IR bins = " , multIR,"\t","(",round(multIR/mult*100), "%)" , "\n",
              "\t\t\t","Alt5'ss bins = ", multAlt5ss, "\t","(",round(multAlt5ss/mult*100), "%)" , "\n", 
              "\t\t\t","Alt3'ss bins = ", multAlt3ss,"\t","(",round(multAlt3ss/mult*100), "%)" , "\n")
      ###########################################################################################
      mcols(exons.introns.unique) <- append(mcols(exons.introns.unique), DataFrame(event=events))
      mcols(exons.introns.unique) <- append(mcols(exons.introns.unique), DataFrame(eventJ=eventsJ))
      ###########################################################################
      cat("* Number of AS bins (not include external) =", totAS,  "\n" )
      cat("* Number of AS bins (include external) =", AsNotExternal,  "\n")
      cat("* Classified as:", "\n", 
              "\t", "ES bins = ", ES, "\t","(",round(ES/totAS*100), "%)" , "\n", 
              "\t","IR bins = " , IR,"\t","(",round(IR/totAS*100), "%)" , "\n",
              "\t","Alt5'ss bins = ", Alt5ss, "\t","(",round(Alt5ss/totAS*100), "%)" , "\n", 
              "\t","Alt3'ss bins = ", Alt3ss,"\t","(",round(Alt3ss/totAS*100), "%)" , "\n",
              "\t", "multiple AS bins = ", mult, "\t","(",round(mult/totAS*100), "%)" , "\n", 
              "classified as:", "\n",
              "\t\t\t", "ES bins = ", multES, "\t","(",round(multES/mult*100), "%)" , "\n", 
              "\t\t\t","IR bins = " , multIR,"\t","(",round(multIR/mult*100), "%)" , "\n",
              "\t\t\t","Alt5'ss bins = ", multAlt5ss, "\t","(",round(multAlt5ss/mult*100), "%)" , "\n", 
              "\t\t\t","Alt3'ss bins = ", multAlt3ss,"\t","(",round(multAlt3ss/mult*100), "%)" , "\n")
      return(exons.introns.unique)
    }
  