library(CNprep)
source("./scripts/helperFunctions.R")

#
# Retrieve information of cytobands (downloaded for genome hg19)
#
retrieveCytobands <- function(dir = "cytoBand.txt"){
  cytobands <- read.table(dir, header=F, sep = "\t", stringsAsFactors = F)
  names(cytobands) <- c("chrom", "start", "end", "cytoloc", "stain")
  return(cytobands)
}

#
# With the chromosome and chromosome location, retrieve the cytoband location
#
findCytolocation <- function(cytobands, chrom, chrom.position){
  chrom = if(chrom == 23) "X" else if (chrom == 24 ) "Y" else chrom
  row <- cytobands[cytobands$chrom == paste("chr", chrom, sep = "") & cytobands$start <= chrom.position & cytobands$end >= chrom.position, ]
  returnme <- data.frame(row)$cytoloc
  return(returnme)
}

#
# Retrieve the norminput argument for CNprep::CNpreprocessing()
#
retrieveNormInput <- function(normalSegments){
  norminput <- data.frame(stringsAsFactors = FALSE)
  for(normalSegments.index in seq(1, nrow(normalSegments))){
    norminput.entry <- data.frame(length = normalSegments[normalSegments.index, ]$end - normalSegments[normalSegments.index, ]$start, segmedian = normalSegments[normalSegments.index, ]$cnlr)
    norminput <- rbind(norminput, norminput.entry)
  }
  return(norminput)
}

#
# Filter norminput from artifacts
#
filterNormInput <- function(norminput, length_threshold=10000000){
  # Determine cutoff
  oNI <- norminput[order(norminput$length), ]
  #plot(oNI$length, oNI$segmedian, ylim = c(-2.25, 2.25), xlim = c(0, 300000000))
  #dev.off()
  nNI <- norminput[norminput$length > length_threshold & abs(norminput$segmedian) < 0.5,]
  #plot(nNI$length, nNI$segmedian, ylim = c(-2.25, 2.25), xlim = c(0, 300000000))
  return(nNI)
}

#
# Retrieve the seginput argument for CNprep::CNpreprocessing()
#
retrieveSegInput <- function(cnvkit_segment_data, segment_probes_assignment, sample, chromosomeSizes, cytobands){
  seginput <- data.frame(stringsAsFactors = FALSE)
  
  # TODO: Determine justification for this - what are these chromosomes?
  accepted_chrom_list <- as.character(chromosomeSizes$chrom)
  remove_extraneous_chrom_indices <- as.numeric(rownames(cnvkit_segment_data[!(cnvkit_segment_data$chromosome %in% accepted_chrom_list), ]))
  cnvkit_segment_data <- cnvkit_segment_data[-remove_extraneous_chrom_indices,]
  segment_probes_assignment <- segment_probes_assignment[-remove_extraneous_chrom_indices,]
  #TODO: Why is segments$probe different than my manual calculation of the probe IDs? 
  num.probes <- segment_probes_assignment$end - segment_probes_assignment$start
  # Iterate through each segment
  for(cnvkit_segment_data.index in seq(1, nrow(cnvkit_segment_data))){
    # Get absolute position of segment
    # TODO: Need to remove extraneous chromosomes beforehand
    raw_chrom <- cnvkit_segment_data[cnvkit_segment_data.index,][[1]]
    substr_chrom <- substr(raw_chrom, 4, nchar(raw_chrom))
    substr_chrom <- if(substr_chrom == "X") 23 else if(substr_chrom == "Y") 24 else substr_chrom
    arg_chrom <- as.numeric(substr_chrom)
    abs_position <- chromsomeToAbsoluteBPConversionForSingleEntry(arg_chrom, cnvkit_segment_data[cnvkit_segment_data.index,][[2]], cnvkit_segment_data[cnvkit_segment_data.index,][[3]], chromosomeSizes)
    
    
    
    cytoband.my.start <- findCytolocation(cytobands = cytobands, chrom = arg_chrom, chrom.position = cnvkit_segment_data[cnvkit_segment_data.index,][[2]])
    cytoband.my.end <- findCytolocation(cytobands = cytobands, chrom = arg_chrom, chrom.position = cnvkit_segment_data[cnvkit_segment_data.index,][[3]])
    
    #
    # TODO: cytoband.my.end sometimes NA - temporarily removing. Check hT1 row 22 end = 198153431
    #
    seginput.entry <- data.frame(ID = sample, start = segment_probes_assignment[cnvkit_segment_data.index, ]$start, end = segment_probes_assignment[cnvkit_segment_data.index, ]$end, 
                                 num.probes = num.probes[cnvkit_segment_data.index], seg.median = cnvkit_segment_data[cnvkit_segment_data.index, ]$log2, 
                                 chrom = arg_chrom, chrom.pos.start = cnvkit_segment_data[cnvkit_segment_data.index, ][[2]], 
                                 chrom.pos.end = cnvkit_segment_data[cnvkit_segment_data.index, ][[3]], abs.pos.start = abs_position$start,
                                 abs.pos.end = abs_position$end)
    
    
    seginput <- rbind(seginput, seginput.entry)
  }
  return(seginput)
}

#
# Retrieve the ratinput argument for CNprep::CNpreprocessing()
#
retrieveRatInput <- function(cnvkit_bins_data, sample){
  ratinput <- data.frame(cnvkit_bins_data$log2)
  names(ratinput) <- c(sample)
  return(ratinput)
}

#
# Wrapper function for CNprep::CNpreprocessing()
#
runCNpreprocessing <- function(seginput, ratinput, norminput,
                               blsize=50, minjoin=0.25, cweight=0.4, bstimes=50,
                               chromrange=1:22, distrib="Rparallel", njobs=4, modelNames="E", ntrial= 10){
  segtable<-CNpreprocessing(segall=seginput,ratall=ratinput,"ID","start","end",
                            chromcol="chrom",bpstartcol="chrom.pos.start",bpendcol="chrom.pos.end",blsize=blsize,
                            minjoin=minjoin,cweight=cweight,bstimes=bstimes,chromrange=chromrange,distrib=distrib,njobs=njobs,
                            modelNames=modelNames,normalength=norminput[,1],normalmedian=norminput[,2], ntrial=ntrial)
  return(segtable)
}
