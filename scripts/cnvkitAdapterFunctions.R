retrieveCNVkitSegments <- function(sample, dir = "resources/testCNVkit/", genes = TRUE){
  print(paste0("Retrieving file: ", paste0(dir, sample, ".ucsc.hg38.bwa.BQSR.cns")))
  cnvkit_segment_data <- as.data.frame(read.table(paste0(dir, sample, ".ucsc.hg38.bwa.BQSR.cns"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  if(genes == FALSE){
    cnvkit_segment_data <- cnvkit_segment_data[,-c(4)]
  }
  return(cnvkit_segment_data)
}

retrieveCNVkitBins <- function(sample, dir = "resources/testCNVkit/"){
  cnvkit_bins_data <- as.data.frame(read.table(paste0(dir, sample, ".ucsc.hg38.bwa.BQSR.cnr"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  return(cnvkit_bins_data)
}

retrieveCNVkitTargetCoverage <- function(sample, dir = "resources/testCNVkit/"){
  cnvkit_bins_data <- as.data.frame(read.table(paste0(dir, sample, ".ucsc.hg38.bwa.BQSR.targetcoverage.cnn"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  return(cnvkit_bins_data)
}

retrieveCNVkitAntiTargetCoverage <- function(sample, dir = "resources/testCNVkit/"){
  cnvkit_bins_data <- as.data.frame(read.table(paste0(dir, sample, ".ucsc.hg38.bwa.BQSR.antitargetcoverage.cnn"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  return(cnvkit_bins_data)
}

cnvKitSegmentsToBedFormat <- function(segments){
  normalSegments <- normalSegments[,c("chromosome", "start", "end", "log2")]
  names(normalSegments) <- c("chrom", "start", "end", "cnlr")
  return(normalSegments)
}

assignProbeNumbers <- function(segments, bins){
  probeNumbers <- do.call(rbind, lapply(seq(nrow(segments)), function(segment.i, bins){
    probeEntry <- data.frame(start = as.numeric(rownames(bins[bins$chromosome == segments[segment.i, ]$chromosome & bins$start == segments[segment.i, ]$start, ]))[1],
                             end = as.numeric(rownames(bins[bins$chromosome == segments[segment.i, ]$chromosome & bins$end == segments[segment.i, ]$end, ]))[1])
    return(probeEntry)
  }, bins))
  return(probeNumbers)
}

sourceTestCode <- function(){
  setwd("~/Git-Projects/Git-Research-Projects/CNVKit-Workflow-Adapter")
  segments_table <- retrieveCNVkitSegments("hF2", dir = "resources/testCNVkit/", genes = FALSE)
  bins_table <- retrieveCNVkitBins("hF2", dir = "resources/testCNVkit/")
  scheme_table <- retrieveCNVkitScheme("hF2", dir = "resources/testCNVkit/") 
  scheme_table <- retrieveCNVkitAntiTargetCoverage("hF2", dir = "resources/testCNVkit/") 
}
