args <- commandArgs(trailingOnly = TRUE)

#
# Set script arguments
#
output_dir <- "" 
mclust_model <- NULL
minjoin <- NULL
ntrial <- NULL
description <- "No description provided by user."
if (length(args) >= 1){
	output_dir <- args[1]
}
if (length(args) >= 2){
	mclust_model <- args[2]	
}
if (length(args) >= 3){
	minjoin <- as.numeric(args[3])
}
if(length(args) >= 4){
	ntrial <- as.numeric(args[4])
}
if (length(args) >= 5){
	description <- args[5]
}

print(paste("PARAMETERS:", "output_dir=", output_dir, "mclust_model=", mclust_model, "minjoin=", minjoin, "ntrial=", ntrial, "description=", description))
#
# Load source libraries
# TODO: Organize dependencies
#
setwd("~/code/cnvkit_cnprep_clustering/")
source("./scripts/helperFunctions.R")
source("./scripts/segmentClusteringLibrary.R")
source("./scripts/cnvkitAdapterFunctions.R")

#
# Load input
#
normal_samples <- load_samples(classes = c("N"), sampleList = "./resources/sampleList.csv")

cytobands <- retrieveCytobands(dir = "./resources/cytoBand.txt")
chromosomeSizes <- readRDS("./resources/chromosomeSizes.rds")

normalSegments <- do.call(rbind, lapply(normal_samples, function(sample){
  return(retrieveCNVkitSegments(sample, dir="./resources/testCNVkitNorm/", genes = FALSE))
}))

normalSegments <- cnvKitSegmentsToBedFormat(normalSegments)

target_samples <- load_samples(classes = c("T", "F", "M"), sampleList = "./resources/sampleList.csv")

# Generate norminput argument
norminput <- retrieveNormInput(normalSegments)
norminput <- filterNormInput(norminput, length_threshold=1.1e7)

# Create folder with output
dir.create(file.path("./output/", output_dir), showWarnings = FALSE)


for(target_samples.i in seq(1, length(target_samples))){
  sample <- target_samples[target_samples.i]
  
  print(paste("Analyzing sample", sample))
  
  #
  # Retrieve sample data
  #
  cnvkit_segment_data <- retrieveCNVkitSegments(sample, dir="./resources/testCNVkit/", genes = FALSE)
  cnvkit_bins_data <- retrieveCNVkitBins(sample, dir="./resources/testCNVkit/")
  
  segment_probes_assignment <- assignProbeNumbers(cnvkit_segment_data, cnvkit_bins_data)
  
  # Generate seginput argument
  seginput <- retrieveSegInput(cnvkit_segment_data = cnvkit_segment_data, segment_probes_assignment = segment_probes_assignment, sample=sample, chromosomeSizes = chromosomeSizes, cytobands = cytobands)
  print(paste("Retrieved segment input for sample", sample))
  
  # Generate ratinput argument
  ratinput <- retrieveRatInput(cnvkit_bins_data, sample)
  print(paste("Retrieved ratio input for sample", sample))  
  
  # Run CNprep:CNpreprocessing
  try({
  	segtable <- runCNpreprocessing(seginput = seginput, ratinput = ratinput, norminput = norminput, modelNames = mclust_model, minjoin = minjoin, ntrial = ntrial) #TODO: Is there a distrib="Grid"?
	  print(paste("Produced segtable for sample", sample))
  
	  write.table(segtable, paste("./output/", output_dir,"/", sample, "_segtable.tsv", sep = ""), row.names = F, sep = "\t", quote = FALSE)
	  print(paste("Wrote output for sample", sample))
  }, silent = TRUE)
}

# TODO: Do this for all HPC scripts. May need to make function in helperFunctions.R to reduce duplicate code
info_filename <- paste("./output/", output_dir, "/JobInformation.txt", sep = "")
file.create(info_filename)
fileConn <- file(info_filename)
writeLines(	c("UGE JOB SUBMISSION NOTES FOR segmentClusteringHPC.R",
	 	paste("User description of job: ", description),
		"The script centers the input segments and clusters the segments using GMM",
		paste("The output files wrote to: ", output_dir),
		paste("This input parameters used mclust_model =", mclust_model, "and minjoin =", minjoin, "and ntrial=", ntrial)),
	 fileConn)
close(fileConn)
