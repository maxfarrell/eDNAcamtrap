# coi mam_ckt_230 pipeline
# for stu

###########################################
## dada2 quality filtering & ASV calling ##
###########################################

# Moving to R
# in order to visualize plots make sure to use ssh -X <server>

# dada2 pipeline tutorial
# https://benjjneb.github.io/dada2/tutorial.html
# adaptation for v3-v4 https://github.com/benjjneb/dada2/issues/227
# additonal comment on v3-v4 https://github.com/benjjneb/dada2/issues/319
# "It's the coi sequencing run of V3/V4 region, so the reads are 250 nts each 
# and the overlap is about 60 nts so they merge to about 440 nts."

# sudo R #needs to be sudo to augment file structure
require("dada2"); packageVersion("dada2") # 1.10.0

path <- "/data/kruger/illumina/coi_mam_ckt_230_fastq" # directory containing the fastq files
list.files(path)

# Set random seed (for reproducibility)
set.seed(4534)

### Filter & Trim ###

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.pick.trim.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.pick.trim.fastq", full.names = TRUE))

# Extract and format sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
kruger.sample.names <- gsub("Max[-]","", sample.names)
kruger.sample.names <- gsub("[-]","_", kruger.sample.names)

### Examine quality profiles of forward and reverse reads ###
forward <- lapply(fnFs, getSequences)
quantile(as.numeric(unlist(lapply(forward,nchar))))
# 0% 13, 25% 230, 50% 250, 75% 251, 100% 301
# hist(as.numeric(unlist(lapply(forward,nchar))))

reverse <- lapply(fnRs, getSequences)
quantile(as.numeric(unlist(lapply(reverse,nchar))))
# 0% 8, 25% 224, 50% 228, 75% 228, 100% 301
# hist(as.numeric(unlist(lapply(reverse,nchar))))

# Looks like we should trim seqs below 200bp
sum(as.numeric(unlist(lapply(reverse,nchar)))<50)# 20,386
sum(as.numeric(unlist(lapply(reverse,nchar)))<100)# 79,832
sum(as.numeric(unlist(lapply(reverse,nchar)))<200)# 152,986
sum(as.numeric(unlist(lapply(reverse,nchar)))<230)# 1,802,426

# Mean length is 460
mean(as.numeric(unlist(lapply(forward,nchar))))# 231.28
mean(as.numeric(unlist(lapply(reverse,nchar))))# 229.25
231 + 229

# Unsure what length of sequence the mini-barcode primers are targeting...
# Gibson 2015 (mam_ckt_230) designed to amplify 310bp fragment
# Should have lots of overlap for merging

pdf("../plots_tables/mam_ckt_230_quality_profiles.pdf")
plotQualityProfile(fnFs, aggregate=TRUE) 
plotQualityProfile(fnRs, aggregate=TRUE) 
dev.off()

# Looks like quality is okay forward reads, but decreases in reverse 

# Moving on with only Forward Reads

### Perform filtering & trimming ###

# Assign the filenames for the filtered fastq.gz files.
filt_path <- file.path(path, "filtered_Fonly") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt_Fonly.fastq.gz"))

# Filter the forward reads
out <- filterAndTrim(fnFs, filtFs, truncLen=c(220),
              maxN=0, truncQ=6, rm.phix=TRUE, maxEE=4,
              compress=TRUE, multithread=3) # On Windows set multithread=FALSE
out
saveRDS(out, "mam_ckt_230_filterAndTrim_Fonly_out.rds")

sum(out[,1])-sum(out[,2])
(sum(out[,1])-sum(out[,2]))/sum(out[,1])
# 489,603 sequences (18.39%) lost 

### Infer Sequence Variants ###
# Infer sequence variants in each sample
# NOTE: This portion of the workflow should be performed on a run-by-run basis, as error rates can differ between runs.
# See - https://benjjneb.github.io/dada2/bigdata_paired.html

# Runs - Batch12, Batch15, Batch54, Batch7
# Get names for each batch and match to pick.trim.fastq files in coi_fastq/
Batch12Fs <- sort(list.files("/data/kruger/illumina/Batch12", pattern="_R1_001.fastq", full.names = TRUE))
Batch15Fs <- sort(list.files("/data/kruger/illumina/Batch15", pattern="_R1_001.fastq", full.names = TRUE))
Batch54Fs <- sort(list.files("/data/kruger/illumina/Batch54", pattern="_R1_001.fastq", full.names = TRUE))
Batch7Fs <- sort(list.files("/data/kruger/illumina/Batch7", pattern="_R1_001.fastq", full.names = TRUE))

# Extract and format sample names
Batch12.sample.names <- sapply(strsplit(basename(Batch12Fs), "_"), `[`, 1)
Batch15.sample.names <- sapply(strsplit(basename(Batch15Fs), "_"), `[`, 1)
Batch54.sample.names <- sapply(strsplit(basename(Batch54Fs), "_"), `[`, 1)
Batch7.sample.names <- sapply(strsplit(basename(Batch7Fs), "_"), `[`, 1)

# Workflow for Big Data #
# Subsetting filtFs and filtRs per Batch
filtFs.sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)

Batch12.filtFs <- filtFs[filtFs.sample.names %in% Batch12.sample.names]
names(Batch12.filtFs) <- Batch12.sample.names
Batch15.filtFs <- filtFs[filtFs.sample.names %in% Batch15.sample.names]
names(Batch15.filtFs) <- Batch15.sample.names
Batch54.filtFs <- filtFs[filtFs.sample.names %in% Batch54.sample.names]
names(Batch54.filtFs) <- Batch54.sample.names
Batch7.filtFs <- filtFs[filtFs.sample.names %in% Batch7.sample.names]
names(Batch7.filtFs) <- Batch7.sample.names

Batches <- c("Batch12","Batch15","Batch54","Batch7")

# output format example: "/path/to/run1/output/seqtab.rds"

learn_derep_save <- function(filt_Fs, Batch.sample.names, output){

	# Learn forward error rates
	errF <- learnErrors(filt_Fs, multithread=TRUE)

	# object to store output (for tracking read numbers)
	ddFs <- vector("list", length(Batch.sample.names))
	names(ddFs) <- Batch.sample.names
	
	# Sample inference and merger of paired-end reads
	for(sam in Batch.sample.names) {
	  cat("Processing:", sam, "\n")
	    derepF <- derepFastq(filt_Fs[[sam]])
	    ddF <- dada(derepF, err=errF, multithread=TRUE, pool=TRUE)
	    ddFs[[sam]] <- ddF
	}
	rm(derepF)

	# Construct and save sequence table
	seqtab <- makeSequenceTable(ddFs)
	saveRDS(seqtab, output)

	# Save ddFs and mergers for tracking read numbers
	saveRDS(ddFs, paste("Batch",length(Batch.sample.names),"ddFs","mam_ckt_230","Fonly.rds",sep="_"))
}

learn_derep_save(Batch12.filtFs, Batch12.sample.names, "./Batch12_seqtab_mam_ckt_230_Fonly.rds")
learn_derep_save(Batch15.filtFs, Batch15.sample.names, "./Batch15_seqtab_mam_ckt_230_Fonly.rds")
learn_derep_save(Batch54.filtFs, Batch54.sample.names, "./Batch54_seqtab_mam_ckt_230_Fonly.rds")
learn_derep_save(Batch7.filtFs, Batch7.sample.names, "./Batch7_seqtab_mam_ckt_230_Fonly.rds")


### Merge Runs & Remove Chimeras ###
# This part is the same as in the single-read version of this workflow.

# Merge multiple runs 
st1 <- readRDS("./Batch12_seqtab_mam_ckt_230_Fonly.rds")
st2 <- readRDS("./Batch15_seqtab_mam_ckt_230_Fonly.rds")
st3 <- readRDS("./Batch54_seqtab_mam_ckt_230_Fonly.rds")
st4 <- readRDS("./Batch7_seqtab_mam_ckt_230_Fonly.rds")
st.all <- mergeSequenceTables(st1, st2, st3, st4)

# Remove chimeras
seqtab.nochim.pooled <- removeBimeraDenovo(st.all, method="pooled", multithread=TRUE)

# Save final seqtab
saveRDS(seqtab.nochim.pooled, "./coi_mam_ckt_230_Fonly_seqtab.rds")


## Track reads
# reading back in "out"
out <- readRDS("mam_ckt_230_filterAndTrim_Fonly_out.rds")
rownames(out) <- sample.names

# Bringing in dada and merger read counts for tracking
ddFs1 <- readRDS("./Batch_12_ddFs_mam_ckt_230_Fonly.rds")
ddFs2 <- readRDS("./Batch_15_ddFs_mam_ckt_230_Fonly.rds")
ddFs3 <- readRDS("./Batch_54_ddFs_mam_ckt_230_Fonly.rds")
ddFs4 <- readRDS("./Batch_7_ddFs_mam_ckt_230_Fonly.rds")
ddFs <- c(ddFs1, ddFs2, ddFs3, ddFs4)

# Must re-order ddFs by out to properly track reads
# Due to batch processing of samples
ddFs <- ddFs[rownames(out)]
# seqtab.nochim.pooled is reordered when calling track

getN <- function(x) sum(getUniques(x))
track <- cbind(out, unlist(lapply(ddFs, getN)), rowSums(seqtab.nochim.pooled[rownames(out),]))
colnames(track) <- c("input", "filtered", "denoised", "nonchim.pooled")
rownames(track) <- sample.names

require(xtable)
print(xtable(track, digits=0),file="./plots_tables/coi_mam_ckt_230_Fonly_track_reads.tex")

# Subset read table to those in camera trap analyis
cam_samples <- c("Max-BLANK-2","Max-BLANK-3","Max-DLP-8A","Max-DLP-8A-S",
				 "Max-DLP-8A-XS","Max-DLP-8B","Max-DLP-8B-S","Max-DLP-8B-XS",
				 "Max-HOY-2A","Max-HOY-2B","Max-HOY-3A","Max-HOY-3A-S",
				 "Max-HOY-3A-XS","Max-HOY-3B","Max-HOY-3B-S","Max-HOY-3B-XS",
				 "Max-KWA-6A","Max-KWA-6A-S","Max-KWA-6A-XS","Max-KWA-6B",
				 "Max-KWA-6B-S","Max-KWA-6B-XS","Max-NGO-2A","Max-NGO-2B",
				 "Max-NGO-3A","Max-NGO-3A-S","Max-NGO-3A-XS","Max-NGO-3B",
				 "Max-NGO-3B-S","Max-NGO-3B-XS","Max-NWA-2A","Max-NWA-2A-S",
				 "Max-NWA-2A-XS","Max-NWA-2B","Max-NWA-2B-S","Max-NWA-2B-XS",
				 "Max-NYA-2A","Max-NYA-2B","Max-NYA-3A","Max-NYA-3B",
				 "Max-NYA-4A","Max-NYA-4A-S","Max-NYA-4A-XS","Max-NYA-4B",
				 "Max-NYA-4B-S","Max-NYA-4B-XS")

track_cam <- track[rownames(track)%in%cam_samples,]
print(xtable(track_cam, digits=0),file="../plots_tables/coi_mam_ckt_230_Fonly_track_reads_cam_samples.tex")
colSums(track_cam)
# input 	1,315,657
# filtered 	1,055,439
# denoised	1,015,899
# nonchim 	973,592

seqtab <- readRDS("../results/coi_mam_ckt_230_Fonly_seqtab.rds")
seqtab_cam <- seqtab[rownames(seqtab)%in%cam_samples,]
seqtab_cam <- seqtab_cam[,colSums(seqtab_cam)>0]
dim(seqtab_cam) # 46 x 2636
