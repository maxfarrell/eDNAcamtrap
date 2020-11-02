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
# 0% 15, 25% 250, 50% 251, 75% 251, 100% 301
# hist(as.numeric(unlist(lapply(forward,nchar))))

reverse <- lapply(fnRs, getSequences)
quantile(as.numeric(unlist(lapply(reverse,nchar))))
# 0% 15, 25% 225, 50% 225, 75% 225, 100% 301
# hist(as.numeric(unlist(lapply(reverse,nchar))))

sum(as.numeric(unlist(lapply(reverse,nchar)))<220)# 15,916
sum(as.numeric(unlist(lapply(reverse,nchar)))<225)# 380,485
sum(as.numeric(unlist(lapply(reverse,nchar)))<230)# 2,051,237

sum(as.numeric(unlist(lapply(forward,nchar)))<200)# 199,104

# Mean length is 473
mean(as.numeric(unlist(lapply(forward,nchar))))# 241.5
mean(as.numeric(unlist(lapply(reverse,nchar))))# 232.5
241 + 232

# Unsure what length of sequence the mini-barcode primers are targeting...

pdf("../plots_tables/mam_ckt_230_quality_profiles.pdf")
plotQualityProfile(fnFs, aggregate=TRUE) 
plotQualityProfile(fnRs, aggregate=TRUE) 
dev.off()

# Looks like quality is okay forward reads, but decreases in reverse 

# Forward: 230
# Reverse: 140
230+140

### Perform filtering & trimming ###

# Assign the filenames for the filtered fastq.gz files.
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter the forward and reverse reads
# remove N bases (dada2 requirement)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=5, 
				truncLen=c(220, 140),
              	maxN=0, truncQ=6, rm.phix=TRUE, maxEE=c(4,4),
              	compress=TRUE, multithread=3) # On Windows set multithread=FALSE
out
saveRDS(out, "mam_ckt_230_filterAndTrim_out.rds")

sum(out[,1])-sum(out[,2])
(sum(out[,1])-sum(out[,2]))/sum(out[,1])
# 470,790 sequences (17.7%) lost with c(220, 140)
# 480,622 sequences (18 %) lost with c(220, 160)
# 516,237 sequences (19 %) lost with c(220, 180)
.

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
Batch12.filtRs <- filtRs[filtFs.sample.names %in% Batch12.sample.names]
names(Batch12.filtFs) <- Batch12.sample.names
names(Batch12.filtRs) <- Batch12.sample.names
Batch15.filtFs <- filtFs[filtFs.sample.names %in% Batch15.sample.names]
Batch15.filtRs <- filtRs[filtFs.sample.names %in% Batch15.sample.names]
names(Batch15.filtFs) <- Batch15.sample.names
names(Batch15.filtRs) <- Batch15.sample.names
Batch54.filtFs <- filtFs[filtFs.sample.names %in% Batch54.sample.names]
Batch54.filtRs <- filtRs[filtFs.sample.names %in% Batch54.sample.names]
names(Batch54.filtFs) <- Batch54.sample.names
names(Batch54.filtRs) <- Batch54.sample.names
Batch7.filtFs <- filtFs[filtFs.sample.names %in% Batch7.sample.names]
Batch7.filtRs <- filtRs[filtFs.sample.names %in% Batch7.sample.names]
names(Batch7.filtFs) <- Batch7.sample.names
names(Batch7.filtRs) <- Batch7.sample.names

Batches <- c("Batch12","Batch15","Batch54","Batch7")

# output format example: "/path/to/run1/output/seqtab.rds"

learn_derep_merge_save <- function(filt_Fs, filt_Rs, Batch.sample.names, output){

	# Learn forward error rates
	errF <- learnErrors(filt_Fs, multithread=TRUE)
	# Learn reverse error rates
	errR <- learnErrors(filt_Rs, multithread=TRUE)

	# save error rates for inspection later	
	saveRDS(errF, paste("Batch",length(Batch.sample.names),"errF","mam_ckt_230.rds",sep="_"))
	saveRDS(errR, paste("Batch",length(Batch.sample.names),"errR","mam_ckt_230.rds",sep="_"))

	# objects to store output (for tracking read numbers)
	mergers <- vector("list", length(Batch.sample.names))
	names(mergers) <- Batch.sample.names
	ddFs <- vector("list", length(Batch.sample.names))
	names(ddFs) <- Batch.sample.names
	
	# Sample inference and merger of paired-end reads
	for(sam in Batch.sample.names) {
	  cat("Processing:", sam, "\n")
	    derepF <- derepFastq(filt_Fs[[sam]])
	    ddF <- dada(derepF, err=errF, multithread=TRUE, pool=TRUE)
	    ddFs[[sam]] <- ddF
	    derepR <- derepFastq(filt_Rs[[sam]])
	    ddR <- dada(derepR, err=errR, multithread=TRUE)
	    # merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = 12,  maxMismatch = 2)
	    merger <- mergePairs(ddF, derepF, ddR, derepR, justConcatenate=TRUE)
	    mergers[[sam]] <- merger
	}
	rm(derepF); rm(derepR)

	# Construct and save sequence table
	seqtab <- makeSequenceTable(mergers)
	saveRDS(seqtab, output)

	# Save ddFs and mergers for tracking read numbers
	saveRDS(ddFs, paste("Batch",length(Batch.sample.names),"ddFs","mam_ckt_230.rds",sep="_"))
	saveRDS(mergers, paste("Batch",length(Batch.sample.names),"mergers","mam_ckt_230.rds",sep="_"))

}

learn_derep_merge_save(Batch12.filtFs, Batch12.filtRs, Batch12.sample.names, "./Batch12_seqtab_mam_ckt_230.rds")
learn_derep_merge_save(Batch15.filtFs, Batch15.filtRs, Batch15.sample.names, "./Batch15_seqtab_mam_ckt_230.rds")
learn_derep_merge_save(Batch54.filtFs, Batch54.filtRs, Batch54.sample.names, "./Batch54_seqtab_mam_ckt_230.rds")
learn_derep_merge_save(Batch7.filtFs, Batch7.filtRs, Batch7.sample.names, "./Batch7_seqtab_mam_ckt_230.rds")

### Merge Runs & Remove Chimeras ###
# This part is the same as in the single-read version of this workflow.

# Merge multiple runs 
st1 <- readRDS("./Batch12_seqtab_mam_ckt_230.rds")
st2 <- readRDS("./Batch15_seqtab_mam_ckt_230.rds")
st3 <- readRDS("./Batch54_seqtab_mam_ckt_230.rds")
st4 <- readRDS("./Batch7_seqtab_mam_ckt_230.rds")
st.all <- mergeSequenceTables(st1, st2, st3, st4)

# Remove chimeras
seqtab.nochim.consensus <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Once the full-study sequence table is created, chimeras can be identified 
# and removed, and taxonomy assigned. For chimera removal, we have found 
# that the "consensus" chimera removal method works better on large studies, 
# but the "pooled" method is also an option.
seqtab.nochim.pooled <- removeBimeraDenovo(st.all, method="pooled", multithread=TRUE)

# Write to disk
saveRDS(seqtab.nochim.pooled, "./coi_mam_ckt_230_seqtab.rds")


## Track reads

# reading back in "out"
out <- readRDS("mam_ckt_230_filterAndTrim_out.rds")
rownames(out) <- sample.names

# Bringing in dada and merger read counts for tracking
ddFs1 <- readRDS("./Batch_12_ddFs_mam_ckt_230.rds")
ddFs2 <- readRDS("./Batch_15_ddFs_mam_ckt_230.rds")
ddFs3 <- readRDS("./Batch_54_ddFs_mam_ckt_230.rds")
ddFs4 <- readRDS("./Batch_7_ddFs_mam_ckt_230.rds")
ddFs <- c(ddFs1, ddFs2, ddFs3, ddFs4)

mergers1 <- readRDS("./Batch_12_mergers_mam_ckt_230.rds")
mergers2 <- readRDS("./Batch_15_mergers_mam_ckt_230.rds")
mergers3 <- readRDS("./Batch_54_mergers_mam_ckt_230.rds")
mergers4 <- readRDS("./Batch_7_mergers_mam_ckt_230.rds")
mergers <- c(mergers1, mergers2, mergers3, mergers4)

# Must re-order ddFs and mergers by out to properly track reads
# Due to batch processing of samples
ddFs <- ddFs[rownames(out)]
mergers <- mergers[rownames(out)]
# seqtab.nochim.pooled is reordered when calling track

getN <- function(x) sum(getUniques(x))
track <- cbind(out, unlist(lapply(ddFs, getN)), sapply(mergers, getN), rowSums(seqtab.nochim.pooled[rownames(out),]))
colnames(track) <- c("input", "filtered", "denoised", "merged", "nonchim.pooled")
rownames(track) <- sample.names

require(xtable)
print(xtable(track, digits=0),file="./coi_mam_ckt_230_track_reads.tex")

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
print(xtable(track_cam, digits=0),file="./coi_mam_ckt_230_track_reads_cam_samples.tex")
colSums(track_cam)
# input 	1315657
# filtered 	1064124
# denoised	1028471
# merged 	1018900
# nonchim 	982248

seqtab <- readRDS("../results/coi_mam_ckt_230_seqtab.rds")
seqtab_cam <- seqtab[rownames(seqtab)%in%cam_samples,]
seqtab_cam <- seqtab_cam[,colSums(seqtab_cam)>0]
dim(seqtab_cam) # 46 x 2963

