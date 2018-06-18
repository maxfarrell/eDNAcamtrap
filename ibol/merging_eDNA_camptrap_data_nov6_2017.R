# Merging Camera Trap eDNA data

# Max Farrell (maxwellfarrell@gmail.com)

# Last updated Noember 6th 2017

rm(list=ls())

require(lubridate)
require(dplyr)
require(tidyr)
require(ggplot2)

### Kruger metadata
env <- read.csv("../../Data/Kruger_2015_Environmental.csv", as.is=T)
names(env) <- tolower(names(env))
camNotes <- read.csv("../../Data/Kruger_2015_Camera Notes.csv", as.is=T)
names(camNotes) <- tolower(names(camNotes))
names(camNotes)[names(camNotes)=="notes"] <- "cam_notes"
labNotes <- read.csv("../../Data/Kruger_2015_Water Samples.csv", as.is=T)
names(labNotes) <- tolower(names(labNotes))

### Camera trap data
trap <- read.csv("../../Data/Trap Annotation/annotation_clean_dec12.csv", as.is=T)
trap$timestamp <- ymd_hms(trap$timestamp)

# Matching Site names
trap$site[trap$site=="Kwaggas"] <- "Kwaggas Pan"
trap$site[trap$site=="Nyamahri"] <- "Nyamarhi"
env$site[env$site=="Nyamahri"] <- "Nyamarhi"
camNotes$site[camNotes$site=="De LaPorte"] <- "De Laporte"
labNotes$site[labNotes$site=="Nhlanguneli"|labNotes$site=="Nhanguleni"] <- "Nhlanguleni"
labNotes$site[labNotes$site=="Nyamahri"] <- "Nyamarhi"

### megablast taxonomy assignments with refLib_v7
# copied from assigned_tax_loop_over_pident.R

# Formatting the taxonomy file
tax <- read.delim("../../Data/Molecular/eDNA_camtrap/classifications/megablast_refLib_v7/KNP_refTaxonomy_v3.tsv", sep="\t", as.is=TRUE, header=FALSE)
str(tax)
tax_list <- strsplit(tax[,2], "[;]")
maxlevels <- 6
taxlevels <- as.data.frame(t(sapply(tax_list,'[',1:maxlevels)), stringsAsFactors=FALSE)
tax <- data.frame(tax,taxlevels)
tax <- tax[,-2]
names(tax) <- c("sseqid","phylum", "class","order","family","genus","species")

# Reading in assignment
assign <- read.csv("../../Data/Molecular/eDNA_camtrap/classifications/megablast_refLib_v7/combined_seqs_KNP_refLib_v7.csv", header=FALSE, stringsAsFactors=FALSE)
names(assign) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart" ,"send" ,"evalue", "bitscore")

# Merging assignments with taxonomy
assign <- left_join(assign, tax)
seqXspec <- unique(assign[,c("qseqid","species")])
length(unique(assign$qseqid))#169548
dim(seqXspec) #203201
dim(seqXspec)[1] - length(unique(assign$qseqid))
# 33653 sequences had discrepancies 
discrep <- seqXspec$qseqid[duplicated(seqXspec$qseqid)]
head(assign[assign$qseqid%in%discrep,])

# Cleaning up qseqid to get sample information
assign_temp <- gsub("[.]","_", assign$qseqid)
assign$site <- sapply(strsplit(assign_temp, "_"), function(x) x[1])
assign$number <- gsub("[A-Z]","",sapply(strsplit(assign_temp, "_"), function(x) x[2]))
assign$AB <- gsub("[0-9]","",sapply(strsplit(assign_temp, "_"), function(x) x[2])) 
assign$S_XS <- sapply(strsplit(assign_temp, "_"), function(x) x[3])
assign$S_XS[grep("[0-9]",assign$S_XS)] <- NA
assign$sample_num <- paste(assign$site, assign$number, sep="_")

str(assign)

## Restricting by pident and length: >97% identity
assign <- assign[assign$pident>=97,]
quantile(assign$length)
# 50% is 29 - there are still many short sequences with high identity

assign <- assign[assign$length>=100,]
quantile(assign$length)
# 0% is 190; 25% is 208

# redoing seqXspec
seqXspec <- unique(assign[,c("qseqid","species")])
dim(seqXspec) #1821
discrep <- seqXspec$qseqid[duplicated(seqXspec$qseqid)]
# now only 6 sequences had discrepancies
unique(assign$species[assign$qseqid%in%discrep])
# can't differentiate between panthera leo and panthera pardus
assign[assign$qseqid%in%discrep,c("qseqid","pident","length","species","evalue")]

# In this case it looks like this is for one sample and discrepancy between lion and leopard
# fortunately all leopard hits are below 98% identity, so we can set a higher threshold for these:
assign[assign$qseqid%in%discrep,c("qseqid","pident","length","species","evalue")]
# To keep
assign[which(assign$qseqid%in%discrep & assign$pident>98),]
# To remove
assign[which(assign$qseqid%in%discrep & assign$pident<98),]
# Removing
assign <- assign[-which(assign$qseqid%in%discrep & assign$pident<98),]

# redoing seqXspec
seqXspec <- unique(assign[,c("qseqid","species")])
dim(seqXspec) #1815
discrep <- seqXspec$qseqid[duplicated(seqXspec$qseqid)]
unique(assign$species[assign$qseqid%in%discrep]) # 0


# Making labNotes sample notation match annotation data
colnames(labNotes)[colnames(labNotes)=="sample_code"] <- "sample"
labNotes$site <- sapply(strsplit(labNotes$sample, "_"), function(x) x[1])
labNotes$number <- gsub("[A-Z]","",sapply(strsplit(labNotes$sample, "_"), function(x) x[2]))
labNotes$AB <- gsub("[0-9]","",sapply(strsplit(labNotes$sample, "_"), function(x) x[2])) 
labNotes$S_XS <- sapply(strsplit(labNotes$sample, "_"), function(x) x[3])
labNotes$sample_num <- paste(labNotes$site, labNotes$number, sep="_")

# Creating file with date ranges to subset 
sample_dates <- subset(labNotes, select=c("site","sample_num","date_collected","time_collected"))
sample_dates$collected <- as.POSIXct(with(sample_dates, paste("2015", date_collected, time_collected, sep=":")), format="%Y:%B_%d:%H:%M", tz="UTC")
sample_dates <- subset(sample_dates, select =-c(date_collected, time_collected))

# Getting a single time of collection for each sample_num
# using min time rather than mean so as to get largest trap length
sample_dates <- unique(sample_dates %>%
				group_by(sample_num) %>%
				summarize(collected = min(collected)))

# Changing trap site name to match seqs site name
sort(unique(trap$site))
sort(unique(assign$site))

lookup <- setNames(sort(unique(assign$site)), sort(unique(trap$site)))
trap$site <- lookup[trap$site]

# Changing env site name to match seqs site name
env$site <- lookup[env$site]
# Subsetting to only include sites for camtrap analysis
env <- env[!is.na(env$site),]

# Joining date collected with seqs
seqs <- left_join(assign, sample_dates)
sample_dates$site <- gsub("_[0-9]","",sample_dates$sample_num)
sample_dates$weekstart <- sample_dates$collected %m+%days(-7)

days_sampled <- with(sample_dates, collected-weekstart)

# Sites / Samples
# Blanks:
# BLANK_2
# BLANK_3

# Main coverage (bookended samples):
# HOY 2_3 (June 22-29) (+ 3AB S/XS)
# KWA 5_6 (June 19-26) (+ 6AB S/XS)
# NGO 2_3 (June 24-July 1) (+ 3AB S/XS)
# NYA 2_3_4 (June 24-July 1; July 1-8) (+ NYA 4AB S/XS)

# Extra coverage (full week):
# DLP 7_8 (July 3-10) (only 8AB sequenced; + S/XS) (14 days of cam trapping)
# NGO 1_2 (June 18-24) # Can be used as a full week of camera data...
# KWA 1_5 (June 15-19)
# NYA 1_2 (June 18-24) # Can be used as a full week of camera data...

# # Full week of camera data + the trough was opened at the beginning of the week (was dry before sampling)
# NWA 1_2 (June 19-26) (only NWA 2AB S/XS sequenced)

# NGO_2 and NYA_2 include only 6 days of sampling instead of 7...
# Because they were set/checked on a Saturday and sampled on a Friday...

camera_samples_fullweek <- c("DLP_8","HOY_3","KWA_6","NGO_2","NGO_3",
								"NWA_2","NYA_2", "NYA_3","NYA_4")

sample_dates <- sample_dates[sample_dates$sample_num%in%camera_samples_fullweek,]

trap$sample_num <- NULL

for (i in seq_along(trap$timestamp)) {

	# print(i)

	if (sum(sample_dates$site==trap$site[i] & 
		sample_dates$weekstart < trap$timestamp[i] & 
		sample_dates$collected > trap$timestamp[i])==1){

	trap$sample_num[i] <- sample_dates$sample_num[which(sample_dates$site==trap$site[i] & 
						sample_dates$weekstart < trap$timestamp[i] & 
						sample_dates$collected > trap$timestamp[i])]
	}else{trap$sample_num[i] <- "REMOVE"}
}

# NGO_2 and NYA_2 still missing only 6 days instead of 7...
# Because they were set/checked on a Saturday and sampled on a Friday...

sp_to_remove <- c("Accipiter_ovampensis","Aegypius_occipitalis",
					"Aegypius_tracheliotos","Buphagus_erythrorhynchus",
					"Corythaixoides_concolor","Homo_sapiens","Necrosyrtes_monachus",
					"Streptopelia_senegaliensis","Struthio_camelus",
					"Terathopius_ecaudataus")

trap <- trap[!trap$binomial%in% sp_to_remove,]
trap <- trap[!is.na(trap$binomial),]
trap <- trap[trap$sample_num!="REMOVE",]

# There are some samples in seqs that are not in trap
seqs <- seqs[seqs$sample_num%in%trap$sample_num,]

rm(list=c("assign","assign_temp","camera_samples_fullweek","days_sampled",
			"discrep","i","lookup","maxlevels","sample_dates","seqXspec","sp_to_remove",
			"taxlevels","tax_list"))

# save.image("merged_eDNA_camtrap_data_nov6_2017.RData")
