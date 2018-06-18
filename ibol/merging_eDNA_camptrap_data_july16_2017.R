# Merging Camera Trap eDNA data

# Max Farrell (maxwellfarrell@gmail.com)

# Last updated July 16th 2017

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

### megablast_refLib_v7_noOTU length > 100 abd cutoffs 90-100 by 0.05
# produced with assigned_tax_loop_over_pident.R
cutoffs_full <- readRDS("../../Data/Molecular/eDNA_camtrap/clean_data/COI_len100_cutoffs90_100.RDS")
str(cutoffs_full[[1]])
sort(unique(cutoffs_full[[1]]$site))

# MISSING BLANKS... WHY?
# BLANKS only have pident ranging from 74.10 - 85.07 
# (see "assigned_tax_loop_over_pident.R")

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
site_names <- sort(unique(unlist(lapply(cutoffs_full, function(x) unique(x$site)))))
lookup <- setNames(site_names, sort(unique(trap$site)))
trap$site <- lookup[trap$site]

# Changing env site name to match seqs site name
env$site <- lookup[env$site]
# Subsetting to only include sites for camtrap analysis
env <- env[!is.na(env$site),]

# Joining date collected with seqs
seqs <- lapply(cutoffs_full, function(x) left_join(x, sample_dates))

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


# test <- trap %>%
# 		group_by(sample_num)%>%
# 		mutate(timesince = as.numeric(difftime(max(timestamp),timestamp, units="secs")))

# test2 <- test %>%
# 			group_by(sample_num)%>%
# 			summarize(camtime=(max(timesince)-min(timesince))/(60*60*24))
# test2

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

rm(list=c("camera_samples_fullweek","site_names","i","lookup","sp_to_remove"))

# save.image("merged_eDNA_camtrap_data_july16_2017.RData")
