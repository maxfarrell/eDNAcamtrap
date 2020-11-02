# Merging Data

require(dada2); packageVersion("dada2") # 1.10.0
require(phyloseq); packageVersion("phyloseq") # 1.26.0
require(ggplot2); packageVersion("ggplot2") # 3.1.0
require(ape); packageVersion("ape")# 5.2
require(dplyr); packageVersion("dplyr")# 0.8.0.1
require(vegan); packageVersion("vegan")# 2.5.3
require(lubridate); packageVersion("lubridate")# 1.7.4
require(tidyr); packageVersion("tidyr")#0.8.2

# Merging data

# COI Sequence Table
seqtab <- readRDS("../results/coi_merged_Fonly_seqtab.rds")
dim(seqtab)# 88 x 9671

# tax assignments
taxKNP_vert_98 <- readRDS("../results/coi_vertebrates_dada2_phy2species_98.rds")
taxKNP_vert_95 <- readRDS("../results/coi_vertebrates_dada2_phy2species_95.rds")
taxKNP_vert_90 <- readRDS("../results/coi_vertebrates_dada2_phy2species_90.rds")
taxKNP_vert_80 <- readRDS("../results/coi_vertebrates_dada2_phy2species_80.rds")
taxKNP_vert_70 <- readRDS("../results/coi_vertebrates_dada2_phy2species_70.rds")
taxKNP_vert_60 <- readRDS("../results/coi_vertebrates_dada2_phy2species_60.rds")
taxKNP_vert_50 <- readRDS("../results/coi_vertebrates_dada2_phy2species_50.rds")
taxTP <- readRDS("../results/coi_TP_chordata_tax_80.rds")
taxMID <- readRDS("../results/coi_MIDORI_tax_80.rds")


# Subsetting to samples in camera trap analysis
cam_samples <- c("Max-BLANK-2","Max-BLANK-3","Max-DLP-8A","Max-DLP-8A-S",
				 "Max-DLP-8A-XS","Max-DLP-8B","Max-DLP-8B-S","Max-DLP-8B-XS",
				 "Max-HOY-3A","Max-HOY-3A-S",
				 "Max-HOY-3A-XS","Max-HOY-3B","Max-HOY-3B-S","Max-HOY-3B-XS",
				 "Max-KWA-6A","Max-KWA-6A-S","Max-KWA-6A-XS","Max-KWA-6B",
				 "Max-KWA-6B-S","Max-KWA-6B-XS","Max-NGO-2A","Max-NGO-2B",
				 "Max-NGO-3A","Max-NGO-3A-S","Max-NGO-3A-XS","Max-NGO-3B",
				 "Max-NGO-3B-S","Max-NGO-3B-XS","Max-NWA-2A","Max-NWA-2A-S",
				 "Max-NWA-2A-XS","Max-NWA-2B","Max-NWA-2B-S","Max-NWA-2B-XS",
				 "Max-NYA-2A","Max-NYA-2B","Max-NYA-3A","Max-NYA-3B",
				 "Max-NYA-4A","Max-NYA-4A-S","Max-NYA-4A-XS","Max-NYA-4B",
				 "Max-NYA-4B-S","Max-NYA-4B-XS")

seqtab <- seqtab[rownames(seqtab)%in%cam_samples,]
seqtab <- seqtab[,colSums(seqtab)>0]
dim(seqtab) # 44 x 4616
# colnames(seqtab)

taxKNP_vert_98 <- taxKNP_vert_98[rownames(taxKNP_vert_98)%in%colnames(seqtab),]
taxKNP_vert_95 <- taxKNP_vert_95[rownames(taxKNP_vert_95)%in%colnames(seqtab),]
taxKNP_vert_90 <- taxKNP_vert_90[rownames(taxKNP_vert_90)%in%colnames(seqtab),]
taxKNP_vert_80 <- taxKNP_vert_80[rownames(taxKNP_vert_80)%in%colnames(seqtab),]
taxKNP_vert_70 <- taxKNP_vert_70[rownames(taxKNP_vert_70)%in%colnames(seqtab),]
taxKNP_vert_60 <- taxKNP_vert_60[rownames(taxKNP_vert_60)%in%colnames(seqtab),]
taxKNP_vert_50 <- taxKNP_vert_50[rownames(taxKNP_vert_50)%in%colnames(seqtab),]
taxTP <- taxTP[rownames(taxTP)%in%colnames(seqtab),]
taxMID <- taxMID[rownames(taxMID)%in%colnames(seqtab),]


### Kruger metadata
env <- read.csv("../data/Kruger_2015_Environmental.csv", as.is=T)
names(env) <- tolower(names(env))
labNotes <- read.csv("../data/Kruger_2015_Water Samples.csv", as.is=T)
names(labNotes) <- tolower(names(labNotes))
sitelist <- read.csv("../data/Kruger_2015_SiteList.csv", as.is=T)
names(sitelist) <- tolower(names(sitelist))
names(sitelist)[3] <- "site"
camNotes <- read.csv("../data/Kruger_2015_Camera Notes.csv", as.is=T)
names(camNotes) <- tolower(names(camNotes))
names(camNotes)[names(camNotes)=="notes"] <- "cam_notes"

### Camera trap data
trap <- read.csv("../data/annotation_clean_dec12.csv", as.is=T)
length(unique(trap$filename))#16,027
trap$timestamp <- ymd_hms(trap$timestamp)

# Matching Site names
trap$site[trap$site=="Kwaggas"] <- "Kwaggas Pan"
trap$site[trap$site=="Nyamahri"] <- "Nyamarhi"
env$site[env$site=="Nyamahri"] <- "Nyamarhi"
camNotes$site[camNotes$site=="De LaPorte"] <- "De Laporte"
labNotes$site[labNotes$site=="Nhlanguneli"|labNotes$site=="Nhanguleni"] <- "Nhlanguleni"
labNotes$site[labNotes$site=="Nyamahri"] <- "Nyamarhi"

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
sample_dates <- sample_dates %>% 
				group_by(sample_num) %>% 
				filter(collected == min(collected)) %>%
				unique()

# Identifying BLANKS
blanks <- seqtab[grep("BLANK",rownames(seqtab)),]
samples.out <- rownames(seqtab)
sample <- gsub("Max-","", samples.out)
sample <- gsub("-","_", sample)
site <- sapply(strsplit(samples.out, "-"), `[`, 2)
num <- as.numeric(gsub("[A-z]","",sample))
AB <- sapply(strsplit(samples.out, "[0-9]"), `[`, 2)
AB <- gsub("-.*$","",AB)
S_XS <- sapply(strsplit(sample, "_"), `[`, 3)

# Adding sample_num to env
labNotes$site_time <- with(labNotes, paste(site, date_collected, sep="_"))
env$site_time <- with(env, paste(site, date, sep="_"))
lookup <- setNames(unique(labNotes$sample_num), unique(labNotes$site_time))
env$sample_num <- lookup[env$site_time]

sampdf <- data.frame(sample=as.character(sample), site=site, sample_num=paste(site,num,sep="_"),num=num, AB=AB, S_XS=S_XS, stringsAsFactors=FALSE)
rownames(sampdf) <- samples.out

# Merge into phyloseq object
colnames(seqtab) <- paste0("SV", 1:ncol(seqtab))

make_ps <- function(seqtab,sampdf,tax){

	rownames(tax) <- paste0("SV", 1:nrow(tax))
	
	ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(sampdf), 
               tax_table(tax))
	
}

ps_98 <- make_ps(seqtab,sampdf,taxKNP_vert_98)
ps_95 <- make_ps(seqtab,sampdf,taxKNP_vert_95)
ps_90 <- make_ps(seqtab,sampdf,taxKNP_vert_90)
ps_80 <- make_ps(seqtab,sampdf,taxKNP_vert_80)
ps_70 <- make_ps(seqtab,sampdf,taxKNP_vert_70)
ps_60 <- make_ps(seqtab,sampdf,taxKNP_vert_60)
ps_50 <- make_ps(seqtab,sampdf,taxKNP_vert_50)
ps_TP <- make_ps(seqtab,sampdf,taxTP)
ps_MID <- make_ps(seqtab,sampdf,taxMID)

sample_dates$site <- gsub("_[0-9]","",sample_dates$sample_num)
sample_dates$weekstart <- sample_dates$collected %m+%days(-7)
days_sampled <- with(sample_dates, collected-weekstart)

# Matching trap$site with sample_dates$site
sort(unique(trap$site))
sort(unique(sample_dates$site))
sort(unique(sample_dates$site))[c(2,4,6,7,9,10)]

lookup <- setNames(sort(unique(sample_dates$site))[c(2,4,6,7,9,10)], sort(unique(trap$site)))
trap$site <- lookup[trap$site]

# Changing env site name to match seqs site name
env$site <- lookup[env$site]
# Subsetting to only include sites for camtrap analysis
env <- env[!is.na(env$site),]

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

# removing birds & humans
# sort(unique(trap$binomial))
sp_to_remove <- c("Accipiter_ovampensis","Aegypius_occipitalis",
					"Aegypius_tracheliotos","Buphagus_erythrorhynchus",
					"Corythaixoides_concolor","Homo_sapiens","Necrosyrtes_monachus",
					"Streptopelia_senegaliensis","Struthio_camelus",
					"Terathopius_ecaudataus")

trap <- trap[!trap$binomial%in% sp_to_remove,]
trap <- trap[!is.na(trap$binomial),]
trap <- trap[trap$sample_num!="REMOVE",]

length(unique(trap$filename))
unique(trap$species)
trap_seen <- trap[trap$n_waterhole>0,]
unique(trap_seen$species)
# fox?
trap_seen %>% filter(species=='fox')
trap_seen %>% filter(species=='fox')

# remove fox and frog



# Adding sample_num to env
labNotes$site_time <- with(labNotes, paste(site, date_collected, sep="_"))
env$site_time <- with(env, paste(site, date, sep="_"))
lookup <- setNames(unique(labNotes$sample_num), unique(labNotes$site_time))
env$sample_num <- lookup[env$site_time]

# to_keep <- c("camera_samples_fullweek","camNotes","env","labNotes","ps","sample_dates","seqtab","sitelist","tax98","trap")
# rm(list=ls()[!ls()%in%to_keep])

save.image("../data/merged_eDNA_camtrap_data_nov20_2019.RData")

