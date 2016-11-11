# Merging Camera Trap eDNA data

# Max Farrell (maxwellfarrell@gmail.com)

# Last updated October 17th 2016

rm(list=ls())

require(plyr)
require(biom)
require(lubridate)

### Kruger metadata
env <- read.csv("../Data/Kruger_2015_Environmental.csv", as.is=T)
names(env) <- tolower(names(env))
camNotes <- read.csv("../Data/Kruger_2015_Camera Notes.csv", as.is=T)
names(camNotes) <- tolower(names(camNotes))
names(camNotes)[names(camNotes)=="notes"] <- "cam_notes"
labNotes <- read.csv("../Data/Kruger_2015_Water Samples.csv", as.is=T)
names(labNotes) <- tolower(names(labNotes))

### Camera trap data
trap <- read.csv("../Data/Trap Annotation/annotation_clean_sept30.csv", as.is=T)
trap$timestamp <- ymd_hms(trap$timestamp)

# Matching Site names
trap$site[trap$site=="Kwaggas"] <- "Kwaggas Pan"
trap$site[trap$site=="Nyamahri"] <- "Nyamarhi"
env$site[env$site=="Nyamahri"] <- "Nyamarhi"
camNotes$site[camNotes$site=="De LaPorte"] <- "De Laporte"
labNotes$site[labNotes$site=="Nhlanguneli"|labNotes$site=="Nhanguleni"] <- "Nhlanguleni"
labNotes$site[labNotes$site=="Nyamahri"] <- "Nyamarhi"

### Classified sequence data

# v6
# blast_KNP_v6 <- read_biom("../Data/Molecular/eDNA_camtrap/classifications/blast_refLib_v6/blast_KNP_refLib_v6.biom")
# blast_KNP_v6_tax <- read.delim("../Data/Molecular/eDNA_camtrap/classifications/blast_refLib_v6/combined_seqs.fna_rep_set_tax_assignments_KNP_refLib_v6.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)

# # Formatting the taxonomy file
# str(blast_KNP_v6_tax)
# blast_KNP_v6_tax_list <- strsplit(blast_KNP_v6_tax[,2], "[;]")
# maxlevels <- 6
# taxlevels <- as.data.frame(t(sapply(blast_KNP_v6_tax_list,'[',1:maxlevels)), stringsAsFactors=FALSE)
# blast_KNP_v6_tax <- cbind(blast_KNP_v6_tax[,-2],taxlevels)
# names(blast_KNP_v6_tax) <- c("seqid","evalue","taxID","phylum", "class","order","family","genus","species")

# comm_OTUs <- t(as.matrix(biom_data(blast_KNP_v6)))
# comm_OTUs <- as.data.frame(comm_OTUs)
# # Confirm the dimensions of your comm file
# dim(comm_OTUs) # 41 samples
# # Change sample names (has . instead of _)
# rownames(comm_OTUs) <- gsub("[.]","_",rownames(comm_OTUs))

# # make a taxonomy-level community matrix
# # species level
# comm_species <- aggregate(t(comm_OTUs), by=list(species=blast_KNP_v6_tax$species), sum)

# length(names(comm_species))
# length(comm_species[names(comm_species)%in%labNotes$Sample_code])
# names(comm_species)[!names(comm_species)%in%labNotes$Sample_code]
# # Only missing DLP_8B_S and DLP_8B_XS
# # Not yet sequenced

### megablast_refLib_v7_noOTU_pident97_length_100
seqs <- read.csv("../Data/Molecular/eDNA_camtrap/classifications/megablast_refLib_v7/seqXspec_megablast_KNP_v7_noOTU_pident97_length100.csv", as.is=T)
seqs$sample <- gsub("[.]","_", seqs$sample)
seqs <- unique(seqs[,c("sample","species")])

seqs$site <- sapply(strsplit(seqs$sample, "_"), function(x) x[1])
seqs$number <- gsub("[A-Z]","",sapply(strsplit(seqs$sample, "_"), function(x) x[2]))
seqs$AB <- gsub("[0-9]","",sapply(strsplit(seqs$sample, "_"), function(x) x[2])) 
seqs$S_XS <- sapply(strsplit(seqs$sample, "_"), function(x) x[3])
seqs$sample_num <- paste(seqs$site, seqs$number, sep="_")

# Making labNotes sample notation match seqs
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
sample_dates <- unique(ddply(sample_dates, ~ sample_num, transform, collected = 
	mean(collected)))

# Changing trap site name to match seqs site name
lookup <- setNames(sort(unique(seqs$site)), sort(unique(trap$site)))
trap$site <- lookup[trap$site]

# Changing env site name to match seqs site name
env$site <- lookup[env$site]
# Subsetting to only include sites for camtrap analysis
env <- env[!is.na(env$site),]

# Joining date collected with seqs
seqs <- join(seqs, sample_dates)

# seen <- NULL
# for (i in seq_along(seqs$sample)){

# 	# Subset by site
# 	trap_sub <- trap[trap$site==seqs$site[i],]
# 	# Subset by time
# 	trap_sub <- trap_sub[trap_sub$timestamp %within% interval(seqs$collected[i] %m+%days(-7), seqs$collected[i]),]
# 	n_waterhole_7 <- ddply(trap_sub, .(site,binomial), summarize, n_waterhole_7 = sum(n_waterhole))	
# 	if (nrow(n_waterhole_7)>0) n_waterhole_7 <- cbind(seqs$sample_num[i], n_waterhole_7)
# 	seen <- rbind(seen, n_waterhole_7)
# }
# View(seen)

# ABOVE WORKS!


# NOW RE-DO DO THAT WE CAN GET INDIVIDUALS SEEN FOR 1:7 DAYS PRIOR

seen_waterhole <- NULL

for (i in seq_along(seqs$sample)){

	# Subset by site
	trap_sub <- trap[trap$site==seqs$site[i],]
	
	for (j in 1:7){
		trap_sub_day <- trap_sub[trap_sub$timestamp %within% interval(seqs$collected[i] %m+%days(-j), seqs$collected[i]),]
		
		if(j==1) {n_waterhole <- ddply(trap_sub_day, .(site,binomial), function(x) {sum(x$n_waterhole)})
			names(n_waterhole)[3] <- paste("n_waterhole_1")}
		
		if(j>1){ 
			n_waterhole_temp <- ddply(trap_sub_day, .(site,binomial), function(x) {sum(x$n_waterhole)})
			names(n_waterhole_temp)[3] <- paste("n_waterhole",j,sep="_")
			n_waterhole <- join(n_waterhole, n_waterhole_temp)
		}
	}		

	if (nrow(n_waterhole)>0) n_waterhole <- cbind(seqs$sample_num[i], n_waterhole)
	seen_waterhole <- rbind(seen_waterhole, n_waterhole)
}

names(seen_waterhole)[1] <- "sample_num"
seen_waterhole <- unique(seen_waterhole)
seen_waterhole <- seen_waterhole[!is.na(seen_waterhole$binomial),]
seen_waterhole <- seen_waterhole[seen_waterhole$binomial!="Homo_sapiens",]
View(seen_waterhole)



# SEEN IN CONTACT

seen_contact <- NULL

for (i in seq_along(seqs$sample)){

	# Subset by site
	trap_sub <- trap[trap$site==seqs$site[i],]
	
	for (j in 1:7){
		trap_sub_day <- trap_sub[trap_sub$timestamp %within% interval(seqs$collected[i] %m+%days(-j), seqs$collected[i]),]
		
		if(j==1) {n_contact <- ddply(trap_sub_day, .(site,binomial), function(x) {sum(x$n_contact)})
			names(n_contact)[3] <- paste("n_contact_1")}
		
		if(j>1){ 
			n_contact_temp <- ddply(trap_sub_day, .(site,binomial), function(x) {sum(x$n_contact)})
			names(n_contact_temp)[3] <- paste("n_contact",j,sep="_")
			n_contact <- join(n_contact, n_contact_temp)
		}
	}		

	if (nrow(n_contact)>0) n_contact <- cbind(seqs$sample_num[i], n_contact)
	seen_contact <- rbind(seen_contact, n_contact)
}

names(seen_contact)[1] <- "sample_num"
seen_contact <- unique(seen_contact)
seen_contact <- seen_contact[!is.na(seen_contact$binomial),]
seen_contact <- seen_contact[seen_contact$binomial!="Homo_sapiens",]
View(seen_contact)


# SEEN PRESENT

seen_present <- NULL

for (i in seq_along(seqs$sample)){

	# Subset by site
	trap_sub <- trap[trap$site==seqs$site[i],]
	
	for (j in 1:7){
		trap_sub_day <- trap_sub[trap_sub$timestamp %within% interval(seqs$collected[i] %m+%days(-j), seqs$collected[i]),]
		
		if(j==1) {n_present <- ddply(trap_sub_day, .(site,binomial), function(x) {sum(x$n_present)})
			names(n_present)[3] <- paste("n_present_1")}
		
		if(j>1){ 
			n_present_temp <- ddply(trap_sub_day, .(site,binomial), function(x) {sum(x$n_present)})
			names(n_present_temp)[3] <- paste("n_present",j,sep="_")
			n_present <- join(n_present, n_present_temp)
		}
	}		

	if (nrow(n_present)>0) n_present <- cbind(seqs$sample_num[i], n_present)
	seen_present <- rbind(seen_present, n_present)
}

names(seen_present)[1] <- "sample_num"
seen_present <- unique(seen_present)
seen_present <- seen_present[!is.na(seen_present$binomial),]
seen_present <- seen_present[seen_present$binomial!="Homo_sapiens",]
View(seen_present)



# Getting amount of time drinking (thereabouts)

str(trap)
trap_sub <- trap[trap$site==seqs$site[1],]
trap_sub_spec <- trap_sub[trap_sub$binomial==seqs$species[1],]
# this causes NAs...
trap_sub_spec <- trap_sub_spec[!is.na(trap_sub_spec$site),]
dim(trap_sub_spec)
head(trap_sub_spec$timestamp)
head(trap_sub_spec$filename)

View(trap_sub_spec)



# Joining 
names(seqs)[names(seqs)=="species"] <- "binomial"
names(seqs)
seqs$edna <- 1
test <- join(seqs[,c("binomial","sample_num","edna")], seen_waterhole, type="full")
test <- unique(test[with(test, order(sample_num, edna, -n_waterhole_1)),])
View(test)
View(test[test$sample_num=="NYA_4",])

View(seqs)
sort(unique(seqs$binomial))









# SCRAPS

# Get date +/- 
# ex:
ymd("2015-07-01") %m+% weeks(-1)

sample_dates$seven <- sample_dates$collected %m+% days(-7)

str(sample_dates)

# NEED TO MATCH TRAP SITE+DATE(S) TO SAMPLE IDS

str(labNotes)
str(trap)
unique(trap$date_range)
unique(labNotes$date)

# create new "sample date" column pasting in the date and time sample was taken
# then turn this into POSIX
# then match each photo in trap based on site and date
# if it falls within the sample date minus one week, mark it...

# might be better way to have this be reproducible and adjustable for different values of cutoff (not one week only)

labNotes$sample_time <- as.POSIXct(with(labNotes, paste("2015", date_collected, time_collected, sep=":")), format="%Y:%B_%d:%H:%M", tz="UTC") 
str(labNotes)



# Creating file with date ranges to subset 
sample_dates <- subset(labNotes, select=c("site","sample_code","date_collected","time_collected"))
sample_dates$collected <- as.POSIXct(with(sample_dates, paste("2015", date_collected, time_collected, sep=":")), format="%Y:%B_%d:%H:%M", tz="UTC")
sample_dates <- subset(sample_dates, select =-c(date_collected, time_collected))
str(sample_dates)

# Get date +/- 
# ex:
ymd("2015-07-01") %m+% weeks(-1)

sample_dates$seven <- sample_dates$collected %m+% days(-7)

str(sample_dates)

# subset trap based on time and site
test <- subset(trap, trap$timestamp %within% interval(labNotes$sample_time[1] %m+%days(-7), labNotes$sample_time[1]))
test <- trap[trap$site=="Imbali" & trap$timestamp %within% interval(sample_dates$collected[sample_dates$site=="Imbali"] %m+%days(-7), sample_dates$collected[sample_dates$site=="Imbali"]),]

trap[trap$timestamp %within% interval(labNotes$sample_time[1] %m+%days(-7), labNotes$sample_time[1]),]


str(test)



# THIS WILL WORK FOR SPECIES LEVEL AGGREGATIONS
# BUT NEED TO MAKE SURE the correct time restricion is in place
test <- ddply(trap, c("site","binomial"), function(x) {sum_n_waterhole <- 
	sum(x$n_waterhole)}
	)


# TRYING TO WRITE FUNCTION TO SUBSET BY time
# WORKS SORT OF BUT CAN"T YET MATCH STARTCOL
split.by.time <- function(dat,timecol,startcol,ndays){
	df.sub <- subset(dat, dat$timecol %within% interval(startcol %m+%days(-ndays), startcol))
}

test <- split.by.time(trap, )

trap.sub <- subset(trap, trap$timestamp %within% interval(sample_dates$collected[1] %m+%days(-7), sample_dates$collected[1]))


head(sort(unique(trap.sub$timestamp)))
tail(sort(unique(trap.sub$timestamp)))



# Maybe best thing to do is temporarily place sample code in trap df
# based on time interval
	
sampleDF <- sample_dates
dfout <- trap
dfout$sample_code <- sampleDF$sample_code[dfout$site%in%sampleDF$site & dfout$timestamp %within% 
		interval(sampleDF$collected[sampleDF$site%in%dfout$site] %m+%days(-7), sampleDF$collected[sampleDF$site%in%dfout$site])]



# NEW STRATEGY:
# for each sample_id, subset trap based on site AND time
# then aggregate based on specie and place this in a list/df





sapply(trap$timestamp, function(x) x %within% interval(sample_dates$collected%m+%days(-7), sample_dates$collected))

for (i in seq_along(trap$timestamp)) {

	site <- trap$site[i]%in%sample_dates$site
	start <- sample_dates$collected[sample_dates$site==site]

	test <- sample_dates$sample_code[ trap$site==site & 
			trap$timestamp[i] %within% interval(start, end) ] 
}



test <- cut(trap$timestamp, by="days")
test <- as.ts(trap)

plot(trap$n_waterhole ~ trap$timestamp)

test <- with(trap, (tapply(n_waterhole, cut(timestamp, breaks="hour"), sum)))

with(trap, aggregate(n_waterhole, by=list(binomial), sum))
with(trap, aggregate(n_waterhole, by=list(cut(timestamp, "hour")), sum))


test <- tapply(trap$n_waterhole, cut(trap$timestamp, breaks="10 min"), sum)
test <-cbind(names(test), test)

test <- aggregate(trap$n_waterhole, by=list(cut(trap$timestamp, breaks="hour")), sum)
str(test)

plot(test)

cut.POSIXt(trap$timestamp, breaks="hour")
class(trap$timestamp)


z <- aggregate(n_waterhole ~ timestamp, 
	seq(from=as.Date("2015-06-01"), to=as.Date("2015-07-01"), by="day"),
	data=trap,FUN=sum)

View(z)

june <- trap$timestamp["2015-06-01"<=trap$timestamp<="2015-07-01"]
june <- subset(trap, timestamp>="2015-06-01" & timestamp<"2015-07-01")
sort(unique(june$timestamp))
# didn't quite work - stopped on july 1st at 6 am...
# after setting to UTC didn't quite work - stopped on july 1st at 4 am...

require(lubridate)

june <- trap[trap$timestamp %within% interval(ymd("2015-06-01"),ymd_hms("2015-06-30 23:59:59")),]
sort(unique(june$timestamp))


july_1_minus_week <- trap[trap$timestamp %within% interval(ymd("2015-07-01"),ymd_hms("2015-06-30 23:59:59")),]



# still off - stopped on july 1st at 2 am...
# must be a timezone issue
# works now that trap$timestamp was set to UTC
tz(ymd("2015-06-01"))

tz(trap$timestamp)

require(ggplot2)
ggplot( data = trap, aes(timestamp, n_waterhole, group=binomial)) + geom_line() 

str(trap)
str(env)
table(env$site, env$date)
str(labNotes)
table(labNotes$site, labNotes$date_collected)


labNotes <- subset(labNotes, select=c(site,sample_code,date_collected,time_collected,sampling_direction,wind_direction))
names(labNotes)[names(labNotes)=="date_collected"] <-"date"
str(labNotes)
str(env)

test <- join(env, labNotes[,c("site","date","sample_code")], type="left")
View(test)


# Testing trap aggregation

test <- table(trap$binomial)
test <- aggregate(trap$n_waterhole, by=list(trap$binomial, trap$site, trap$date_range), sum)
View(test)
str(test)
table(test)



# NEED TO MATCH TRAP SITE+DATE(S) TO SAMPLE IDS

str(labNotes)
str(trap)
unique(trap$date_range)
unique(labNotes$date)

# create new "sample date" column pasting in the date and time sample was taken
# then turn this into POSIX
# then match each photo in trap based on site and date
# if it falls within the sample date minus one week, mark it...

# might be better way to have this be reproducible and adjustable for different values of cutoff (not one week only)