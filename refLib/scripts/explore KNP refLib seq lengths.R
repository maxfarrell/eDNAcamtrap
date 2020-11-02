# explore KNP refLib seq lengths
# in mamm_lib

# dat <- read.table("KNP_refLib_seq_lengths.txt", as.is=T)
dat <- read.table("../output/KNP_seq_lengths.txt", as.is=T, sep="\t")
str(dat)
range(dat[,2])
names(dat) <- c("processid","length")

head(unique(dat[,1]))
tail(unique(dat[,1]))

BOLD <- dat[grep("-",dat[,1]),]
GB <- dat[grep("-",dat[,1], invert=TRUE),]

range(dat$length)#100 - 2868
range(GB$length)#100 -2868 

dim(dat)#8249
dim(GB)#4764

dim(dat[dat$length>600,])#7105
dim(BOLD[BOLD[,2]>600,])#3197
dim(GB[GB[,2]>600,])#3908


dim(dat[dat[,2]<500,])#689
dim(BOLD[BOLD[,2]<500,])#245
dim(GB[GB[,2]<500,])#444

# Identifying which species have long / short sequences
# tax <- read.delim("../KNP_refTaxonomy_v3.tsv", header=FALSE, sep=";", stringsAsFactors=FALSE)
tax <- read.delim("../output/KNP_refTaxonomy.tsv", header=FALSE, sep=";", stringsAsFactors=FALSE)
id_split <- data.frame(tax[,1], do.call(rbind, strsplit(tax[,1], "[\t]")), stringsAsFactors=FALSE)
tax <- cbind(id_split[,-1],tax[,-1])
names(tax) <- c("processid","phylum", "class","order","family","genus","species")
length(unique(tax$species))#410


# # length(unique(tax$order))#11

# spec_order <- table(unique(tax[,c("species","order")]))
# colSums(spec_order)
# 15 Carnivora
# 23 Cetartiodactyla
# 20 Chiroptera
# 1 Hyracoidea
# 1 Lagomorpha
# 3 Perrisodactyla
# 1 Pholifota
# 5 Primates
# 7 Rodentia
# 1 Tubulidentata

require(plyr)

str(tax)
str(dat)
names(dat) <- c("processid","length")
dat <- join(dat,tax, type="left")

table(dat$species)
dat$length[dat$species=="Panthera_pardus"]

dat$species[dat$length==max(dat$length)]
max(dat$length)
dat$processid[dat$length==max(dat$length)]
dat$species[dat$length==max(dat$length)]
# Cuculus_canorus has longest COI seq of 2868
# check processid AY274060.1
# This is the COI gene for C. canorus, 
# but COI should be about 1500 bp, so not sure what extra seq is...


sort(unique(dat$length))
hist(dat$length)

min_length <- 250
table(dat$species[dat$length < min_length])
sort(unique(dat$species[dat$length < min_length]))
short_spec <- unique(dat$species[dat$length < min_length])
sum(table(dat$species[dat$species%in%short_spec])==1)
# Would lose 8 species by cutting off below min_length
sum(table(dat$species[dat$species%in%short_spec]))
table(dat$species[dat$species%in%short_spec])

table(dat$species[dat$length > 1500])
dat$species[dat$length > 1500]
short_spec <- unique(dat$species[dat$length > 1500])
table(dat$species[dat$species%in%short_spec])

table(dat$species)







