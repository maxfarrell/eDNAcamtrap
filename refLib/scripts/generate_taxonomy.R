# Formatting GB_BOLD reference library from R

# Maxwell Farrell
# maxwellfarrell@gmail.com

# in shell make sure to cat all Mammal_FASTAs AND mito_genomes *fastas into one file before
# cat * > Mammal_FASTAs_cat.fasta 
# in mito_genomes: cat *.fasta > mito_CO1s.fasta
# cat both sets (mv first)

# Final FASTA
fas <- read.table("../output/FASTAs_cat.fasta", as.is=T, sep="\n")
headers <- fas[grep(">", fas[,1]),]

# KNP species taxonomy
# Species Lists
mammList <- read.csv("../data/species_lists/kruger_mamm.csv", as.is=T, header=T)
mammList <- cbind("Mammalia", mammList)
names(mammList)[1] <- "class"

amphList <- read.csv("../data/species_lists/kruger_amph.csv", as.is=T, header=T)
amphList <- cbind("Amphibia", amphList)
names(amphList)[1] <- "class"

birdList <- read.csv("../data/species_lists/kruger_birds.csv", as.is=T, header=T)
birdList <- cbind("Aves", birdList)
names(birdList)[1] <- "class"

reptList <- read.csv("../data/species_lists/kruger_rept.csv", as.is=T, header=T)
reptList <- cbind("Reptilia", reptList)
names(reptList)[1] <- "class"

fishList <- read.csv("../data/species_lists/kruger_fish.csv", as.is=T, header=T)

KNPtax <- rbind(mammList, amphList, birdList, reptList, fishList)
write.csv(KNPtax, "../output/KNP_inventory_taxonomy.csv", row.names=F)

# length(unique(KNPtax$binomial))# 854
# length(unique(c(KNPtax$binomial,KNPtax$synonym)))# 1,283

tax <- unique(as.data.frame(cbind(KNPtax$binomial, with(KNPtax, noquote(paste(
	   "Chordata",class, order, family, genus, binomial, 
	   sep=";")))), stringsAsFactors=FALSE))
names(tax) <- c("binomial","taxonomy")
tax[,2] <- gsub(" ","_", tax[,2])


# Matching headers to binomials
refLib <- as.data.frame(as.character(NULL), stringsAsFactors=FALSE)

for (i in seq_along(KNPtax$binomial)){

	a <- headers[grep(KNPtax$binomial[i], headers)]
	if (length(a)>0) refLib <- rbind(refLib, cbind(rep(KNPtax$binomial[i],length(a)),a))
}

# str(refLib)# 8254 observations

# Adding synonyms
for (i in seq_along(KNPtax$synonym)){

	a <- headers[grep(KNPtax$synonym[i], headers)]
	if (length(a)>0) refLib <- rbind(refLib, cbind(rep(KNPtax$synonym[i],length(a)),a))
}
refLib <- unique(refLib)
names(refLib) <- c("binomial","header")
refLib$binomial <- as.character(refLib$binomial)
# str(refLib)#8430 obs.

# length(unique(refLib$binomial))#432
# length(unique(tax$binomial))#854
# length(intersect(tax$binomial, refLib$binomial))#383
# length(intersect(KNPtax$synonym, refLib$binomial))#50

# Some binomials in the refLib are actually synonyms.
# To join with tax, these synonyms in refLib$binomial file 
# need to be replaced with the KNPtax$binomial variants

synonyms <- setNames(KNPtax$binomial, KNPtax$synonym)
refLib$binomial[refLib$binomial%in%KNPtax$synonym] <- synonyms[refLib$binomial[refLib$binomial%in%KNPtax$synonym]]

# length(unique(refLib$binomial))#411
# length(unique(tax$binomial))#854
# length(intersect(tax$binomial, refLib$binomial))#410
# length(intersect(KNPtax$synonym, refLib$binomial))#1
# intersect(KNPtax$synonym, refLib$binomial)# - this is NA

# Joining refLib with tax
require(dplyr)
refLib <- left_join(refLib, tax)
refLib <- subset(refLib, select=-c(binomial))

# Removing any headers or taxonomy strings that are NA
refLib <- refLib[!is.na(refLib$header),]
refLib <- refLib[!is.na(refLib$taxonomy),]

# Cleaning header names
refLib$header <- sub("(.*?[|].*?)[|].*", "\\1", refLib$header)
refLib$header <- sub("(.*?) .*", "\\1", refLib$header)
refLib$header <- gsub("gi[|]", "", refLib$header)
refLib$header <- gsub(">", "", refLib$header)

write.table(refLib, "../output/KNP_refTaxonomy.tsv", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
