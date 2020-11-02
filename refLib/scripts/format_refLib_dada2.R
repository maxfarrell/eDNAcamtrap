# Formatting refLib for dada2

require(seqinr);packageVersion("seqinr")# 3.4.5

# Data
lib <- read.fasta("../output/refLib_final_Oct15_2019.fasta")
tax <- read.table("../output/KNP_refTaxonomy.tsv", sep="\t", as.is=T)

# 
lib_original <- lib


## Formatting for dada2
# https://benjjneb.github.io/dada2/training.html

# Example of dada2 formatting for Kingdom to Genus reference database (6 levels)
# >Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Candidatus_Regiella;
# TTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGGCAGCGGGGAGTAGCTTGCTACT [...]# 

tax_list <- strsplit(tax[,2], "[;]")
maxlevels <- 6
taxlevels <- as.data.frame(t(sapply(tax_list,'[',1:maxlevels)), stringsAsFactors=FALSE)
tax <- data.frame(tax,taxlevels)
tax <- tax[,-2]
names(tax) <- c("sseqid","phylum", "class","order","family","genus","species")

# Order taxonomy by order, family, genus, species
tax <- tax[with(tax, order(order, family, genus, species)),]

# There are sequences in tax that are not in refLib (due to cleaning)
tax[!tax$sseqid%in%names(lib),]
dim(tax[!tax$sseqid%in%names(lib),])# 3464...
tax <- tax[tax$sseqid%in%names(lib),]

# Remove volant mammals (Chiroptera) from tax and lib
# tax <- tax[tax$order!="CHIROPTERA",]
# nrow(tax)# 779 sequences
# length(unique(tax$species))# 57

# Order lib by taxonomy
lib <- lib[order(match(names(lib), tax$sseqid))]

# Swap sseqiq for binomial + sseqid
lookup <- setNames(paste(tax$species, tax$sseqid, sep="_"), tax$sseqid)
names(lib) <- lookup[names(lib)]

# Remove sequences with name of NA (if any)
lib <- lib[!is.na(names(lib))]

# Keeping copy of lib
lib_original <- lib

# Reformatting names
dada2_tax <- tax
dada2_tax$species_sseqid <- paste(dada2_tax$species, dada2_tax$sseqid, sep="_")
dada2_tax <- dada2_tax[dada2_tax$species_sseqid%in%names(lib),]
names(lib) <- with(dada2_tax, paste("Eukarya",phylum,class,order,family,genus,sep=";"))

# Reformatting sequences: remove "-" & convert to UPPER CASE GCTA
reFASTA <- function(x){
	toupper(x[x!="-"])
}

lib <- lapply(lib, reFASTA)
write.fasta(sequences=lib,names=names(lib),file.out="../output/Kruger_Vertebrates_refLib_dada2.fasta")

# Example of species level assignment
# >ID Genus species
# ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGCTGTAGCATGA [...]
names(lib) <- with(dada2_tax, paste(sseqid,species,sep=" "))
names(lib) <- gsub("_"," ",names(lib))
write.fasta(sequences=lib,names=names(lib),file.out="../output/Kruger_Vertebrates_refLib_dada2_species.fasta")


# Formatting to 6 levels, but from phylum to species
# Reformatting names
dada2_tax <- tax
dada2_tax$species_sseqid <- paste(dada2_tax$species, dada2_tax$sseqid, sep="_")
dada2_tax <- dada2_tax[dada2_tax$species_sseqid%in%names(lib_original),]
names(lib_original) <- with(dada2_tax, paste(phylum,class,order,family,genus,species,sep=";"))

# Reformatting sequences: remove "-" & convert to UPPER CASE GCTA
reFASTA <- function(x){
	toupper(x[x!="-"])
}

lib_original <- lapply(lib_original, reFASTA)

write.fasta(sequences=lib_original,names=names(lib_original),file.out="../output/Kruger_Vertebrates_refLib_dada2_phy2species.fasta")

# reflib diagnostics

length(unique(dada2_tax$sseqid))# 4965
length(unique(dada2_tax$species))# 391
head(with(dada2_tax, table(class)))

# Number species per class
dada2_tax %>% group_by(class) %>% summarize(nspecies=length(unique(species)))

# Number species per Order
dada2_tax %>% filter(class=="Mammalia" ) %>% group_by(order) %>% summarize(nspecies=length(unique(species)))

knp_mammals <- 
