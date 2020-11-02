# KNP Mammal COI reference library diagnostics
# Plot sequence alignment
# Create and plot neighbour joining tree

# https://github.com/sjspielman/alignfigR
# http://sjspielman.org/alignfigR/
require(alignfigR)

lib <- read_alignment("../Mammal_FASTAs_final_v3.trim.good.good.align.muscle.fasta")
seq_lengths <- sapply(lib, length)
range(seq_lengths)

is._ <- function(x){
	return(sum(x=="-"))
}
# seq_NAs <- sapply(lib, is._)
# hist(seq_NAs)

tax <- read.table("../KNP_refTaxonomy_v3.tsv", sep="\t", as.is=T)
tax_list <- strsplit(tax[,2], "[;]")
maxlevels <- 6
taxlevels <- as.data.frame(t(sapply(tax_list,'[',1:maxlevels)), stringsAsFactors=FALSE)
tax <- data.frame(tax,taxlevels)
tax <- tax[,-2]
names(tax) <- c("sseqid","phylum", "class","order","family","genus","species")

# Order taxonomy by order, family, genus, species
tax <- tax[with(tax, order(order, family, genus, species)),]

# There are some sequences in tax that are not in refLib
tax[!tax$sseqid%in%names(lib),]
dim(tax[!tax$sseqid%in%names(lib),])# 166...
# not sure why this is occuring... 
# maybe cleaning the references removed some and were not trimmed from tax
tax <- tax[tax$sseqid%in%names(lib),]

# Remove volant mammals (Chiroptera) from tax and lib
tax <- tax[tax$order!="CHIROPTERA",]
nrow(tax)# 779 sequences
length(unique(tax$species))# 57

# Order lib by taxonomy
lib <- lib[order(match(names(lib), tax$sseqid))]

# Swap sseqiq for binomial + sseqid
lookup <- setNames(paste(tax$species, tax$sseqid, sep="_"), tax$sseqid)
names(lib) <- lookup[names(lib)]

# Remove sequences with name of NA (Chiroptera) 
lib <- lib[!is.na(names(lib))]

require(seqinr)
# write.fasta(sequences=lib,names=names(lib),file.out="../Mammal_FASTAs_final_v3.trim.good.good.align.muscle.nonVolant.taxOrdered.fasta")

# Alignment of Capra hircus 658bp full length COI barcode with AliView
# Trim of alignment to include only this region
# Deletion of gap-only columns
# Output = Mammal_FASTAs_final_v3.trim.good.good.align.muscle.nonVolant.taxOrdered.folmer.fasta 

lib <- read.fasta("../Mammal_FASTAs_final_v3.trim.good.good.align.muscle.nonVolant.taxOrdered.folmer.fasta")
# Questionable sequences to remove:
# (Viewing sequences in AliView)
# Canis_adustus_GBMIN41760-14
# Canis_adustus_GBMIN43481-14
# Canis_adustus_695275528
# Canis_adustus_2826546
# Canis_mesomelas_2826560
# Canis_mesomelas_2826558
# Lycaon_pictus_GBMIN41771-14
# Lycaon_pictus_2826568
# Otocyon_megalotis_2826574 # Note this is the only sequence for this species (bat-eared fox)
# KEEP ^
# Diceros_bicornis_530539552 # Short sequence
# Diceros_bicornis_530539550 # Short sequence
# Thryonomys_swinderianus_695275928

seqs_to_remove <- c("Canis_adustus_GBMIN41760-14","Canis_adustus_GBMIN43481-14",
	"Canis_adustus_695275528","Canis_adustus_2826546","Canis_mesomelas_2826560",
	"Canis_mesomelas_2826558","Lycaon_pictus_GBMIN41771-14","Lycaon_pictus_2826568",
	"Diceros_bicornis_530539552","Diceros_bicornis_530539550","Thryonomys_swinderianus_695275928")

lib <- lib[!names(lib)%in%seqs_to_remove]
# write.fasta(sequences=lib,names=names(lib),file.out="../Kruger_Mammals_refLib.fasta")

### Formatting for dada2
# https://benjjneb.github.io/dada2/training.html

# Example of dada2 formatting for Kingdom to Genus reference database (6 levels)
# >Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Candidatus_Regiella;
# TTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGGCAGCGGGGAGTAGCTTGCTACT [...]

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
# write.fasta(sequences=lib,names=names(lib),file.out="../Kruger_Mammals_refLib_dada2.fasta")

# Example of species level assignment
# >ID Genus species
# ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGCTGTAGCATGA [...]
names(lib) <- with(dada2_tax, paste(sseqid,species,sep=" "))
names(lib) <- gsub("_"," ",names(lib))
# write.fasta(sequences=lib,names=names(lib),file.out="../Kruger_Mammals_refLib_dada2_species.fasta")

# Trying with 6 levels, but from phylum to species
# Reformatting names
dada2_tax <- tax
dada2_tax$species_sseqid <- paste(dada2_tax$species, dada2_tax$sseqid, sep="_")
dada2_tax <- dada2_tax[dada2_tax$species_sseqid%in%names(lib),]
names(lib_original) <- with(dada2_tax, paste(phylum,class,order,family,genus,species,sep=";"))

# Reformatting sequences: remove "-" & convert to UPPER CASE GCTA
reFASTA <- function(x){
	toupper(x[x!="-"])
}

lib_original <- lapply(lib_original, reFASTA)
# write.fasta(sequences=lib_original,names=names(lib_original),file.out="../Kruger_Mammals_refLib_dada2_phy2species.fasta")






######################
## NJ Tree building ##
######################
require(ape)
require(adegenet)

# tree for lib after removing volant mammals and ordering by taxonomy
bin <- fasta2DNAbin("../Mammal_FASTAs_final_v3.trim.good.good.align.muscle.nonVolant.taxOrdered.fasta", chunkSize=10)
tree <- njs(dist.dna(bin, model="K80"))
# write.tree(tree,"../KNP_refLib_nonvolant_nj.tree")

# tree for lib after trimming to folmer region
bin <- fasta2DNAbin("../Mammal_FASTAs_final_v3.trim.good.good.align.muscle.nonVolant.taxOrdered.folmer.fasta", chunkSize=10)
tree <- njs(dist.dna(bin, model="K80"))
# write.tree(tree,"../KNP_refLib_nonvolant_folmer_nj.tree")

# tree for lib after removing misaligned sequences
bin <- fasta2DNAbin("../Kruger_Mammals_refLib.fasta", chunkSize=10)
tree <- njs(dist.dna(bin, model="K80"))
write.tree(tree,"../Kruger_Mammals_refLib_nj.tree")

pdf(height=220, width=72)
plot(tree)
# plot(tree, type="fan")
dev.off()