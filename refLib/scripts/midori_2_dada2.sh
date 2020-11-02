## Converting MIDORI_UNIQUE_20180221_COI.fasta to dada2 format

# merge fasta to singe lines
cat MIDORI_UNIQUE_20180221_COI.fasta | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' > temp.fasta

# add semicolon to end of headers
sed 's/^\(>.*\)$/\1;/' -i temp.fasta

# To greate a species level database
# remove everything after ">" and before second last ";"
sed 's/[^>]root.*;\(.*;\)/ \1/' temp.fasta > MIDORI_UNIQUE_COI_dada2_Species.fasta

## Creating Kingdom -> Genus and Phylum -> Species databases

# remove everything after ">" and before first ";"
sed 's/[^>][^;]*;//' -i temp.fasta
# This results in the FASTA headers starting at Eukarya and going to latin binomial

# To create a six level database that goes from Phylum to Genus_species
# repeat above line once more to go remove Kingdom distinction
sed 's/[^>][^;]*;//' temp.fasta > MIDORI_UNIQUE_COI_dada2_Phylum_Species.fasta

# To create a six level Kingdom to Genus database
# remove the species distinction before the final ";" and replace with ";"
sed 's/[^>][^;]*;$/;/' temp.fasta > MIDORI_UNIQUE_COI_dada2_Kingdom_Genus.fasta

rm temp.fasta