#!/bin/bash
# Grab BOLD and GB references based on species KNP species lists
# cat individual fastas
# generate taxonomy mapping file from mammal list and fasta
# format headers in fasta

# Maxwell Farrell 
# maxwellfarrell@gmail.com

# first grab sequences
# Rscript --vanilla GB_BOLD_seq_download.R
# Rscript --vanilla CO1_from_GB_mito_genomes.R

# cat individual fastas per species into one refLib
cat ../data/mito_genomes/*.fasta ../data/FASTAs/* > ../output/FASTAs_cat.fasta

# REMOVE SEQUENCES WITH AMBIGUOUS BASES (N or -)
# this does not remove other ambig bases (MKYSR)

# Change multi-line sequences to single-line
gawk -i inplace '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ../output/FASTAs_cat.fasta

# Identify all sequences that contain trailing n or N or - one or more times
# grep -B 1 '[nN\-]\+$' ../output/FASTAs_cat.fasta
# Identify all sequences that contain leading n or N or - one or more times
# grep -B 1 '^[nN\-]\+' ../output/FASTAs_cat.fasta
 
# Remove trailing ambiguous bases
sed -i 's/[nN\-]\+$//g' ../output/FASTAs_cat.fasta 
# Remove leading ambiguous bases
sed -i 's/^[nN\-]\+//g' ../output/FASTAs_cat.fasta 

# Remove sequences that have ambiguous bases and the preceding header line
grep -B 1 "^[ACTGactg].*[nN\-]" ../output/FASTAs_cat.fasta | grep -v -f - ../output/FASTAs_cat.fasta > ../output/FASTAs_cat_temp.fasta
mv -f ../output/FASTAs_cat_temp.fasta ../output/FASTAs_cat.fasta

# Verify that there are no sequences with ambig bases
grep -B 1 '^[ACTGactg].*[nN\-]' ../output/FASTAs_cat.fasta


# Generate taxonomy file ("KNP_refTaxonomy.tsv") if it doesnt exist
if [ ! -f ../output/KNP_refTaxonomy.tsv ]; then
	Rscript "generate_taxonomy.R"
fi

# Trim taxonomy to remove extra trailing ; - comment out this line if using mothur for tax assignment
sed 's/[;_]*$//' -i ../output/KNP_refTaxonomy.tsv

# Remove eveything after second instance of "|" in FASTA header (returns gb number)
cut -d '|' -f1,2 ../output/FASTAs_cat.fasta > ../output/temp.fasta

# Remove eveything before first instance of "|" in FASTA header
sed 's/gi|//g' -i ../output/temp.fasta 

# Remove blank lines
sed '/^$/d' -i ../output/temp.fasta
sed '/^ $/d' -i ../output/temp.fasta

# Remove everything after the first " " in FASTA header
cut -d ' ' -f 1 ../output/temp.fasta > ../output/FASTAs_cat.fasta

# Remove temp.fasta
rm ../output/temp.fasta

# Change multi-line sequences to single-line
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ../output/FASTAs_cat.fasta > ../output/FASTAs_cat_singleline.fasta

rm ../output/FASTAs_cat.fasta

# Prune FASTA sequences to only those in the taxonomy
# NEED TO SUBSET processid from specData
cut -f1 ../output/KNP_refTaxonomy.tsv > ../output/taxnomy_ids.txt
grep -Fwf ../output/taxnomy_ids.txt -A 1  ../output/FASTAs_cat_singleline.fasta | grep -v '^--$'  > ../output/FASTAs_taxMatch.fasta
# The grep -v '^--$' simply filters out the lines with -- that grep adds between groups of output lines when using the -A option.

rm ../output/FASTAs_cat_singleline.fasta

# The following removes any duplicated sequence IDs in the fasta (they seem to come from GB)
sed -e '/^>/s/$/@/' -e 's/^>/#/' ../output/FASTAs_taxMatch.fasta | tr -d '\n' | tr "#" "\n" | tr "@" "\t" | sort -u -t ' ' -f -k1,1 | sed -e 's/^/>/' -e 's/\t/\n/' > ../output/FASTAs_final.fasta

# For some reason the first line is a lone ">"
# Remove lines that are only ">"
sed '/^>$/d' -i ../output/FASTAs_final.fasta

rm ../output/FASTAs_taxMatch.fasta

# Length distributions of FASTAs
cat ../output/FASTAs_final.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > ../output/KNP_seq_lengths.txt


##
# Attempting to clean based on alginment to reference

# muscle alginment for huge alignments
muscle -in ../output/FASTAs_final.fasta -out ../output/FASTAs_final_align.afa -maxiters 1 -diags1 -sv

# Add goat reference to alignment to identify folmer region
muscle -profile -in1 ../output/FASTAs_final_align.afa -in2 ../data/Capra_hircus_GBMA3852-12_658bp_ref.fas -out ../output/refLib_plusGoat.afa

# Alignment manually viewed with AliView and sequences trimmed to only include the range of the reference sequence
# output: refLib_plusGoat_trimmed.afa

# Removing duplicated sequences (different headers but same sequence)
# awk '/>/{s=$0;next}!a[$0]++{print s;print}' example_duplicated_seqs.fasta
awk '/>/{s=$0;next}!a[$0]++{print s;print}' ../output/refLib_plusGoat_trimmed.afa > ../output/refLib_plusGoat_trimmed_uniq.afa

# Re-alignment with muscle
muscle -in ../output/refLib_plusGoat_trimmed_uniq.afa -out ../output/refLib_plusGoat_trimmed_uniq_aligned.afa -maxiters 1 -diags1 -sv

# Removing poorly aligning sequences manually in AliView
# AY274060.1 - Cuculus canorus
# JN543896.1 - Pachydactylus punctatus
# KY177138.1 - Kassina senegalensis
# KY177150.1 - Hyperolius pusillus
# KY177122.1 - Pyxicephalus edulis
# KF665515.1 - Schismaderma carens
# KY177086.1 - Hemisus marmoratus (genbank title even says UNVERIFIED)
# EU349707.1 - Dasymys incomtus
# AF028191.1 - Canis mesomelas elongae
# GBMIN41767-14 - Canis mesomelas elongae
# GBIR7653-19 - Coturnix coturnix (in BOLD, but accession # says UNVERIFIED in genbank)
# JX178002.1 - Coturnix coturnix UNVERIFIED
# KY765588.1 - Anguilla marmota
# AF028190.1 - Canis mesomelas elongae
# GBMIN41766-14 - Canis mesomelas elongae
# AY030184.1 - Passer domesticus

# End of cleaning round 1

# Re-alignment with muscle
muscle -in ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned.afa -out ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned.afa -maxiters 1 -diags1 -sv

# Removing poorly aligning sequences manually in AliView
# KY177139.1 - UNVERIFIED: Kassina maculata
# DQ249073.1 - Acontias plumbeus
# KY177165.1 - Arthroleptis stenodactylus
# KT728354.1 - UNVERIFIED: Anguilla marmorata
# EU541346.1 - Acanthopagrus berda
# AF028184.1 - Canis adustus
# GBMIN41760-14 - Canis adustus
# KY177164.1 - Arthroleptis stenodactylus
# AF407487.1 - Eurystomus glaucurus
# KJ179594.1 - Nothobranchius orthonotus

# Trimmed sequences to length of goat folmer region

# Removed alignment gaps
sed '/^[^>]/s/[-]//g' refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed.afa > refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap.afa

# Check length
cat ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap.afa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > KNP_seq_lengths.txt

# We will lose 1009 sequences and representation of 9 species by trimming to sequences > 250
# Trimming to sequence length >250 bp

# Minimum length trimming
awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 250 {print ">"$0}' ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap.afa > refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap_noshort.afa

# Double check lengths - worked
# cat ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap_noshort.afa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'

# Re-alignment with muscle
muscle -in ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap_noshort.afa -out ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap_noshort_aligned.afa -maxiters 1 -diags1 -sv

# Removing poorly aligning sequences manually in AliView
# MH176327.1 - Hypophthalmichthys molitrix
# HNDS001-18 - Anguilla marmorata
# KX657711.1 - Cyprinus carpio
# KY938000.1 - Cyprinus carpio
# EU621620.1 - Falco peregrinus
# GBIR7654-19 - UNVERIFIED: Coturnix coturnix
# JX178044.1 - UNVERIFIED: Coturnix coturnix

# Removing duplicated sequences (different headers but same sequence)
# awk '/>/{s=$0;next}!a[$0]++{print s;print}' example_duplicated_seqs.fasta
awk '/>/{s=$0;next}!a[$0]++{print s;print}' ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap_noshort_aligned_cleaned_nogap.afa > ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap_noshort_aligned_cleaned_nogap_uniq.afa

# Re-alignment with muscle
muscle -in ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap_noshort_aligned_cleaned_nogap_uniq.afa -out ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap_noshort_aligned_cleaned_nogap_uniq_aligned.afa -maxiters 1 -diags1 -sv

# Removing poorly aligning sequences manually in AliView
# AF279732.1 - Columba guinea
# KT023333.1 - Oena capensis

# Removed alignment gaps
sed '/^[^>]/s/[-]//g' ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap_noshort_aligned_cleaned_nogap_uniq_aligned_cleaned.afa > ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap_noshort_aligned_cleaned_nogap_uniq_aligned_cleaned_nogap.afa

# Saving this version as "final_Oct15_2019"
cp ../output/refLib_plusGoat_trimmed_uniq_aligned_cleaned_aligned_cleaned_trimmed_nogap_noshort_aligned_cleaned_nogap_uniq_aligned_cleaned_nogap.afa ../output/refLib_final_Oct15_2019.fasta


