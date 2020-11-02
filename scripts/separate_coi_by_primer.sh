# COI Pipeline - Separated by pcr step / primers

# Separate COI sequences

# ### MiSeq output file structure ###

# Example from Kruger 2015: Max-DLP-8A_S7_L001_R1_001.fastq.gz
# All 2015 samples have: 
# 	- been run on individual lanes (L001)
# 	- have two reads (R1 / R2)
# 	- have one one set number (001) 

#### Pair reads and make contigs with MOTHUR ####

# 1. make.contigs and pull out 16S

# First thing is to combine two sets of reads using make.contigs.
# To do this across all samples, make.contigs can have an input file
# listing the sample name and the files corresponging to each read.
# You can make this with function make.file(inputdir=<folder>)
	
	# mothur	"#make.file(inputdir=../raw_data)"
	# creates file fileList.paired.file)
	# filled in sample names manually
	# Rename and move to folder above (used for 18S and COI):
	# mv fileList.paired.file paired.file

	# Add sample IDs
	# went in manually to add sample name as first column 

	mothur "#make.contigs(file=paired.file, processors=4)"
	# It took 2909 secs to process 11359597 sequences.

	# Group count: 
	# BLANK_2	92452
	# BLANK_3	28916
	# DLP_8A	156722
	# DLP_8A_S	116172
	# DLP_8A_XS	122298
	# DLP_8B	267599
	# DLP_8B_S	126642
	# DLP_8B_XS	83140
	# GIR_1A	115777
	# GIR_1B	110178
	# GIR_2A	104587
	# GIR_2B	109795
	# GIR_3A	110961
	# GIR_3B	116341
	# HOY_2A	148482
	# HOY_2B	217004
	# HOY_3A	60728
	# HOY_3A_S	126805
	# HOY_3A_XS	118097
	# HOY_3B	166235
	# HOY_3B_S	49523
	# HOY_3B_XS	132499
	# HOY_4A	79259
	# HOY_4B	123464
	# IMB_2A	109091
	# IMB_2B	94283
	# IMB_3A	157646
	# IMB_3B	101227
	# IMB_4A	107872
	# IMB_4B	107254
	# KWA_5A	192231
	# KWA_5B	213107
	# KWA_6A	100388
	# KWA_6A_S	124897
	# KWA_6A_XS	124562
	# KWA_6B	143380
	# KWA_6B_S	93689
	# KWA_6B_XS	164370
	# NGO_2A	145781
	# NGO_2B	191240
	# NGO_3A	173089
	# NGO_3A_S	125472
	# NGO_3A_XS	143578
	# NGO_3B	152986
	# NGO_3B_S	127717
	# NGO_3B_XS	138996
	# NGO_4A	103019
	# NGO_4B	99879
	# NHL_2A	161846
	# NHL_2B	92307
	# NHL_3A	91166
	# NHL_3B	91277
	# NHL_4A	70946
	# NHL_4B	105627
	# NWA_2A	231498
	# NWA_2A_S	131982
	# NWA_2A_XS	107384
	# NWA_2B	163278
	# NWA_2B_S	134078
	# NWA_2B_XS	179778
	# NWA_3A	112315
	# NWA_3B	108002
	# NWA_4A	121824
	# NWA_4B	109247
	# NWA_5A	85784
	# NWA_5B	113344
	# NWA_6A	135955
	# NWA_6B	103073
	# NWA_7A	111252
	# NWA_7B	121396
	# NWA_8A	135722
	# NWA_8B	138010
	# NYA_2A	166566
	# NYA_2B	174504
	# NYA_3A	153468
	# NYA_3B	182601
	# NYA_4A	135390
	# NYA_4A_S	120521
	# NYA_4A_XS	135578
	# NYA_4B	185691
	# NYA_4B_S	154440
	# NYA_4B_XS	160274
	# WIT_2A	114485
	# WIT_2B	85673
	# WIT_3A	189929
	# WIT_3B	106032
	# WIT_4A	92051
	# WIT_4B	125873

	# Total of all groups is 11359597

	# Output File Names: 
	# paired.trim.contigs.fasta
	# paired.trim.contigs.qual
	# paired.contigs.report
	# paired.scrap.contigs.fasta
	# paired.scrap.contigs.qual
	# paired.contigs.groups

		# [WARNING]: your sequence names contained ':'.  
		# I changed them to '_' to avoid problems in your 
		# downstream analysis.

	mothur "#summary.seqs(fasta=paired.trim.contigs.fasta, processors=4)"

	#                Start   End     NBases  Ambigs  Polymer NumSeqs
	# Minimum:        1       35      35      0       2       1
	# 2.5%-tile:      1       76      76      0       3       283990
	# 25%-tile:       1       277     277     0       4       2839900
	# Median:         1       353     353     0       5       5679799
	# 75%-tile:       1       410     410     1       5       8519698
	# 97.5%-tile:     1       481     481     14      8       11075608
	# Maximum:        1       602     602     111     300     11359597
	# Mean:   1       332.281 332.281 1.67656 4.89861
	# # of Seqs:      11359597

	# Output File Names: 
	# paired.trim.contigs.summary

	# It took 19 secs to summarize 11359597 sequences.


# 1.2 Including trim with oliogs file... 

	# Using oligos file that sets each F + R pair as a forward primer 
	# Using checkorient=T to make sure the reverse complement of each is checked
	# Using keepforward=F to remove primers

	mothur "#trim.seqs(fasta=paired.trim.contigs.fasta, oligos=/data/kruger/illumina/oligos_forwardFR.txt, keepforward=F, checkorient=T, pdiffs=2, processors=4)"

	# Group count: 
	# 16S_v3v4	2164262
	# 18S_euk	3080303
	# 18S_nem	1129948
	# BR5	2312760
	# mam_ckt	11110
	# mam_ckt_230	2661214
	# Total of all groups is 11359597

	# Output File Names: 
	# paired.trim.contigs.trim.fasta
	# paired.trim.contigs.scrap.fasta
	# paired.trim.contigs.groups

	# ALL SEQUENCES GROUPED...


# Now pull out sequences for each primer...
# The steps completed so far are the same for each gene and primer, 
# so do not delete paired.trim.contigs.groups or paired.trim.contigs.trim.fasta


# 1.3 Pull out COI sequences (separately for each primer)

	# nohup mothur "#get.groups(group=paired.trim.contigs.groups, groups=BR5, fasta=paired.trim.contigs.trim.fasta)" &
	mothur "#get.groups(group=paired.trim.contigs.groups, groups=BR5, fasta=paired.trim.contigs.trim.fasta)"
	# Outputs *.pick* files: 
	
	# Selected 2312760 sequences from your fasta file.
	# Selected 2312760 sequences from your group file.

	# Output File names: 
	# paired.trim.contigs.trim.pick.fasta
	# paired.trim.contigs.pick.groups
	mv paired.trim.contigs.trim.pick.fasta coi_BR5.pick.fasta
	mv paired.trim.contigs.pick.groups coi_BR5.pick.groups


	mothur "#get.groups(group=paired.trim.contigs.groups, groups=mam_ckt, fasta=paired.trim.contigs.trim.fasta)"
	# Outputs *.pick* files: 
	
	# Selected 11110 sequences from your fasta file.
	# Selected 11110 sequences from your group file.

	# Output File names: 
	# paired.trim.contigs.trim.pick.fasta
	# paired.trim.contigs.pick.groups
	mv paired.trim.contigs.trim.pick.fasta coi_mam_ckt.pick.fasta
	mv paired.trim.contigs.pick.groups coi_mam_ckt.pick.groups


	mothur "#get.groups(group=paired.trim.contigs.groups, groups=mam_ckt_230, fasta=paired.trim.contigs.trim.fasta)"
	# Outputs *.pick* files: 
	
	# Selected 2661214 sequences from your fasta file.
	# Selected 2661214 sequences from your group file.

	# Output File names: 
	# paired.trim.contigs.trim.pick.fasta
	# paired.trim.contigs.pick.groups
	mv paired.trim.contigs.trim.pick.fasta coi_mam_ckt_230.pick.fasta
	mv paired.trim.contigs.pick.groups coi_mam_ckt_230.pick.groups



# Get original fasta using named sequences from pick file

	mothur "#list.seqs(fasta=coi_BR5.pick.fasta)"
	# Output File Names: 
	# coi_BR5.pick.accnos

	mothur "#get.seqs(accnos=coi_BR5.pick.accnos, group=paired.contigs.groups)"
	# Selected 2312760 sequences from your group file.
	# Output File Names: 
	# paired.contigs.pick.groups

	# rename group file
	mv paired.contigs.pick.groups coi_BR5.groups

	# rename fasta to *_raw.fasta
	mv coi_BR5.pick.fasta coi_BR5_raw.fasta


	mothur "#list.seqs(fasta=coi_mam_ckt.pick.fasta)"
	# Output File Names: 
	# coi_mam_ckt.pick.accnos

	mothur "#get.seqs(accnos=coi_mam_ckt.pick.accnos, group=paired.contigs.groups)"
	# Selected 11110 sequences from your group file.
	# Output File Names: 
	# paired.contigs.pick.groups

	# rename group file
	mv paired.contigs.pick.groups coi_mam_ckt.groups

	# rename fasta to *_raw.fasta
	mv coi_mam_ckt.pick.fasta coi_mam_ckt_raw.fasta


	mothur "#list.seqs(fasta=coi_mam_ckt_230.pick.fasta)"
	# Output File Names: 
	# coi_mam_ckt_230.pick.accnos

	mothur "#get.seqs(accnos=coi_mam_ckt_230.pick.accnos, group=paired.contigs.groups)"
	# Selected 2661214 sequences from your group file.
	# Output File Names: 
	# paired.contigs.pick.groups

	# rename group file
	mv paired.contigs.pick.groups coi_mam_ckt_230.groups

	# rename fasta to *_raw.fasta
	mv coi_mam_ckt_230.pick.fasta coi_mam_ckt_230_raw.fasta


	# remove unecessary files
	rm paired.contigs*
	rm paired.trim*
	rm paired.scrap*


## Format for dada2

### BR5 ###

	# coi_BR5.pick.accnos

# looping over all files

	# pulling out coi_BR5
	for f in ../raw_data/*.fastq
		do
 		echo "Processing $f"
		mothur "#get.seqs(fastq=$f, accnos=coi_BR5.pick.accnos)"
 	done

	# separating fastq
	for f in ../raw_data/*.pick.fastq
		do
 		echo "Processing $f"
		mothur "#fastq.info(fastq=$f)"
		# output: *.fasta, *.fastq
		rm $f
	done

	# trim primers off
	for f in ../raw_data/*.fasta
		do
 		echo "Processing $f"
		mothur "#trim.seqs(fasta=$f, qfile=${f%.*}.qual,  oligos=/data/kruger/illumina/oligos_forwardFR.txt, keepforward=F, checkorient=T, pdiffs=2, processors=4)"
		# output: *.pick.trim.fasta, *pick.trim.qual
	done

	rm ../raw_data/*.scrap*
	rm ../raw_data/*.pick.fasta
	rm ../raw_data/*.pick.qual
	rm ../raw_data/*.pick.groups

	# re-merge into fastq	
	for f in ../raw_data/*.fasta
		do
 		echo "Processing $f"
		mothur "#make.fastq(fasta=$f, qfile=${f%.*}.qual)"
		# output: *.pick.trim.fastq
	done

	rm ../raw_data/*.fasta
	rm ../raw_data/*.qual

# Remove mothur log files
	rm mothur*

# Move *.pick.trim.fastq files to coi_BR5_fastq folder
	
	sudo mkdir /data/kruger/illumina/coi_BR5_fastq
	sudo mv ../raw_data/*pick.trim.fastq /data/kruger/illumina/coi_BR5_fastq/

# Make coi_BR5_fastq folder and symlink to original files
	sudo mkdir coi_BR5_fastq
	sudo ln -s /data/kruger/illumina/coi_BR5_fastq/* coi_BR5_fastq/	

ls ../raw_data* | wc -l
ls coi_BR5_fastq/* | wc -l
# raw_data and coi_BR5_fastq have 176 files



### mam_ckt ###

	# coi_mam_ckt.pick.accnos

# looping over all files

	# pulling out coi_mam_ckt
	for f in ../raw_data/*.fastq
		do
 		echo "Processing $f"
		mothur "#get.seqs(fastq=$f, accnos=coi_mam_ckt.pick.accnos)"
 	done

	# separating fastq
	for f in ../raw_data/*.pick.fastq
		do
 		echo "Processing $f"
		mothur "#fastq.info(fastq=$f)"
		# output: *.fasta, *.fastq
		rm $f
	done

	# trim primers off
	for f in ../raw_data/*.fasta
		do
 		echo "Processing $f"
		mothur "#trim.seqs(fasta=$f, qfile=${f%.*}.qual,  oligos=/data/kruger/illumina/oligos_forwardFR.txt, keepforward=F, checkorient=T, pdiffs=2, processors=4)"
		# output: *.pick.trim.fasta, *pick.trim.qual
	done

	rm ../raw_data/*.scrap*
	rm ../raw_data/*.pick.fasta
	rm ../raw_data/*.pick.qual
	rm ../raw_data/*.pick.groups

	# re-merge into fastq	
	for f in ../raw_data/*.fasta
		do
 		echo "Processing $f"
		mothur "#make.fastq(fasta=$f, qfile=${f%.*}.qual)"
		# output: *.pick.trim.fastq
	done

	rm ../raw_data/*.fasta
	rm ../raw_data/*.qual

# Remove mothur log files
	rm mothur*

# Move *.pick.trim.fastq files to coi_mam_ckt_fastq folder
	
	sudo mkdir /data/kruger/illumina/coi_mam_ckt_fastq
	sudo mv ../raw_data/*pick.trim.fastq /data/kruger/illumina/coi_mam_ckt_fastq/

# Make coi_mam_ckt_fastq folder and symlink to original files
	sudo mkdir coi_mam_ckt_fastq
	sudo ln -s /data/kruger/illumina/coi_mam_ckt_fastq/* coi_mam_ckt_fastq/	

ls ../raw_data* | wc -l
ls coi_mam_ckt_fastq/* | wc -l
# raw_data and coi_mam_ckt_fastq have 176 files



### mam_ckt_230 ###

	# coi_mam_ckt_230.pick.accnos

# looping over all files

	# pulling out coi_mam_ckt_230
	for f in ../raw_data/*.fastq
		do
 		echo "Processing $f"
		mothur "#get.seqs(fastq=$f, accnos=coi_mam_ckt_230.pick.accnos)"
 	done

	# separating fastq
	for f in ../raw_data/*.pick.fastq
		do
 		echo "Processing $f"
		mothur "#fastq.info(fastq=$f)"
		# output: *.fasta, *.fastq
		rm $f
	done

	# trim primers off
	for f in ../raw_data/*.fasta
		do
 		echo "Processing $f"
		mothur "#trim.seqs(fasta=$f, qfile=${f%.*}.qual,  oligos=/data/kruger/illumina/oligos_forwardFR.txt, keepforward=F, checkorient=T, pdiffs=2, processors=4)"
		# output: *.pick.trim.fasta, *pick.trim.qual
	done

	rm ../raw_data/*.scrap*
	rm ../raw_data/*.pick.fasta
	rm ../raw_data/*.pick.qual
	rm ../raw_data/*.pick.groups

	# re-merge into fastq	
	for f in ../raw_data/*.fasta
		do
 		echo "Processing $f"
		mothur "#make.fastq(fasta=$f, qfile=${f%.*}.qual)"
		# output: *.pick.trim.fastq
	done

	rm ../raw_data/*.fasta
	rm ../raw_data/*.qual

# Remove mothur log files
	rm mothur*

# Move *.pick.trim.fastq files to coi_mam_ckt_230_fastq folder
	
	sudo mkdir /data/kruger/illumina/coi_mam_ckt_230_fastq
	sudo mv ../raw_data/*pick.trim.fastq /data/kruger/illumina/coi_mam_ckt_230_fastq/

# Make coi_mam_ckt_230_fastq folder and symlink to original files
	sudo mkdir coi_mam_ckt_230_fastq
	sudo ln -s /data/kruger/illumina/coi_mam_ckt_230_fastq/* coi_mam_ckt_230_fastq/	

ls ../raw_data* | wc -l
ls coi_mam_ckt_230_fastq/* | wc -l
# raw_data and coi_mam_ckt_230_fastq have 176 files
