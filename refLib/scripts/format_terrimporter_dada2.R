# Format the terrimporter COI RDP refLib for dada2
# https://github.com/terrimporter/CO1Classifier/releases/tag/v3.0

# Download + unzip training data
# wget "https://github.com/terrimporter/CO1Classifier/releases/download/v3.0-ref/CO1v3_training.tar.gz"
# tar -xzf CO1v3_training.tar.gz

# cd to "mydata_training" folder

# Starting with a separate file of headers speeds things up 
# (compared to reading from a system call directly)
system("grep \">\" mytrainseq.fasta > headers.txt")
tax <- read.table("headers.txt", sep=";")
seqid <- sub(" .*$","", tax[,1])
tax[,1] <- sub("^.* ","", tax[,1])
tax <- cbind(seqid,tax)
names(tax) <- c("seqid","cellular","domain","kingdom","phylum","class","order","family","genus","species") 

seqnames <- with(tax, paste0(">",phylum,";",class,";",order,";",family,";",genus,";",species))
write.table(seqnames,file="dada2_headers.txt",sep="\n",row.names=FALSE, quote = FALSE, col.names=FALSE)

# generating dada2trainseq.fasta via bash
# modified from https://www.biostars.org/p/103089/ to also convert to uppercase nucleotides
system("awk 'NR%2==0' mytrainseq.fasta | tr [a-z] [A-Z] | paste -d'\\n' dada2_headers.txt - > dada2trainseq.fasta")

# Bash command
# awk 'NR%2==0' mytrainseq.fasta | tr [a-z] [A-Z] | paste -d'\n' dada2_headers.txt - > dada2trainseq.fasta
# sort(unique(tax[,2])) # cellularOrganisms
# sort(unique(tax[,3])) # Domain
# sort(unique(tax[,4])) # Kingdom (least informative)
# sort(unique(tax[,5])) # Phylum
# sort(unique(tax[,6])) # Class
# sort(unique(tax[,7])) # Order
# sort(unique(tax[,8])) # Family
# sort(unique(tax[,9])) # Genus
# sort(unique(tax[,10])) # Species

# Bash to subset to chordata
# Select lines that include Chordata and the following line (sequence)
# grep -A 1 "Chordata" dada2trainseq.fasta > dada2trainseq_chordata.fasta
