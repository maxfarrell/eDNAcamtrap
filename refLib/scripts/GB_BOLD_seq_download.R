# Combined download of Vertebrate Reference Sequences 
# from BOLD (via 'bold') and GB (via 'rentrez')

# Maxwell J. Farrell
# maxwellfarrell@gmail.com

# Libraries 
require(rentrez); packageVersion("rentrez") # 1.1.0
require(bold); packageVersion("bold") # 0.4.0

# Species List
mammList <- read.csv("../data/species_lists/kruger_mamm.csv", as.is=T, header=T)
amphList <- read.csv("../data/species_lists/kruger_amph.csv", as.is=T, header=T)
birdList <- read.csv("../data/species_lists/kruger_birds.csv", as.is=T, header=T)
fishList <- read.csv("../data/species_lists/kruger_fish.csv", as.is=T, header=T)
reptList <- read.csv("../data/species_lists/kruger_rept.csv", as.is=T, header=T)

specList <- as.vector(sort(unique(c(
			mammList$binomial,mammList$synonym,
			amphList$binomial,amphList$synonym,
			birdList$binomial,birdList$synonym,
			fishList$binomial,fishList$synonym,
			reptList$binomial,reptList$synonym
			))))

#Creating variables to collect searches with no hits
notFound_GB <- NULL
notFound_BOLD <- NULL

# Looping over specList
for (i in seq_along(specList)){

	name <- specList[i]

	print(paste0("Species ",i," of ",length(specList)))

	# Rentrez
	searchResult <- entrez_search(db="nuccore", use_history= TRUE, term=paste("(COX1[gene] OR cox1[gene] OR coxI[gene] OR CO1[gene] OR COI[gene] OR Cytochrome c oxidase subunit I[gene] OR cytochrome c oxidase subunit I[gene] OR cytochrome oxidase subunit 1[gene] OR Cytochrome oxidase subunit 1[gene]) AND 0:5000[Sequence Length] AND ", name,"[ORGN]", sep=""))

	if(searchResult$count > 0){
	result_GB <- entrez_fetch(db="nuccore", rettype="fasta", web_history=searchResult$web_history)

	} else { notFound_GB <- c(notFound_GB, name)}



	# BOLD
	# Re-setting result_BOLD so that if bold_seqspec returns nothing, the previous result is not evalutated
	result_BOLD <- NULL

	#Accessing BOLD
	try(result_BOLD <- bold_seqspec(taxon=name, sepfasta=F, format="tsv"))

	#Proceeding only if search is returned and contains COI sequences	
	if (!is.null(dim(result_BOLD))){

		if (isTRUE(any(result_BOLD$markercode=="COI-5P"))) { 

		result_BOLD <- result_BOLD[result_BOLD$markercode=="COI-5P",]

		#Reformatting output of result_BOLD to FASTA
		seqs <- as.vector(NULL)

			for (k in 1:nrow(result_BOLD)){
			seqs[k] <- as.vector(paste0(">",result_BOLD$processid[k]," ",result_BOLD$species_name[k]," ",result_BOLD$markercode[k]," ","\n", 
		 	paste(sapply(seq(from=1, to=nchar(result_BOLD$nucleotides[k]), by=70), 
		 		function(j) substr(result_BOLD$nucleotides[k], j, j+70)), collapse="\n")))
			}
		}

	} else { notFound_BOLD <- c(notFound_BOLD, name) }


	# Combining GB and BOLD output writing to file ../data/FASTAs/*.fasta

	if (!name %in% notFound_GB | !name %in% notFound_BOLD) {
		
		filename <- paste("../data/FASTAs/", gsub(" ","_", name), ".fasta", sep="")	
		file.create(filename)
		fileConn <- file(filename)
		if (!name %in% notFound_BOLD) writeLines(seqs, sep="\n \n", fileConn)
		close(fileConn)
		
		if (!name %in% notFound_GB) write (result_GB, file=filename, append=TRUE) 

	}
	

	# Wait 10 seconds every 10 searches, otherwise wait 2 second
		# Modulus operation
   		if(i %% 10==0) {
   			Sys.sleep(10)
   		}else(Sys.sleep(2))

}

write.csv(notFound_GB, "../output/withSynonyms_notFound_GB.csv", row.names=F)
write.csv(notFound_BOLD, "../output/withSynonyms_notFound_BOLD.csv", row.names=F)
write.csv(intersect(notFound_BOLD,notFound_GB), "../output/withSynonyms_notFound_either.csv", row.names=F)

# intersect(notFound_BOLD,notFound_GB)