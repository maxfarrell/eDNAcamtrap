# Download mitochondrial genomes and excise CO1 gene

require(rentrez); packageVersion("rentrez") # 1.1.0
source("PrimerMiner/Download_mito.R")
source("PrimerMiner/Mito_GB2fasta_MF.R")

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


# Wait 0.5 seconds between each query to avoid timeout (API has a maximum of three requests per second)
for (i in seq_along(specList)){

	Download_mito(specList[i], folder="../data/mito_genomes")
	Sys.sleep(0.5)
}

Mito_GB2fasta("../data/mito_genomes")
