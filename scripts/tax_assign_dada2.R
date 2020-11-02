## Kruger COI Taxonomy assignment

require("dada2"); packageVersion("dada2") # 1.10.0

# COI Sequence Table
if (!file.exists("../results/coi_merged_Fonly_seqtab.rds")){
	
	st_BR5 <- readRDS("../results/coi_BR5_Fonly_seqtab.rds")
	st_mam_ckt <- readRDS("../results/coi_mam_ckt_Fonly_seqtab.rds")
	st_mam_ckt_230 <- readRDS("../results/coi_mam_ckt_230_Fonly_seqtab.rds")
	seqtab <- mergeSequenceTables(st_BR5, st_mam_ckt, st_mam_ckt_230, repeats="sum")
	saveRDS(seqtab, "../results/coi_merged_Fonly_seqtab.rds")	
	
	} else {seqtab <- readRDS("../results/coi_merged_Fonly_seqtab.rds")
}

# Taxnomy Assignment -- Kruger Vertebrate RefLib
if (!file.exists("../results/coi_vertebrates_dada2_phy2species_98.rds")){

	taxKNP_vert_98 <- assignTaxonomy(seqtab, "./refLib/output/Kruger_Vertebrates_refLib_dada2_phy2species.fasta", tryRC=TRUE, multithread=TRUE, minBoot=98)	
	colnames(taxKNP_vert_98) <- c("Phylum","Class","Order","Family","Genus","Species")
	saveRDS(taxKNP_vert_98, "../results/coi_vertebrates_dada2_phy2species_98.rds")

	} else {taxKNP_vert_98 <- readRDS("../results/coi_vertebrates_dada2_phy2species_98.rds")
}

if (!file.exists("../results/coi_vertebrates_dada2_phy2species_95.rds")){

	taxKNP_vert_95 <- assignTaxonomy(seqtab, "./refLib/output/Kruger_Vertebrates_refLib_dada2_phy2species.fasta", tryRC=TRUE, multithread=TRUE, minBoot=95)	
	colnames(taxKNP_vert_95) <- c("Phylum","Class","Order","Family","Genus","Species")
	saveRDS(taxKNP_vert_95, "../results/coi_vertebrates_dada2_phy2species_95.rds")

	} else {taxKNP_vert_95 <- readRDS("../results/coi_vertebrates_dada2_phy2species_95.rds")
}

if (!file.exists("../results/coi_vertebrates_dada2_phy2species_90.rds")){

	taxKNP_vert_90 <- assignTaxonomy(seqtab, "./refLib/output/Kruger_Vertebrates_refLib_dada2_phy2species.fasta", tryRC=TRUE, multithread=TRUE, minBoot=90)	
	colnames(taxKNP_vert_90) <- c("Phylum","Class","Order","Family","Genus","Species")
	saveRDS(taxKNP_vert_90, "../results/coi_vertebrates_dada2_phy2species_90.rds")

	} else {taxKNP_vert_90 <- readRDS("../results/coi_vertebrates_dada2_phy2species_90.rds")
}

if (!file.exists("../results/coi_vertebrates_dada2_phy2species_80.rds")){

	taxKNP_vert_80 <- assignTaxonomy(seqtab, "./refLib/output/Kruger_Vertebrates_refLib_dada2_phy2species.fasta", tryRC=TRUE, multithread=TRUE, minBoot=80)	
	colnames(taxKNP_vert_80) <- c("Phylum","Class","Order","Family","Genus","Species")
	saveRDS(taxKNP_vert_80, "../results/coi_vertebrates_dada2_phy2species_80.rds")

	} else {taxKNP_vert_80 <- readRDS("../results/coi_vertebrates_dada2_phy2species_80.rds")
}

if (!file.exists("../results/coi_vertebrates_dada2_phy2species_70.rds")){

	taxKNP_vert_70 <- assignTaxonomy(seqtab, "./refLib/output/Kruger_Vertebrates_refLib_dada2_phy2species.fasta", tryRC=TRUE, multithread=TRUE, minBoot=70)	
	colnames(taxKNP_vert_70) <- c("Phylum","Class","Order","Family","Genus","Species")
	saveRDS(taxKNP_vert_70, "../results/coi_vertebrates_dada2_phy2species_70.rds")

	} else {taxKNP_vert_70 <- readRDS("../results/coi_vertebrates_dada2_phy2species_70.rds")
}

if (!file.exists("../results/coi_vertebrates_dada2_phy2species_60.rds")){

	taxKNP_vert_60 <- assignTaxonomy(seqtab, "./refLib/output/Kruger_Vertebrates_refLib_dada2_phy2species.fasta", tryRC=TRUE, multithread=TRUE, minBoot=60)	
	colnames(taxKNP_vert_60) <- c("Phylum","Class","Order","Family","Genus","Species")
	saveRDS(taxKNP_vert_60, "../results/coi_vertebrates_dada2_phy2species_60.rds")

	} else {taxKNP_vert_60 <- readRDS("../results/coi_vertebrates_dada2_phy2species_60.rds")
}

if (!file.exists("../results/coi_vertebrates_dada2_phy2species_50.rds")){

	taxKNP_vert_50 <- assignTaxonomy(seqtab, "./refLib/output/Kruger_Vertebrates_refLib_dada2_phy2species.fasta", tryRC=TRUE, multithread=TRUE, minBoot=50)	
	colnames(taxKNP_vert_50) <- c("Phylum","Class","Order","Family","Genus","Species")
	saveRDS(taxKNP_vert_50, "../results/coi_vertebrates_dada2_phy2species_50.rds")

	} else {taxKNP_vert_50 <- readRDS("../results/coi_vertebrates_dada2_phy2species_50.rds")
}

sort(unique(unname(taxKNP_vert_98)[,6])) # 20 species detected
sort(unique(unname(taxKNP_vert_95)[,6])) # 20 species detected
sort(unique(unname(taxKNP_vert_90)[,6])) # 20 species detected
sort(unique(unname(taxKNP_vert_80)[,6])) # 21 species detected
sort(unique(unname(taxKNP_vert_70)[,6])) # 23 species detected
sort(unique(unname(taxKNP_vert_60)[,6])) # 28 species detected
sort(unique(unname(taxKNP_vert_50)[,6])) # 42 species detected

# Taxonomy Assignment with terrimporter COI RDP training database
# https://github.com/terrimporter/CO1Classifier/releases/tag/v3.0
# Formatted to be used with dada2
# taxTP <- assignTaxonomy(seqtab, "./refLib/mydata_training/dada2trainseq.fasta", tryRC=TRUE, multithread=TRUE)
# Memmory allocation limit with full refLib...

if (!file.exists("../results/coi_TP_chordata_tax_80.rds")){

	# Default minBoot=50; using minBoot=80 as this is as good as 98 for curated KNP reference
	taxTP <- assignTaxonomy(seqtab, "./refLib/mydata_training/dada2trainseq_chordata.fasta", tryRC=TRUE, multithread=TRUE, minBoot=80)
	colnames(taxTP) <- c("Phylum","Class","Order","Family","Genus","Species")
	saveRDS(taxTP, "../results/coi_TP_chordata_tax_80.rds")

	} else {taxTP <- readRDS("../results/coi_TP_chordata_tax_80.rds")
}

sort(unique(unname(taxTP)[,6])) # 20 species detected


# Taxonomy Assignment with MIDORI UNIQUE COI RDP training database

if (!file.exists("../results/coi_MIDORI_tax_80.rds")){

	# Default minBoot=50; using minBoot=80 as this is as good as 98 for curated KNP reference
	taxMID <- assignTaxonomy(seqtab, "./MIDORI_UNIQUE_COI_dada2_Phylum_Species_chordata.fasta", tryRC=TRUE, multithread=TRUE, minBoot=80)
	colnames(taxMID) <- c("Phylum","Class","Order","Family","Genus","Species")
	saveRDS(taxMID, "../results/coi_MIDORI_tax_80.rds")

	} else {taxMID <- readRDS("../results/coi_MIDORI_tax_80.rds")
}

sort(unique(unname(taxMID)[,6]))# 18 species detected

