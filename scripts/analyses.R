## Kruger COI community analyses ##

require(dada2); packageVersion("dada2") # 1.10.0
require(phyloseq); packageVersion("phyloseq") # 1.26.0
require(ggplot2); packageVersion("ggplot2") # 3.1.0
require(ape); packageVersion("ape")# 5.2
require(dplyr); packageVersion("dplyr")# 0.8.0.1
require(vegan); packageVersion("vegan")# 2.5.3
require(lubridate); packageVersion("lubridate")# 1.7.4
require(rstan); packageVersion("rstan") # 2.18.2
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
require(bayesplot);packageVersion("bayesplot")# 1.5.0
require(kableExtra)
require(xtable)
require(tidyr);packageVersion("tidyr")# 1.0.2

# seqtab
seqtab <- readRDS("../results/coi_merged_Fonly_seqtab.rds")

# tax assignments
taxKNP_vert_98 <- readRDS("../results/coi_vertebrates_dada2_phy2species_98.rds")
taxKNP_vert_95 <- readRDS("../results/coi_vertebrates_dada2_phy2species_95.rds")
taxKNP_vert_90 <- readRDS("../results/coi_vertebrates_dada2_phy2species_90.rds")
taxKNP_vert_80 <- readRDS("../results/coi_vertebrates_dada2_phy2species_80.rds")
taxKNP_vert_70 <- readRDS("../results/coi_vertebrates_dada2_phy2species_70.rds")
taxKNP_vert_60 <- readRDS("../results/coi_vertebrates_dada2_phy2species_60.rds")
taxKNP_vert_50 <- readRDS("../results/coi_vertebrates_dada2_phy2species_50.rds")
taxTP <- readRDS("../results/coi_TP_chordata_tax_80.rds")
taxMID <- readRDS("../results/coi_MIDORI_tax_80.rds")


# Subsetting to samples in camera trap analysis
cam_samples <- c("Max-BLANK-2","Max-BLANK-3","Max-DLP-8A","Max-DLP-8A-S",
         "Max-DLP-8A-XS","Max-DLP-8B","Max-DLP-8B-S","Max-DLP-8B-XS",
         "Max-HOY-2A","Max-HOY-2B",
         "Max-HOY-3A","Max-HOY-3A-S",
         "Max-HOY-3A-XS","Max-HOY-3B","Max-HOY-3B-S","Max-HOY-3B-XS",
         "Max-KWA-6A","Max-KWA-6A-S","Max-KWA-6A-XS","Max-KWA-6B",
         "Max-KWA-6B-S","Max-KWA-6B-XS","Max-NGO-2A","Max-NGO-2B",
         "Max-NGO-3A","Max-NGO-3A-S","Max-NGO-3A-XS","Max-NGO-3B",
         "Max-NGO-3B-S","Max-NGO-3B-XS","Max-NWA-2A","Max-NWA-2A-S",
         "Max-NWA-2A-XS","Max-NWA-2B","Max-NWA-2B-S","Max-NWA-2B-XS",
         "Max-NYA-2A","Max-NYA-2B","Max-NYA-3A","Max-NYA-3B",
         "Max-NYA-4A","Max-NYA-4A-S","Max-NYA-4A-XS","Max-NYA-4B",
         "Max-NYA-4B-S","Max-NYA-4B-XS")

seqtab_cam <- seqtab[rownames(seqtab)%in%cam_samples,]
seqtab_cam <- seqtab_cam[,colSums(seqtab_cam)>0]
dim(seqtab_cam) # 46 x 5501
# colnames(seqtab_cam)

taxKNP_vert_98 <- taxKNP_vert_98[rownames(taxKNP_vert_98)%in%colnames(seqtab_cam),]
taxKNP_vert_95 <- taxKNP_vert_95[rownames(taxKNP_vert_95)%in%colnames(seqtab_cam),]
taxKNP_vert_90 <- taxKNP_vert_90[rownames(taxKNP_vert_90)%in%colnames(seqtab_cam),]
taxKNP_vert_80 <- taxKNP_vert_80[rownames(taxKNP_vert_80)%in%colnames(seqtab_cam),]
taxKNP_vert_70 <- taxKNP_vert_70[rownames(taxKNP_vert_70)%in%colnames(seqtab_cam),]
taxKNP_vert_60 <- taxKNP_vert_60[rownames(taxKNP_vert_60)%in%colnames(seqtab_cam),]
taxKNP_vert_50 <- taxKNP_vert_50[rownames(taxKNP_vert_50)%in%colnames(seqtab_cam),]
taxTP <- taxTP[rownames(taxTP)%in%colnames(seqtab_cam),]
taxMID <- taxMID[rownames(taxMID)%in%colnames(seqtab_cam),]

# Taxonomic assignment richness plots for KNP
plot_dat <- data.frame(var=c("KNP Vertebrates - minBoot 50",
               "KNP Vertebrates - minBoot 60",
               "KNP Vertebrates - minBoot 70",
               "KNP Vertebrates - minBoot 80",
               "KNP Vertebrates - minBoot 90",
               "KNP Vertebrates - minBoot 95",
               "KNP Vertebrates - minBoot 98",
               "Porter Chordata - minBoot 80",
               "MIDORI COI - minBoot 80"),
            spec_r=c(sum(!is.na(unique(unname(taxKNP_vert_50)[,6]))),
                sum(!is.na(unique(unname(taxKNP_vert_60)[,6]))),
                sum(!is.na(unique(unname(taxKNP_vert_70)[,6]))),
                sum(!is.na(unique(unname(taxKNP_vert_80)[,6]))),
                sum(!is.na(unique(unname(taxKNP_vert_90)[,6]))),
                sum(!is.na(unique(unname(taxKNP_vert_95)[,6]))),
                sum(!is.na(unique(unname(taxKNP_vert_98)[,6]))),
                sum(!is.na(unique(unname(taxTP)[,6]))),
                sum(!is.na(unique(unname(taxMID)[,6])))),
            gen_r=c(sum(!is.na(unique(unname(taxKNP_vert_50)[,5]))),
                sum(!is.na(unique(unname(taxKNP_vert_60)[,5]))),
                sum(!is.na(unique(unname(taxKNP_vert_70)[,5]))),
                sum(!is.na(unique(unname(taxKNP_vert_80)[,5]))),
                sum(!is.na(unique(unname(taxKNP_vert_90)[,5]))),
                sum(!is.na(unique(unname(taxKNP_vert_95)[,5]))),
                sum(!is.na(unique(unname(taxKNP_vert_98)[,5]))),
                sum(!is.na(unique(unname(taxTP)[,5]))),
                sum(!is.na(unique(unname(taxMID)[,5])))),
            ord_r=c(sum(!is.na(unique(unname(taxKNP_vert_50)[,3]))),
                sum(!is.na(unique(unname(taxKNP_vert_60)[,3]))),
                sum(!is.na(unique(unname(taxKNP_vert_70)[,3]))),
                sum(!is.na(unique(unname(taxKNP_vert_80)[,3]))),
                sum(!is.na(unique(unname(taxKNP_vert_90)[,3]))),
                sum(!is.na(unique(unname(taxKNP_vert_95)[,3]))),
                sum(!is.na(unique(unname(taxKNP_vert_98)[,3]))),
                sum(!is.na(unique(unname(taxTP)[,3]))),
                sum(!is.na(unique(unname(taxMID)[,3])))))
View(plot_dat)

blank_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), 
                     axis.title=element_text(size=14))

sp_rich_plot <-  ggplot(plot_dat, aes(x = var, y=spec_r)) + 
              geom_bar(stat="identity",colour=FALSE) +
              blank_theme + coord_flip() + 
              scale_x_discrete(limits = rev(levels(plot_dat$var)))+
              theme(axis.title.y=element_blank()) + 
            labs(y="Number of species identified")
sp_rich_plot
ggsave("../plots_tables/sp_rich_plot.pdf",sp_rich_plot, width=7, height=4)


gen_rich_plot <-  ggplot(plot_dat, aes(x = var, y=gen_r)) + 
              geom_bar(stat="identity",colour=FALSE) +
              blank_theme + coord_flip() + 
              scale_x_discrete(limits = rev(levels(plot_dat$var)))+
              theme(axis.title.y=element_blank()) + 
            labs(y="Number of genera identified")
gen_rich_plot
ggsave("../plots_tables/gen_rich_plot.pdf",gen_rich_plot, width=7, height=4)


ord_rich_plot <-  ggplot(plot_dat, aes(x = var, y=ord_r)) + 
              geom_bar(stat="identity",colour=FALSE) +
              blank_theme + coord_flip() + 
              scale_x_discrete(limits = rev(levels(plot_dat$var)))+
              theme(axis.title.y=element_blank()) + 
            labs(y="Number of orders identified")
ord_rich_plot
ggsave("../plots_tables/ord_rich_plot.pdf",ord_rich_plot, width=7, height=4)


# Table
names(plot_dat) <- c("Approach","# Species","# Genera","# Orders")
# require(xtable)
# print(xtable(plot_dat, digits=0),file="../plots_tables/tax_assignments.tex")


# Taxonomic groups
vert_80_tax <- unique(data.frame(class = taxKNP_vert_80[,2], order = taxKNP_vert_80[,3],
            species = taxKNP_vert_80[,6]))
dim(vert_80_tax)
# assigned to species
vert_80_tax <- vert_80_tax[!is.na(vert_80_tax$species),]
table(vert_80_tax$class)


# Stacked barplots across classes and Jaccard plots for KNP 50, KNP 80, MIDORI, PORTER

KNP50 <-  data.frame(tax="KNP50",unique(unname(taxKNP_vert_50)[,c(2,3,4,6)]))
names(KNP50) <- c("tax","Class","Order","Family","Species")
# KNP50 <- KNP50[!is.na(KNP50$species),]

KNP80 <-  data.frame(tax="KNP80",unique(unname(taxKNP_vert_80)[,c(2,3,4,6)]))
names(KNP80) <- c("tax","Class","Order","Family","Species")

TP <-  data.frame(tax="TP",unique(unname(taxTP)[,c(2,3,4,6)]))
names(TP) <- c("tax","Class","Order","Family","Species")

MID <-  data.frame(tax="MID",unique(unname(taxMID)[,c(2,3,4,6)]))
names(MID) <- c("tax","Class","Order","Family","Species")
MID$Class[MID$Class=="Actinopteri"] <- "Actinopterygii"

tax_dat <- rbind(KNP50, KNP80, TP, MID)

# Remove NA species
tax_dat <- tax_dat[!is.na(tax_dat$Species),]

# Make latex table of this for SM (species only)
tax_dat %>% filter(Class=="Actinopterygii")

tax_ident_tabl <- tax_dat %>% select(!c(Order, Family)) %>% arrange(tax, Class, Species)
tax_ident_tabl$Species <- gsub("_"," ",tax_ident_tabl$Species)
tax_ident_tabl$Species <- paste0("textit{",tax_ident_tabl$Species,"}")
require(xtable)
print(xtable(tax_ident_tabl, digits=0),file="../plots_tables/tax_assignments_class_species.tex",include.rownames=FALSE,
                      tabular.environment = "longtable")

class_counts <- tax_dat %>% 
  count(Class, tax) %>% 
  group_by(Class, tax)

# expanding to include zero counts
class_counts <- expand.grid(Class=unique(class_counts$Class), tax=unique(class_counts$tax)) %>% data.frame %>% left_join(class_counts)

class_counts$Class <- as.factor(class_counts$Class)

class_count_plot <- 
  ggplot(class_counts, aes(fill=Class, y=n, x=tax)) + 
    geom_bar(position="dodge", stat="identity") + 
    scale_fill_manual(values = c(
                 "#56B4E9",#Fish
                   "#0072B2",#Amphibians
                 "#E69F00",#Birds
                 "#D55E00",#Mammals
                 "#009E73",#green
                 "#672262"
                 ),drop=FALSE) +
    xlab("Taxonomy") + ylab("Number of unique species") +
    scale_y_continuous(expand = c(0, 0)) + theme_classic()

ggsave("../plots_tables/class_tax_plot.pdf",class_count_plot, width=6, height=4)


# BetaDiv plot data for species identified
KNP50_species <-  unique(unname(taxKNP_vert_50)[,6])
KNP50_species <- KNP50_species[!is.na(KNP50_species)]
KNP50_species <- data.frame(tax="KNP50", species=KNP50_species)

KNP80_species <-  unique(unname(taxKNP_vert_80)[,6])
KNP80_species <- KNP80_species[!is.na(KNP80_species)]
KNP80_species <- data.frame(tax="KNP80", species=KNP80_species)

TP_species <-  unique(unname(taxTP)[,6])
TP_species <- TP_species[!is.na(TP_species)]
TP_species <- data.frame(tax="TP", species=TP_species)

MID_species <-  unique(unname(taxMID)[,6])
MID_species <- MID_species[!is.na(MID_species)]
MID_species <- data.frame(tax="MID", species=MID_species)

beta_dat <- rbind(KNP50_species, KNP80_species, TP_species, MID_species)

require(betapart); packageVersion("betapart")#1.5.0
com_pa <- table(beta_dat)
com.dist <- vegdist(com_pa)
com_pa[com_pa>0] <-1
pair.s <- beta.pair(com_pa)


pdf("../plots_tables/betapart.pdf", width=8, height=4)
par(mfrow=c(1,3))

plot(hclust(pair.s$beta.sor, method="average"), main="", sub="", xlab="")
title(xlab=expression(beta[sor]), line=0.3)

plot(hclust(pair.s$beta.sim, method="average"), main="", sub="", xlab="")
title(xlab=expression(beta[sim]), line=0.3)

plot(hclust(pair.s$beta.sne, method="average"), main="", sub="", xlab="")
title(xlab=expression(beta[sne]), line=0.3)

dev.off()



#############################
# Activity overlap analysis #
#############################

require(ggplot2); packageVersion("ggplot2") # 3.3.0
require(dplyr); packageVersion("dplyr")# 0.8.5
require(lubridate); packageVersion("lubridate")# 1.7.8
require(overlap);packageVersion("overlap")# 0.3.3
require(Hmisc);packageVersion("Hmisc")# 4.4.0
require(dada2); packageVersion("dada2") # 1.16.0
require(phyloseq); packageVersion("phyloseq") # 1.32.0
require(ggcorrplot);packageVersion("ggcorrplot")#0.1.3


# loading data
load("../../data/merged_eDNA_camtrap_data_nov20_2019.RData")

# merging & cleaning species
trap$species[trap$species=="whiterhino"] <- "rhino"
trap <- trap[trap$species!="frog",]
trap <- trap[trap$species!="fox",]


# CALCULATE CAMERA OVERLAPS AS COUNTS PER SITE

# Subset to sample x species
trap_spec <- trap %>% select(sample_num, species)

a <- t(table(trap_spec))
a[a>1] <- 1
a <- a%*%t(a)
cam_mat <- a
diag(cam_mat) <- NA
cam_mat[upper.tri(cam_mat)] <- NA
cam_mat <- apply(cam_mat, 1, rev)

grad_len <- max(cam_mat, na.rm=T)

p <- ggcorrplot(t(cam_mat), outline.col="grey", show.diag=FALSE) + 
        scale_fill_gradientn(limits=c(0,grad_len), breaks=c(0:grad_len), colors = c("#fffff2", "#8dadc4", "#E46726"), guide="legend")
p <- p + labs(fill = "overlap via \ncamera samples \n")
p 

ggsave("../plots_tables/camera_overlap.pdf", p, height=8, width=8)


###################################
## ggcorrplot for eDNA detection ##
###################################

# setting which taxonomy we want to use for community analyses
ps <- ps_80

# subset ps to samples with camera data
ps_trap <- subset_samples(ps, sample_num%in%camera_samples_fullweek)
dim(otu_table(ps_trap))# 42 x 4616

# make seqs object (data frame of taxonomy, sample names, and time collected.)
names(sample_dates)
tax <- tax_table(ps_trap)
edna_taxa <- rownames(tax_table(ps_trap)[!is.na(tax_table(ps_trap)[,6])])
ps_trap <- prune_taxa(edna_taxa, ps_trap)

# Converting ps_trap to dataframe
edna <- psmelt(ps_trap)

seqs <- left_join(edna, sample_dates)

# Removing abundances of 0
seqs <- seqs[seqs$Abundance>0,]

# Subset to sample x species
sample_num_seq <- seqs %>% select(sample_num, Species)
sample_num_seq$Species <- as.character(sample_num_seq$Species)

sp_names <- trap %>% select(species, binomial) %>% unique()
name_translation <- setNames(sp_names$species, sp_names$binomial)
sample_num_seq$species <- name_translation[sample_num_seq$Species]
sample_num_seq <- sample_num_seq[,-2]

sample_num_seq <- sample_num_seq[!is.na(sample_num_seq$species),]

a <- t(table(sample_num_seq))
a[a>1] <- 1
a <- a%*%t(a)
edna_mat <- a
diag(edna_mat) <- NA

intersect(rownames(cam_mat), rownames(edna_mat))
edna_mat <- edna_mat[intersect(rownames(cam_mat), rownames(edna_mat)),intersect(colnames(cam_mat), colnames(edna_mat))]
edna_mat <- edna_mat[nrow(edna_mat):1,]
edna_mat[lower.tri(edna_mat)] <- NA
edna_mat <- edna_mat[nrow(edna_mat):1,]

grad_len <- max(edna_mat, na.rm=T)

p2 <- ggcorrplot(t(edna_mat), outline.col="grey", show.diag=FALSE) + 
        scale_fill_gradientn(limits=c(0,grad_len), breaks=c(0:grad_len), colors = c("#fffff2", "#8dadc4", "#E46726"), guide="legend")
p2 <- p2 + labs(fill = "overlap via \neDNA samples \n")
p2 

ggsave("../plots_tables/edna_overlap.pdf", p2, height=8, width=8)


cam_mat_small <- cam_mat[intersect(rownames(edna_mat), rownames(cam_mat)),intersect(colnames(edna_mat), colnames(cam_mat))]

rownames(cam_mat_small)
rownames(edna_mat)
colnames(cam_mat_small)
colnames(edna_mat)

cor(c(cam_mat_small), c(edna_mat), use="complete.obs")# 0.39

grad_len <- max(cam_mat_small, na.rm=T)

p3 <- ggcorrplot(t(cam_mat_small), outline.col="grey", show.diag=FALSE) + 
        scale_fill_gradientn(limits=c(0,grad_len), breaks=c(0:grad_len), colors = c("#fffff2", "#8dadc4", "#E46726"), guide="legend")
p3 <- p3 + labs(fill = "overlap via \ncamera samples \n")
p3 

ggsave("../plots_tables/camera_overlap_small.pdf", p3, height=8, width=8)


#####################################
## Re-do for 36 hour camera window ##
#####################################

trap_36  <- unique(left_join(trap,seqs) %>% select(species, n_waterhole, timestamp, collected, sample_num))
trap_36$tminus <- trap_36$collected - trap_36$timestamp
trap_36$tminus <- as.numeric(trap_36$tminus/(24*60)) #now in days
trap_36 <- trap_36[trap_36$tminus<1.5,]

# Subset to sample x species
trap_spec <- trap_36 %>% select(sample_num, species)

a <- t(table(trap_spec))
a[a>1] <- 1
a <- a%*%t(a)
cam_mat <- a
diag(cam_mat) <- NA
cam_mat[upper.tri(cam_mat)] <- NA
cam_mat <- apply(cam_mat, 1, rev)

grad_len <- max(cam_mat, na.rm=T)

p4 <- ggcorrplot(t(cam_mat), outline.col="grey", show.diag=FALSE) + 
        scale_fill_gradientn(limits=c(0,grad_len), breaks=c(0:grad_len), colors = c("#fffff2", "#8dadc4", "#E46726"), guide="legend")
p4 <- p4 + labs(fill = "overlap via \ncamera samples \n")
p4 

ggsave("../plots_tables/camera_overlap_36h.pdf", p4, height=8, width=8)

cam_mat_small <- cam_mat[intersect(rownames(edna_mat), rownames(cam_mat)),intersect(colnames(edna_mat), colnames(cam_mat))]
cor(c(cam_mat_small), c(edna_mat), use="complete.obs")# 0.57

# Plotting pruned 36h cam_mat

grad_len <- max(cam_mat_small, na.rm=T)

p5 <- ggcorrplot(t(cam_mat_small), outline.col="grey", show.diag=FALSE) + 
        scale_fill_gradientn(limits=c(0,grad_len), breaks=c(0:grad_len), colors = c("#fffff2", "#8dadc4", "#E46726"), guide="legend")
p5 <- p5 + labs(fill = "overlap via \ncamera samples \n")
p5 

ggsave("../plots_tables/camera_overlap_36h_small.pdf", p5, height=8, width=8)


####################################
# Comparing eDNA & Camera trapping #
####################################

# loading data
load("../data/merged_eDNA_camtrap_data_nov20_2019.RData")

# setting which taxonomy we want to use for community analyses
ps <- ps_80

###################################################################
## Calculating numbers of ASVs and sequences across main samples ##
###################################################################
# subset ps to samples with camera data
ps_trap <- subset_samples(ps, sample_num%in%camera_samples_fullweek)
dim(otu_table(ps_trap))# 42 x 4616

sort(unique(sample_dates$sample_num))

ps_main <- prune_samples(sample_data(ps_trap)$sample_num%in%sample_dates$sample_num, ps_trap)
ps_main <- prune_samples(is.na(sample_data(ps_main)$S_XS), ps_main)

# remove ASVs with zero abundance
ps_main <- prune_taxa(taxa_sums(ps_main) >= 1, ps_main)

# Total ASVs and sequence reads
ncol(otu_table(ps_main))# 2,986
sum(otu_table(ps_main))# 921,086
ps_main_mamms <- subset_taxa(ps_main, Class == "Mammalia")
ncol(otu_table(ps_main_mamms))# 32
sum(otu_table(ps_main_mamms))# 4,255

ps_main_aves <- subset_taxa(ps_main, Class == "Aves")
ncol(otu_table(ps_main_aves))# 39
sum(otu_table(ps_main_aves))# 1,701

# 921086 - 4255 - 1701 # 915,130 not assigned to a class

########################################################################

# make seqs object (data frame of taxonomy, sample names, and time collected.)
names(sample_dates)
tax <- tax_table(ps_trap)
edna_taxa <- rownames(tax_table(ps_trap)[!is.na(tax_table(ps_trap)[,6])])
ps_trap <- prune_taxa(edna_taxa, ps_trap)

# Converting ps_trap to dataframe
edna <- psmelt(ps_trap)

seqs <- left_join(edna, sample_dates)

# Removing abundances of 0
seqs <- seqs[seqs$Abundance>0,]

# including time of sample collection to calclate time since last visit for each species per sample
collected <- seqs %>% select(sample_num, collected) %>% unique()

dat1 <- trap %>% full_join(collected)%>%
        filter(n_waterhole!=0)%>%
        group_by(sample_num, binomial)%>%
        mutate(n_waterhole=sum(n_waterhole), 
            seen=((n_waterhole>0)*1),
            last_visit=as.numeric(unique(collected)-max(timestamp), unit="hours")) %>%
        select(site,sample_num,binomial,n_waterhole,last_visit,seen)%>%
        unique()%>%
        arrange(site, sample_num, binomial)     

names(seqs) <- tolower(names(seqs))
seqs <- subset(seqs, select=-c(sample))
# dim(seqs)

dat2 <- seqs %>%
        filter(is.na(s_xs))%>%
        # group_by(sample_num, AB, species)%>%
        group_by(sample_num, species)%>%
        mutate(detected=1, abundance=sum(abundance))%>%
        select(site, sample_num, species, detected, abundance) %>%
        unique() %>%
        arrange(site, sample_num, species)      
# View(dat2)
names(dat2)[names(dat2)=="species"]<- "binomial"
# head(dat2)

dat3 <- env %>%
        group_by(sample_num)%>%
        summarize(
            temp = mean(temp_c),
            ms = mean(ms_cm),
            do_p = mean(do_percent),
            do_ch = mean(do_ch),
            ft = mean(ft),
            ph = mean(ph))

dat <- full_join(dat1,dat2) %>%
        # arrange(site, sample_num, binomial, AB)
        arrange(site, sample_num, binomial)

dat <- left_join(dat, dat3)

dat$detected[is.na(dat$detected)] <- 0
dat <- dat[dat$binomial!="unknown",]

dim(dat[is.na(dat$binomial),])
# no NA binomial names

# renaming abundance to seq_count
names(dat)[names(dat)=="abundance"] <- "seq_count"

# replacing NA in seq_count with zero
dat$seq_count[is.na(dat$seq_count)] <- 0

# making sure seq_count == 0 when detection == 0
sum(dat$seq_count[dat$detected==0])

# adding body size
pan <- read.table("../data/pantheria.txt", sep="\t", as.is=T, header=TRUE)
masses <- pan[,c("MSW05_Binomial","AdultBodyMass_g","MSW05_Family","MSW05_Order")]
names(masses) <- c("binomial","mass","family","order") 
masses$mass[masses$binomial=="Chlorocebus_pygerythrus"] <- masses$mass[masses$binomial=="Chlorocebus_aethiops"]

dat <- left_join(dat,masses)
unique(dat$binomial[is.na(dat$mass)])

View(dat)


####################################################
### eDNA detection without camera trap detection ###
####################################################

dat$seen[is.na(dat$seen)] <- 0

# 4 instances of detecting an animal with eDNA but not seeing it in the traps
dim(dat[dat$seen==0,])
sort(unique(dat$binomial[dat$seen==0]))
dat[dat$seen==0,] %>% arrange(binomial) %>% print(n=Inf)#

# KNP_VERT_MINBOOT_80
# Gyps africanus
# Numida meleagris
# Syncerus caffer (3 reads) - sample NYA_4
# Tragelaphus scriptus (34 reads) - sample KWA_6

# Syncerus caffer (3 reads) - sample NYA_4
trap %>% filter(sample_num=="NYA_4") %>% select(species)%>% unique()
# lion, hyena, jackal
trap %>% filter(sample_num=="NYA_4") %>% filter(species=="lion") %>% select(timestamp) %>% tail(1) - trap %>% filter(sample_num=="NYA_4") %>% select(timestamp) %>% tail(1)
# lion came 6.8 days before sampling...
trap %>% filter(sample_num=="NYA_4") %>% filter(species=="hyena") %>% select(timestamp) %>% tail(1) - trap %>% filter(sample_num=="NYA_4") %>% select(timestamp) %>% tail(1)
# hyena came 31 minutes before sampling - could be deposition
trap %>% filter(sample_num=="NYA_4") %>% filter(species=="jackal") %>% select(timestamp) %>% tail(1) - trap %>% filter(sample_num=="NYA_4") %>% select(timestamp) %>% tail(1)
# jackal came immediately before sampling - could be deposition


# Tragelaphus scriptus (34 reads) - sample KWA_6
trap %>% filter(sample_num=="KWA_6") %>% select(species)%>% unique()
# lion, hyena
trap %>% filter(sample_num=="KWA_6") %>% filter(species=="lion") %>% select(timestamp) %>% tail(1) - trap %>% filter(sample_num=="KWA_6") %>% select(timestamp) %>% tail(1)
# came immediately before sampling
trap %>% filter(sample_num=="KWA_6") %>% filter(species=="hyena") %>% select(timestamp) %>% tail(1) - trap %>% filter(sample_num=="KWA_6") %>% select(timestamp) %>% tail(1)
# $ came 5 hours before sampling


#####################
# Overall detection #
#####################

dat$seen[is.na(dat$seen)] <- 0

head(dat)

d_summ <- dat %>% group_by(binomial) %>% summarise(detections=sum(detected),seens=sum(seen), rate = detections/seens)
d_summ %>% print(n=23)
d_summ <- d_summ %>% filter(!is.infinite(rate)) %>% arrange(-rate)

# Not correct because sometimes you get detection as 1 and seen as 0...
# need to filter to remove these first - or move to separate separate table...

detected_not_seen <- dat %>% filter(detected==1 && seen==0)
# detected_not_seen$binomial

dat_seen <- dat %>% filter(seen==1)
d_summ <- dat_seen %>% group_by(binomial) %>% summarise(detections=sum(detected),seens=sum(seen), rate = detections/seens)
d_summ %>% print(n=23)
d_summ <- d_summ %>% filter(!is.infinite(rate)) %>% arrange(-rate, -seens)
d_summ %>% print(n=21)

names(d_summ) <- c("Species", "eDNA detections","Camera detections","Rate")

print(xtable(d_summ, type = "latex",caption="Per species rates of eDNA and camera detections by sample",
                        digits=c(0,0,0,0,2)), file = "../plots_tables/detect_rates_species.tex",include.rownames=FALSE)


# Detections are the same for ps_80 and ps_98
# but they dramatically increase with ps_50...

# Species detected by eDNA
sort(unique(dat$binomial[dat$detected==1]))

dat <- dat_seen

sum(dat$detected)# 26
nrow(dat)# 017
sum(dat$detected)/nrow(dat)# 0.243 of species x sample combinations were detected

#subset to last 48 hours
dat_48h <- dat[dat$last_visit<48 | is.na(dat$last_visit),]
sum(dat_48h$detected)# 25
nrow(dat_48h)# 81
sum(dat_48h$detected)/nrow(dat_48h)# 0.31 of species x sample compbinations were detected

d_summ <- dat_48h %>% group_by(binomial) %>% summarise(detections=sum(detected),seens=sum(seen), rate = detections/seens) %>% filter(!is.na(seens))%>% arrange(-rate, -seens)
d_summ

#subset to last 36 hours
dat_36h <- dat[dat$last_visit<36 | is.na(dat$last_visit),]
sum(dat_36h$detected)# 25
nrow(dat_36h)# 74
sum(dat_36h$detected)/nrow(dat_36h)# 0.34 of species x sample compbinations were detected

d_summ <- dat_36h %>% group_by(binomial) %>% summarise(detections=sum(detected),seens=sum(seen), rate = detections/seens) %>% filter(!is.na(seens))%>% arrange(-rate, -seens)
d_summ

#subset to last 24 hours
dat_24h <- dat[dat$last_visit<24 | is.na(dat$last_visit),]
sum(dat_24h$detected)# 23
nrow(dat_24h)# 66
sum(dat_24h$detected)/nrow(dat_24h)# 0.348 of species x sample compbinations were detected

d_summ <- dat_24h %>% group_by(binomial) %>% summarise(detections=sum(detected),seens=sum(seen), rate = detections/seens) %>% filter(!is.na(seens))%>% arrange(-rate, -seens)
d_summ

#subset to last 12 hours
dat_12h <- dat[dat$last_visit<12 | is.na(dat$last_visit),]
sum(dat_12h$detected)# 15
nrow(dat_12h)# 30
sum(dat_12h$detected)/nrow(dat_12h)# 0.5 of species x sample compbinations were detected

d_summ <- dat_12h %>% group_by(binomial) %>% summarise(detections=sum(detected),seens=sum(seen), rate = detections/seens) %>% filter(!is.na(seens)) %>% arrange(-rate, -seens)
d_summ
names(d_summ) <- c("Species", "eDNA detections","Camera detections","Rate")

print(xtable(d_summ, type = "latex",caption="Per species rates of eDNA and camera detections by sample, subset to camera detections within 12 hours prior to samlping",
                        digits=c(0,0,0,0,2)), file = "../plots_tables/detect_rates_species_12h.tex",include.rownames=FALSE)


#subset to last 6 hours
dat_6h <- dat[dat$last_visit<6 | is.na(dat$last_visit),]
sum(dat_6h$detected)# 9
nrow(dat_6h)# 18
sum(dat_6h$detected)/nrow(dat_6h)# 0.5 of species x sample compbinations were detected

d_summ <- dat_6h %>% group_by(binomial) %>% summarise(detections=sum(detected),seens=sum(seen), rate = detections/seens) %>% filter(!is.na(seens))%>% arrange(-rate, -seens)
d_summ
colSums(d_summ[,2:4])

#subset to last 3 hours
dat_3h <- dat[dat$last_visit<3 | is.na(dat$last_visit),]
sum(dat_3h$detected)# 7
nrow(dat_3h)# 14
sum(dat_3h$detected)/nrow(dat_3h)# 0.5 of species x sample compbinations were detected

d_summ <- dat_3h %>% group_by(binomial) %>% summarise(detections=sum(detected),seens=sum(seen), rate = detections/seens) %>% filter(!is.na(seens))%>% arrange(-rate, -seens)
d_summ


############
## Figure ##
############

plot_breaks <- c(3,6,12,24,48,96,168)

p_dat <- NULL

for( i in seq_along(plot_breaks)) {

  p_dat$plot_break[i] <- plot_breaks[i]
  df_small <- dat[dat$last_visit < plot_breaks[i] | is.na(dat$last_visit),]
  p_dat$edna_detect[i] <- sum(df_small$detected)/nrow(df_small)
}
p_dat <- as.data.frame(p_dat)


# calc site by species richness of eDNA and camera traps in the same way
p2_dat <- NULL

for( i in seq_along(plot_breaks)) {

  df_small <- dat[dat$last_visit < plot_breaks[i] | is.na(dat$last_visit),]
  df_small <- df_small %>% group_by(sample_num) %>% 
  transmute(edna_rich = sum(detected, na.rm=T), cam_rich=sum(seen, na.rm=T), prop=edna_rich/cam_rich)
  df_small$hour <- plot_breaks[i]
  p2_dat <- rbind(p2_dat, unique(df_small))
  
}
p2_dat
with(p2_dat, plot(prop ~ hour))

# plotting
names(p2_dat) <- c("sample_num","eDNA","Camera traps","prop","hour")

p3_dat <- gather(p2_dat, variable, value, -sample_num, -hour)
p3_dat
names(p3_dat) <- c("sample_num","hour", "Method","value")

p3_dat <- p3_dat %>% group_by(hour, Method) %>% 
        mutate(mean=mean(value), var=var(value), max=max(value), min=min(value))

p3_dat <- p3_dat[p3_dat$Method!="prop",]

colors <- c("#E7B800", "#00AFBB", "#FC4E07")

comparison_plot <- ggplot(p3_dat, aes(x=hour, y=mean, colour=Method)) + 
    geom_ribbon(aes(ymin=min,ymax=max,fill=Method),alpha=0.4) +
    geom_line() +
    geom_point() +
    scale_x_continuous(name="sampling window (hours)", 
                        # breaks=plot_breaks, trans="none") +
                        breaks=plot_breaks, trans="sqrt")+
    scale_y_continuous(name="Species richness", 
                        breaks=seq(0,16,by=2), expand = c(0, 0)) +
    scale_color_manual(values = colors)+
    scale_fill_manual(values = colors) +
    theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))

comparison_plot

ggsave("plots_tables/detection_by_method.pdf", plot = comparison_plot, device = "pdf",
  scale = 1, width = 6, height = 4, units = "in",
  dpi = 300, limitsize = TRUE)


##########################################################


# Subset to only entries where animals were seen
dat <- dat[dat$seen==1,]

require(ape)
tree <- read.tree("../data/mammals.tre")
unique(dat$binomial[!dat$binomial%in%tree$tip.label])
dat$binomial[dat$binomial=="Equus_quagga"] <- "Equus_burchellii"
dat$binomial[!dat$binomial%in%tree$tip.label]#0

# Raw Data Plots
blank_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), 
                     axis.title=element_text(size=14),
                     legend.text=element_text(size=14))

blank_theme_large <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=18), 
                     axis.title=element_text(size=18),
                     legend.text=element_text(size=18))

# Getting sum of visits over whole trap set
sum_species <- dat %>% 
                group_by(binomial) %>%
                summarise(sum_waterhole=sum(n_waterhole))

# Number of photos used in this study
with(trap, length(unique(paste(site,date_range,camera,folder,filename))))#11,138 photos

# Number of individuals identified to species level
with(trap, sum(n_present))#55,628 individuals
with(trap, sum(n_waterhole))#30,694 individuals at waterhole
with(trap, sum(n_contact))#15,011 individuals in contact

sum_species
sum(sum_species$sum_waterhole)
unique(sum_species$binomial)

# Identifying missing based on sequences
missed <-setdiff(sum_species$binomial,dat$binomial[dat$detected==1])
sum_species$missed <- as.factor(sum_species$binomial%in%missed*1)
sum_species$binomial <- gsub("_"," ", sum_species$binomial)

# Reorder by descending and make binomial a factor
sum_species <- sum_species %>% 
    arrange((sum_waterhole)) %>%
    mutate(binomial = factor(binomial,binomial))

head(sum_species)


total <- ggplot(sum_species, aes(x = binomial, y=(sum_waterhole+0.03),
        fill=missed)) + 
        scale_fill_brewer(palette = "Set2", 
            breaks=c("0", "1"), labels=c("Detected by eDNA", "Missed by eDNA"))+
        geom_bar(stat="identity",colour=FALSE) +
        blank_theme_large + 
        theme(
            axis.title.y=element_blank(), legend.title=element_blank(),
            legend.text=element_text(size=18)) +
        scale_y_log10(breaks = c(0,1,2, 5,20,50,100,250,750,2000,5000,13000),
            expand = c(0,0), limits=c(0.9,20000)) +
        expand_limits(y = 0) +
        ylab("\n Total number of individuals in all photos") +
        theme(legend.position=c(.80, .25)) +
        coord_flip()
total 

ggsave("plots_tables/total_individuals_detected_barplot_large.pdf", plot = total, device = "pdf",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)

ggsave("../plots_tables/total_individuals_detected_barplot_large.png", plot = total, device = "png",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)

total_bw <- ggplot(sum_species, aes(x = binomial, y=(sum_waterhole+0.03), fill=factor(0))) + 
        geom_bar(stat="identity",colour=FALSE) +
        blank_theme_large + 
        scale_fill_brewer(palette = "Set2", guide=FALSE)+
        theme(
            axis.title.y=element_blank()) +
        scale_y_log10(breaks = c(0,1,2, 5,20,50,100,250,750,2000,5000,13000),
            expand = c(0,0), limits=c(0.9,20000)) +
        expand_limits(y = 0) +
        ylab("\n Total number of individuals in all photos") +
        theme(legend.position=c(.80, .25)) +
        coord_flip()
total_bw 

ggsave("../plots_tables/total_individuals_barplot_large.pdf", plot = total_bw, device = "pdf",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)

ggsave("../plots_tables/total_individuals_barplot_large.png", plot = total_bw, device = "png",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)



# Removing species what we do not have barcode sequences for
dat <- dat[!dat$binomial%in%c("Mungos_mungo","Paraxerus_cepapi"),]
sort(unique(dat$binomial))
sort(unique(dat$sample_num))


r1 <- ggplot(dat, aes(mass/1000, n_waterhole)) + 
    blank_theme_large + geom_point(alpha = 1, size=3, 
        aes(shape=factor(detected))) +
    # scale_shape_identity()+ 
    scale_shape_manual(values = c(1,19), name="",labels=c("Missed by eDNA","Detected by eDNA"))+
    scale_x_continuous(trans = "log", breaks = c(5,10,20,50,100,250,500,1000,2000,4000)) +
    scale_y_continuous(trans = "log", breaks = c(1, 2,5,15,35,100,200,500,1000,2500,5000))+ 
    # scale_x_continuous(expand = c(0, 0.01), trans = "log", breaks = c(0,1,2,6,12,24,36,48,72,48*2,48*3)) +
    # scale_y_continuous(expand = c(0, 0.006), breaks = c(0,0.2,0.4,0.6,0.8,1.0)) + 
    xlab("\n Mass (kg)") + ylab("Visitation per species per sample\n")+
    geom_vline(xintercept=50,show.legend = FALSE,lty=2) +
    theme(legend.position=c(0.16,0.80)) +
    annotate("text", x = 30, y = 5500, label = "<50 kg", size=5)+
    annotate("text", x = 80, y = 5500, label = ">50 kg", size=5)
r1
ggsave("../plots_tables/visitation_vs_mass_detected_large.pdf", plot = r1, device = "pdf",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)

ggsave("../plots_tables/visitation_vs_mass_detected_large.png", plot = r1, device = "png",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)


dat$common_name <-NA
dat$common_name[dat$binomial=="Hippopotamus_amphibius"] <- "Hippo"

r2 <- ggplot(dat, aes(last_visit, n_waterhole, label=common_name)) + 
    blank_theme_large + geom_point(alpha = 1, size=3, 
        aes(shape=factor(detected))) +
    # scale_shape_identity()+ 
    scale_shape_manual(values = c(1,19), name="",labels=c("Missed by eDNA","Detected by eDNA"))+
    scale_x_continuous(trans = "sqrt", breaks = c(1,2,6,12,24,36,48,72,48*2,48*3)) +
    scale_y_continuous(trans = "log", breaks = c(1, 2,5,15,35,100,200,500,1000,2500,5000))+ 
    # scale_x_continuous(expand = c(0, 0.01), trans = "log", breaks = c(0,1,2,6,12,24,36,48,72,48*2,48*3)) +
    # scale_y_continuous(expand = c(0, 0.006), breaks = c(0,0.2,0.4,0.6,0.8,1.0)) + 
    xlab("\n Last visit (hours)") + ylab("Visitation per species per sample\n")+
    geom_vline(xintercept=40,show.legend = FALSE,lty=2) +
    theme(legend.position=c(0.75,0.80)) +
    annotate("text", x = 28, y = 5500, label = "<40 hrs", size=6)+
    annotate("text", x = 55, y = 5500, label = ">40 hrs", size=6)
r2 <- r2 + geom_text(aes(label=ifelse(last_visit>40 & detected==1,as.character(common_name),'')),hjust=0,vjust=-1, size=5)
r2

dat[dat$last_visit>40 & dat$detected==1,]
# Hippos in NYA_3 (collected july1) apparently last_visit was 106 hours before...
# NYA is pipeline trough, so perhaps coming in from olifants?

ggsave("../plots_tables/visitation_vs_last_visit_detected_large.pdf", plot = r2, device = "pdf",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)

ggsave("../plots_tables/visitation_vs_last_visit_detected_large.png", plot = r2, device = "png",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)


##########################
## TIME SERIES ANALYSES ##
##########################

# Plotting time series

require(scales);packageVersion("scales")# 1.0.0

p_by_species$species
sp <- p_by_species$species[c(1,2,4,16,14)]

# For KWA let's try these
sp <- p_by_species$species[c(1,2,14,4,5)]

trap$binomial[trap$binomial=="Equus_quagga"] <- "Equus_burchellii"

trap_sub <- trap[trap$binomial %in% sp,]
trap_sub <- trap_sub[trap_sub$sample_num %in% camera_samples_fullweek,]

# Subset to single site
trap_sub <- trap_sub[trap_sub$site == "KWA",]

# adding detection
trap_sub <- left_join(trap_sub, dat[,c("sample_num","binomial","detected")])

trap_sub$binomial <- factor(trap_sub$binomial, levels = sp)
# renaming to common names
levels(trap_sub$binomial) <- c("White Rhino", "Hyena", "Lion", "Elephant", "Zebra")

trap_sub$site <- as.factor(trap_sub$site)
trap_sub$detected <- as.factor(trap_sub$detected)

require(RColorBrewer);packageVersion("RColorBrewer")# 1.1.2

trap_sub <- trap_sub[trap_sub$n_waterhole > 0,]

ts <- ggplot(trap_sub, aes(x=timestamp, y=n_waterhole, shape=site)) + 
      geom_bar(aes(color=detected, fill=detected), stat="identity") +
      xlab("") + ylab("Abundance at waterhole") + 
      theme_minimal() + 
      guides(color=FALSE, fill=FALSE) + 
      scale_x_datetime(labels = date_format("%a %Hh"), date_breaks="6 hours") +
      facet_grid(binomial ~ ., switch="y") + 
      theme(strip.text.y = element_text(angle = 180), axis.text.x=element_text(angle=90)) + 
      scale_y_continuous(position = "right", breaks= pretty_breaks()) + 
      scale_color_manual(values = c("#FC8D62","#66C2A5")) + 
      geom_vline(aes(xintercept=collected$collected[collected$sample_num=="KWA_6"]), color="red", linetype="dashed")
ts

ggsave("../plots_tables/timeseries_KWA.pdf", plot = ts, device = "pdf",
  scale = 1, width = 10, height = 6, units = "in",
  dpi = 300, limitsize = TRUE)


dat %>% filter(site=="NYA") %>% 
        filter(detected==1) %>%
        print(n=Inf)

dat %>% filter(site=="KWA") %>% 
        filter(detected==1) %>%
        print(n=Inf)

unique(trap$binomial[trap$sample_num=="KWA_6"])




