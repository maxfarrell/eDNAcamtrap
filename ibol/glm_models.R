# eDNA camtrap
# GLM MODELS

rm(list=ls())

require(dplyr)
require(tidyr)

# Using merged edna camtrap data from july 16th 2017 (after all samples sequenced)
load("merged_eDNA_camtrap_data_july16_2017.RData")

# Using merged edna camtrap data from december 14th 2016 (before samples sequenced)
# load("../../Data/merged_eDNA_camtrap_data_dec14.RData")

sort(unique(seqs[[1]]$sample_num))
# There are some samples in seqs that are not in trap
seqs <- lapply(seqs, function(x) x[x$sample_num%in%trap$sample_num,])

cutoff_values <- seq(90,100, by=0.05)

seqXspec_cutoffs <- lapply(seqs, function(x) unique(x[,c("qseqid","species")]))
discrep <- lapply(seqXspec_cutoffs, function(x) x$qseqid[duplicated(x$qseqid)])
nodiscrep <-  lapply(seqXspec_cutoffs, function(x) x$qseqid[!duplicated(x$qseqid)])
nodiscrep_spec <-  lapply(seqXspec_cutoffs, function(x) unique(x$species[!duplicated(x$qseqid)]))

n_nodiscrep_seqs <- unlist(lapply(nodiscrep, length))
n_discrep_seqs <- unlist(lapply(discrep, length))
nodiscrep_richness <- unlist(lapply(nodiscrep_spec, length))


# Make sure species richness is correct sample numbers and seqs have no discrep...
sort(nodiscrep_spec[[which(cutoff_values==96)]])
sort(nodiscrep_spec[[which(cutoff_values==97)]])
# Lose Tragelaphus_scriptus

eDNA97 <- seqs[[which(cutoff_values==97)]]
eDNA97 <- eDNA97[eDNA97$sample_num%in%trap$sample_num,]
eDNA97 <- eDNA97[eDNA97$qseqid %in% nodiscrep[[which(cutoff_values==97)]],]

str(eDNA97)
unique(eDNA97$site[eDNA97$species=="Tragelaphus_scriptus"])

str(trap)
range(trap$n_waterhole)

dat1 <- trap %>% 
		filter(n_waterhole!=0)%>%
		group_by(sample_num, binomial)%>%
		summarize(n_waterhole=sum(n_waterhole), 
			seen=((n_waterhole>0)*1))

dat2 <- eDNA97 %>%
		group_by(sample_num, species)%>%
		summarise(detected=(length(qseqid)>1)*1)
names(dat2)[2] <- "binomial"

dat <- full_join(dat1,dat2)
dat$detected[is.na(dat$detected)] <- 0
dat <- dat[dat$n_waterhole>0,]
dat <- dat[dat$binomial!="unknown",]
dat <- dat[!is.na(dat$binomial),]

pan <- read.table("../../Data/pantheria.txt", sep="\t", as.is=T, header=TRUE)
masses <- pan[,c("MSW05_Binomial","AdultBodyMass_g","MSW05_Order")]
names(masses) <- c("binomial","mass","order") 
masses$mass[masses$binomial=="Chlorocebus_pygerythrus"] <- masses$mass[masses$binomial=="Chlorocebus_aethiops"]

dat <- left_join(dat,masses)
dat$binomial[is.na(dat$mass)]

require(ape)
tree <- read.tree("../../Data/mammals.tre")

unique(dat$binomial[!dat$binomial%in%tree$tip.label])
dat$binomial[dat$binomial=="Equus_quagga"] <- "Equus_burchellii"
dat$binomial[!dat$binomial%in%tree$tip.label]#0

# Removing species what we do not have barcode sequences for
dat <- dat[!dat$binomial%in%c("Mungos_mungo","Paraxerus_cepapi"),]
sort(unique(dat$binomial))
sort(unique(dat$sample_num))

(dat$detected*18)+1
# Plot of visitation ~ mass (line indicates min detected bodysize)
with(dat, summary(lm(log(n_waterhole) ~ log(mass))))

# pdf("../../Data/Plots/mass_visitation_plain.pdf")
# jpeg("../../Data/Plots/mass_visitation_plain.jpg", width = 7, height = 7, units = 'in', res = 600)
# with(dat, plot(log(n_waterhole) ~ log(mass), pch=1,
# 		ylab="log( Visitation par site par espèce )",
# 		xlab="log( Masse d'espèce )"))
# dev.off()

# # pdf("../../Data/Plots/mass_visitation_detected.pdf")
# jpeg("../../Data/Plots/mass_visitation_detected.jpg", width = 7, height = 7, units = 'in', res = 600)
# with(dat, plot(log(n_waterhole) ~ log(mass), pch=(detected*18)+1,
# 		ylab="log( Visitation par site par espèce )",
# 		xlab="log( Masse d'espèce )"))
# legend(x=8.5,y=8, legend=c("Manqué","Détecté"),  pch=c(1,19))
# dev.off()

# # pdf("../../Data/Plots/mass_visitation_detected_mass.pdf")
# jpeg("../../Data/Plots/mass_visitation_detected_mass.jpg", width = 7, height = 7, units = 'in', res = 600)
# with(dat, plot(log(n_waterhole) ~ log(mass), pch=(detected*18)+1,
# 		ylab="log( Visitation par site par espèce )",
# 		xlab="log( Masse d'espèce )"))
# legend(x=8.5,y=8, legend=c("Manqué","Détecté"),  pch=c(1,19))
# abline(v=min(log(dat$mass[dat$detected==1]))-0.1, lty=2)
# text(10,8.2, "<53 kg")
# text(11.8,8.2, ">53 kg")
# dev.off()

# with(dat, text(log(n_waterhole) ~ log(mass), labels=binomial, cex= 0.7, pos=1))

exp(min(log(dat$mass[dat$detected==1])))/1000
exp(min(log(dat$mass[dat$detected==1 & dat$n_waterhole>2])))/1000

log(unique(dat$mass[dat$binomial=="Papio_ursinus"]))

# but next biggest is papio ursiunus
# ~18kg

# jpeg("../../Data/Plots/mass_visitation_detected_mass_papio.jpg", width = 7, height = 7, units = 'in', res = 600)
# with(dat, plot(log(n_waterhole) ~ log(mass), pch=(detected*18)+1,
# 		ylab="log( Visitation par site par espèce )",
# 		xlab="log( Masse d'espèce )"))
# legend(x=8.5,y=8, legend=c("Manqué","Détecté"),  pch=c(1,19))
# abline(v=log(unique(dat$mass[dat$binomial=="Papio_ursinus"]))+0.1, lty=2)
# text(9,6, "< 18 kg")
# dev.off()



dat[dat$n_waterhole==1,]
dat[dat$n_waterhole==1 & dat$detected==1,]

dat3<- dat %>%
		filter(detected==1)%>%
		group_by(binomial)%>%
		mutate(min_n_detected = min(n_waterhole),mass==mass)

dat3 <- unique(dat3[,c("binomial","min_n_detected","mass")])
# with(dat3, plot(log(min_n_detected) ~ log(mass)))

spec_never_detected <- dat%>%
						group_by(binomial)%>%
						summarize(ever_detected=sum(detected))
spec_never_detected <- spec_never_detected$binomial[spec_never_detected$ever_detected==0]

dat4<- dat %>%
		filter(binomial%in%spec_never_detected)%>%
		filter(detected==0)%>%
		group_by(binomial)%>%
		mutate(max_n_undetected = max(n_waterhole),mass==mass)

dat4 <- unique(dat4[,c("binomial","max_n_undetected","mass")])
# with(dat4, plot(log(max_n_undetected) ~ log(mass)))
# abline(v=min(log(dat$mass[dat$detected==1])))
# with(dat4, text(log(max_n_undetected) ~ log(mass), labels=binomial, cex= 0.7, pos=3))

dat4 %>% ungroup() %>% arrange(mass)


library(ape)
library(geiger)
library(nlme)
library(phytools)

hist(log(dat$n_waterhole))
hist(log(dat$mass))

test <- dat%>%
		group_by(binomial)%>%
		summarise(sd_waterhole = sd(n_waterhole),
			n_waterhole = mean(n_waterhole),
			mass=mean(mass))
test

test$sd_waterhole[is.na(test$sd_waterhole)] <- 0

with (test, plot(log(n_waterhole) ~ log(mass)))
with (test, plot(log(n_waterhole) ~ log(mass), cex=log(sd_waterhole+2)))
# redo with error bars

geom_errorbar(aes(ymin=len-se, ymax=len+se), width=.1, position=pd)
limits <- aes(ymax = n_waterhole + sd_waterhole, ymin=n_waterhole - sd_waterhole)

ggplot(test, aes(x=mass/1000, y=n_waterhole)) + geom_point() + 
	# geom_pointrange(aes(ymax=n_waterhole+sd_waterhole, ymin=n_waterhole-sd_waterhole)) + 
	scale_x_log10() + scale_y_log10() +
	geom_vline(xintercept = min((dat$mass[dat$detected==1])/1000)-2.5)

exp(min(log(dat$mass[dat$detected==1])))/1000

exp(min(log(dat$mass[dat$binomial=="Papio_ursinus"])))/1000
exp(min(log(dat$mass[dat$binomial=="Panthera_pardus"])))/1000

hist(log(unique(dat$mass)), breaks=10)
abline(v=min(log(dat$mass[dat$detected==1])), col=2)


# CALCULATING TIME SINCE LAST VISIT

str(trap)
dat1 <- trap %>% 
		group_by(sample_num)%>%
		mutate(lastphoto = max(timestamp))%>%
		filter(n_waterhole!=0)%>%
		group_by(sample_num, binomial)%>%
		summarize(n_waterhole=sum(n_waterhole), 
			seen=((n_waterhole>0)*1),
			lastseen=as.numeric(unique(lastphoto)-max(timestamp), unit="hours"))
dat2 <- eDNA97 %>%
		group_by(sample_num, species)%>%
		summarise(detected=(length(qseqid)>1)*1)
names(dat2)[2] <- "binomial"

dat <- full_join(dat1,dat2)
dat$detected[is.na(dat$detected)] <- 0
dat <- dat[dat$n_waterhole>0,]
dat <- dat[dat$binomial!="unknown",]
dat <- dat[!is.na(dat$binomial),]

pan <- read.table("../../Data/pantheria.txt", sep="\t", as.is=T, header=TRUE)
masses <- pan[,c("MSW05_Binomial","AdultBodyMass_g","MSW05_Order")]
names(masses) <- c("binomial","mass","order") 
masses$mass[masses$binomial=="Chlorocebus_pygerythrus"] <- masses$mass[masses$binomial=="Chlorocebus_aethiops"]

dat <- left_join(dat,masses)
dat$binomial[is.na(dat$mass)]

with(dat, plot(detected ~ lastseen))
with(dat, plot((lastseen) , log(mass), ))


# jpeg("../../Data/Plots/lastseen_visitation_detected.jpg", width = 7, height = 7, units = 'in', res = 600)
# with(dat, plot((lastseen) , log(n_waterhole), pch=(detected*18)+1,
# 		ylab="log( Visitation par site par espèce )",
# 		xlab="Dernière observation d'une espèce avant l'échantillonnage (heures)"))
# legend(x=120,y=8, legend=c("Manqué","Détecté"),  pch=c(1,19))
# abline(v=max((dat$lastseen[dat$detected==1]))+1.5, lty=2)
# text(55,8, ">36 heures")
# dev.off()

max((dat$lastseen[dat$detected==1]))

with(dat, plot((lastseen) , log(n_waterhole), pch=(detected*18)+1, cex=log10((mass)/1000)))
# with(dat, text(lastseen , log(n_waterhole), labels=binomial, cex= 0.7, pos=3))


max((dat$lastseen[dat$detected==1]))
# 104 hours
max((dat$lastseen[dat$detected==1]))/24
# 4.3 days

dat$binomial[dat$lastseen==max(dat$lastseen[dat$detected==1])]
# Hippo (96)
# Rhino

sort(unique(dat$lastseen[dat$detected==1]))
length(sort(unique(dat$lastseen[dat$detected==1])))-1

dat$binomial[dat$lastseen==sort(unique(dat$lastseen[dat$detected==1]))[28]]
dat$binomial[dat$lastseen==sort(unique(dat$lastseen[dat$detected==1]))[28]]

unique(dat$lastseen)
str(dat)


# with(dat, plot((lastseen) , log(n_waterhole), pch=(detected)*19))
# with(dat, plot((lastseen) , log(n_waterhole), pch=(detected)*19, cex=log10(mass/100)))
# abline(v=max((dat$lastseen[dat$detected==1]))-0.1)
# abline(v=median((dat$lastseen[dat$detected==1]))-0.1)
# with(dat, text(lastseen , log(n_waterhole), labels=binomial, cex= 0.7, pos=3))


## CALCULATING TOTAL EDNA
str(dat)
dat$edna <- dat$n_waterhole*dat$mass
with(dat, plot((lastseen) , (edna), pch=(detected)*19))
with(dat, plot((lastseen) , log(edna), pch=(detected)*19))
with(dat, plot((lastseen) , log(edna), pch=(detected)*19, cex=log10(mass/100)))

edna1 <- trap %>% 
		group_by(sample_num)%>%
		mutate(lastphoto = max(timestamp),
			time=as.numeric(lastphoto-timestamp,units="secs"))%>%
		filter(n_waterhole!=0)%>%
		group_by(sample_num, binomial)%>%
		mutate(edna=n_waterhole*exp(-time/(60*24*3)))%>%
		summarise(edna=sum(edna),
			n_waterhole=sum(n_waterhole),
			seen=((n_waterhole>0)*1),
			lastseen=as.numeric(unique(lastphoto)-max(timestamp), unit="hours"))

edna2 <- eDNA97 %>%
		group_by(sample_num, species)%>%
		summarise(detected=(length(qseqid)>1)*1)
names(edna2)[2] <- "binomial"

edna <- full_join(edna1,edna2)
edna$detected[is.na(edna$detected)] <- 0
edna <- edna[edna$n_waterhole>0,]
edna <- edna[edna$binomial!="unknown",]
edna <- edna[!is.na(edna$binomial),]
edna <- left_join(edna,masses)
edna$binomial[is.na(dat$mass)]
edna$total_edna <- edna$edna * edna$mass

with(edna, plot((lastseen) , log(total_edna), pch=(detected)*19, cex=(log10(mass/10000))))
with(edna, plot(log(n_waterhole) , log(total_edna), pch=(detected)*19, cex=(log10(mass/10000))))

with(edna, hist(log(total_edna*mass)))
abline(v=min(edna$total_edna[dat$detected==1]))

with(edna, plot(as.factor(detected), log(total_edna), pch=(detected)*19, cex=(log10(mass/10000))))
with(edna, summary(glm(detected ~ log(total_edna))))



# Calculating how many elephants + other species came to drink

ele1 <- trap %>% 
		group_by(sample_num)%>%
		mutate(lastphoto = max(timestamp),
			time=as.numeric(lastphoto-timestamp,units="secs"),
			elephants= sum(n_waterhole[binomial=="Loxodonta_africana"]),
			all_visits=sum(n_waterhole))%>%
		filter(n_waterhole!=0)%>%
		group_by(sample_num, binomial)%>%
		mutate(edna=n_waterhole*exp(-time/(60*24*3)))%>%
		summarise(edna=sum(edna),
			n_waterhole=sum(n_waterhole),
			seen=((n_waterhole>0)*1),
			lastseen=as.numeric(unique(lastphoto)-max(timestamp), unit="hours"),
			elephants=unique(elephants), all_visits=unique(all_visits))

ele <- full_join(ele1,edna2)
ele$detected[is.na(ele$detected)] <- 0
ele <- ele[ele$n_waterhole>0,]
ele <- ele[ele$binomial!="unknown",]
ele <- ele[!is.na(ele$binomial),]
ele <- left_join(ele,masses)
ele$binomial[is.na(dat$mass)]
ele$total_ele <- ele$ele * ele$mass

# with(ele, plot(as.factor(detected), log(elephants), pch=(detected)*19, cex=(log10(mass/10000))))
# with(ele, plot(as.factor(detected), log(total_ele), pch=(detected)*19, cex=(log10(mass/10000))))
# with(ele, plot(as.factor(detected), log(all_visits), pch=(detected)*19, cex=(log10(mass/10000))))
# with(ele, plot(as.factor(detected), log(all_visits*mass), pch=(detected)*19, cex=(log10(mass/10000))))
# means that as visitation of all species increases, detection decreases
# OBVIOUSLY CONFOUNDING MYSELF

# MAYBE NEED TO MEASURE ELEPHANT VISITATION 
# AFTER VISITATION OF OTHER SPECIES
# TOTAL VISITATION DOES NOT MATTER, 
# ONLY THE AMOUNT AFTER A PARTICULAR VISITATION EVENT


# MODELS

# BAYESIAN MODEL AVERAGING

require("BMA")
str(ele)

require(car)
fit <- glm(detected ~ log2(n_waterhole) + log2(mass) + lastseen, data=ele, family=binomial(link="logit"))
vif(fit)

# Removing species what we do not have barcode sequences for
ele <- ele[!ele$binomial%in%c("Mungos_mungo","Paraxerus_cepapi"),]
sort(unique(ele$binomial))

predictors <- c("n_waterhole","mass", "lastseen")
# predictors <- c("mass","n_waterhole", "lastseen")
dat.x <- data.frame(ele[,predictors])
dat.x[,predictors[predictors!="lastseen"]] <- log2(dat.x[,predictors[predictors!="lastseen"]])
head(dat.x)

str(dat.x)

glm.out <- bic.glm(x=dat.x, y=ele$detected, strict = FALSE, OR = 200,
glm.family="binomial", factor.type=FALSE)
summary(glm.out)
glm.out$postmean

# glm.out$condpostmean

# to help interpretation:
# http://www.ats.ucla.edu/stat/mult_pkg/faq/general/odds_ratio.htm


# LOR of mass
glm.out$postmean[3]

# Can use inverse logit to transform logodds back to probablity
inv_logit <- function(p) {exp(p)/(1+exp(p))}

# On probability scale
inv_logit(glm.out$postmean[3])

# inv_logit of intercept is probability of detection when all else is zero
# should be zero, but is 0.07 meaning fitting is not perfect
inv_logit(glm.out$postmean[1])

# n_waterhole - for a doubling in mass (becuase we scaled mass as log2)
inv_logit(2^(glm.out$postmean[2]))
exp(2^(glm.out$postmean[2]))
# Why is this so high?

# On probability scale - for a doubling in mass (becuase we scaled mass as log2)
inv_logit(2^(glm.out$postmean[3]))
exp(2^(glm.out$postmean[3]))

# for Lastseen (On probability scale)
inv_logit((glm.out$postmean[4]))
exp((glm.out$postmean[4]))
# this makes sense
1-exp((glm.out$postmean[4]))
# ~8 % decrease in detection probability per hour

exp((glm.out$postsd[4]))
1-exp((glm.out$postsd[4]))


# jpeg("../../Data/Plots/postmean_plots.jpg", width = 12, height = 4, units = 'in', res = 600)
plot(glm.out, mfrow=c(1,3))
# dev.off()
(100-glm.out$probne0)/100


glm.out$postsd
# Removing species what we do not have barcode sequences for
ele <- ele[!ele$binomial%in%c("Mungos_mungo","Paraxerus_cepapi"),]
sort(unique(ele$binomial))

predictors <- c("mass", "lastseen")
# predictors <- c("mass","n_waterhole", "lastseen")
dat.x <- data.frame(ele[,predictors])
dat.x[,predictors[predictors!="lastseen"]] <- log2(dat.x[,predictors[predictors!="lastseen"]])
head(dat.x)

glm.out <- bic.glm(x=dat.x, y=ele$detected, strict = FALSE, OR = 200,
glm.family="binomial", factor.type=FALSE)
summary(glm.out)
glm.out$postmean
glm.out$postsd
plot(glm.out, mfrow=c(1,3))


# Claculating optimal threshold / ROC
# http://stackoverflow.com/questions/23240182/deciding-threshold-for-glm-logistic-regression-model-in-r
require(ROCR)

predict(glm.out, newdata=dat.x)

pred <- prediction(predict(glm.out, newdata=dat.x), labels=ele$detected)

str(pred)
pred@tp[[1]][round(pred@cutoffs[[1]], 3)==0.364]
sum(ele$detected)#34
sum(ele$detected==1)#34
length(ele$detected)#105
26/34



length(pred@tp[[1]])
length(pred@cutoffs[[1]])

plot(pred@tp[[1]]/max(pred@tp[[1]]) ~ pred@cutoffs[[1]], cex=0)
lines((1-(pred@tp[[1]])/max(pred@tp[[1]])) ~ pred@cutoffs[[1]], col="red", )
lines(pred@tp[[1]]/max(pred@tp[[1]]) ~ pred@cutoffs[[1]], col="blue")

perf <- performance(pred,"tpr","fpr")
plot(perf)

perf1 <- performance(pred, "prec", "rec")
plot(perf1)

perf2 <- performance(pred, "sens", "cutoff")
plot(perf2)

perf3 <- performance(pred, "spec","cutoff")
plot(perf3)
lines(((pred@tp[[1]])/max(pred@tp[[1]])) ~ pred@cutoffs[[1]], col="red", )
abline(v=0.365)




# Prediction - demonstrating where LOR comes from 
newdata <- with(ele, 
			data.frame(n_waterhole = mean(log(n_waterhole)), 
						mass=c(mean(log(mass)),mean(log(mass))+1),
						lastseen=mean(lastseen), elephants=mean(log(elephants))))

predict(glm.out, newdata=newdata)
# Prediction is on what scale?
# Probability?

# use logit to get back to logit scale
logit(predict(glm.out, newdata=newdata))
LOR <- logit(predict(glm.out, newdata=newdata)[[2]]) - logit(predict(glm.out, newdata=newdata)[[1]])
LOR
#YES! - got the answer back
inv_logit(2^LOR)
# There is approximately 77.8% increase in detection as body size doubles


newdata <- with(ele, 
			data.frame(n_waterhole = c(mean(log(n_waterhole)),mean(log(n_waterhole))+1), 
						mass=mean(log(mass)),
						lastseen=mean(lastseen), elephants=mean(log(elephants))))

predict(glm.out, newdata=newdata)
# Prediction is on what scale?
# Probability?

# use logit to get back to logit scale
logit(predict(glm.out, newdata=newdata))
LOR <- logit(predict(glm.out, newdata=newdata)[[2]]) - logit(predict(glm.out, newdata=newdata)[[1]])
LOR
#YES! very close answer... off by 0.00001

# to get back to odds, we exponentiate this
exp(2^LOR)

# There is approximately 77.8% increase in detection as body size doubles



plot(ele$detected, predict(glm.out, newdata=dat.x))


imageplot.bma(glm.out)



glm.out$se
glm.out$mle

str(glm.out)

View(dat)




visits <- trap %>%









# SINGLE GLMS

glm1 <- glm(detected ~ n_waterhole, data=dat, family=binomial(link="logit"))
summary(glm1)
# plot(glm1)

glm2 <- glm(detected ~ log(n_waterhole), data=dat, family=binomial(link="logit"))
summary(glm2)
# plot(glm2)
# plot(glm2$fitted,glm2$residuals)
pseudorR2 <- 1 - (glm2$deviance/glm2$null.deviance)# 0.09

glm3 <- glm(detected ~ log(n_waterhole) + log(mass), data=dat, family=binomial(link="logit"))
summary(glm3)
# plot(glm3)

glm4 <- glm(detected ~ log(mass), data=dat, family=binomial(link="logit"))
summary(glm4)
# plot(glm4)

glm5 <- glm(detected ~ order + log(mass) + log(n_waterhole), data=dat, family=binomial(link="logit"))
summary(glm5)




glm4 <- glm(detected ~ log(n_waterhole) * log(mass), data=dat, family=binomial(link="logit"))
summary(glm4)

## SCRAPS ###

# CALCULATING N_WATERHOLE AT LASTSEEN (per species per sample)
# dat1 <- trap %>% 
# 		group_by(sample_num)%>%
# 		mutate(lastphoto = max(timestamp))%>%
# 		filter(n_waterhole!=0)%>%
# 		group_by(sample_num, binomial)%>%
# 		summarize(lastseen=as.numeric(unique(lastphoto)-max(timestamp), unit="hours"),
# 			last_n = unique(n_waterhole[timestamp==max(timestamp)]),
# 			n_waterhole=sum(n_waterhole), 
# 			seen=((n_waterhole>0)*1))
# dat2 <- eDNA97 %>%
# 		group_by(sample_num, species)%>%
# 		summarise(detected=(length(qseqid)>1)*1)
# names(dat2)[2] <- "binomial"

# dat <- full_join(dat1,dat2)
# dat$detected[is.na(dat$detected)] <- 0
# dat <- dat[dat$n_waterhole>0,]
# dat <- dat[dat$binomial!="unknown",]
# dat <- dat[!is.na(dat$binomial),]

# pan <- read.table("../../Data/pantheria.txt", sep="\t", as.is=T, header=TRUE)
# masses <- pan[,c("MSW05_Binomial","AdultBodyMass_g","MSW05_Order")]
# names(masses) <- c("binomial","mass","order") 
# masses$mass[masses$binomial=="Chlorocebus_pygerythrus"] <- masses$mass[masses$binomial=="Chlorocebus_aethiops"]

# dat <- left_join(dat,masses)
# dat$binomial[is.na(dat$mass)]

# with(dat, plot(last_n ~ lastseen, pch=(detected)*19))
# with(dat, plot(last_n ~ lastseen, pch=(detected)*19, cex=log10(mass/100)))
# # not great...





# require(pez)

# phy <- drop.tip(tree, tree$tip.label[!tree$tip.label%in%dat$binomial])
# # plot(phy)
# nspp <- length(unique(dat$binomial))
# Vphy <- vcv(phy)
# Vphy <- Vphy/(det(Vphy)^(1/nspp))

# r.intercept.spp.phy <- list(1, sp = phy$tip.label, covar = Vphy)
# str(dat)

# model <- communityPGLMM(detected ~ n_waterhole, data=dat, sp=as.factor(dat$binomial), 
# 	site=as.factor(dat$sample_num), family = "binomial", 
# 	random.effects = r.intercept.spp.phy)

