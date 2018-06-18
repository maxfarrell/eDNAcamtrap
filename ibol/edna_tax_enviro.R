# Animal detection via waterhole eDNA
# Basic model
# MJ Farrell

# Last updated November 6th 2017

rm(list=ls())

require(dplyr)

require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# loading data
load("merged_eDNA_camtrap_data_nov6_2017.RData")

# Adding sample_num to env
labNotes$site_time <- with(labNotes, paste(site, date_collected, sep="_"))
env$site_time <- with(env, paste(site, date, sep="_"))
lookup <- setNames(unique(labNotes$sample_num), unique(labNotes$site_time))
env$sample_num <- lookup[env$site_time]

# Adding A/B to trap for merging with sequence data
# dim(trap)
# trap2 <- rbind(trap, trap)
# trap2$AB <- rep(c("A","B"),each=nrow(trap))

# including time of sample collection to calclate time since last visit for each species per sample
collected <- seqs %>% select(sample_num, collected) %>% unique()

dat1 <- trap %>% full_join(collected)%>%
        filter(n_waterhole!=0)%>%
        group_by(sample_num, binomial)%>%
        mutate(n_waterhole=sum(n_waterhole), 
            seen=((n_waterhole>0)*1),
            last_visit=as.numeric(unique(collected)-max(timestamp), unit="hours"))%>%
        # select(site,sample_num,AB,binomial,n_waterhole,last_visit,seen)%>%
        select(site,sample_num,binomial,n_waterhole,last_visit,seen)%>%
        unique()%>%
        arrange(site, sample_num, binomial)     
# View(dat1)

str(seqs)
dat2 <- seqs %>%
        filter(is.na(S_XS))%>%
        # group_by(sample_num, AB, species)%>%
        group_by(sample_num, species)%>%
        mutate(detected=1)%>%
        select(site, sample_num, species, detected) %>%
        unique() %>%
        arrange(site, sample_num, species)      
# View(dat2)
names(dat2)[names(dat2)=="species"]<- "binomial"

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
dat <- dat[!is.na(dat$binomial),]
# one instance of seeing an animal but not detecting it -
# dat[dat$seen==0,]
# Syncerus caffer in NYA_4
# Unsure what to to with lastseen in this case... set as NA currently
dat$n_waterhole[is.na(dat$n_waterhole)] <- 0
dat$seen[is.na(dat$seen)] <- 0
dat <- dat[dat$seen==1,]

pan <- read.table("../../Data/pantheria.txt", sep="\t", as.is=T, header=TRUE)
masses <- pan[,c("MSW05_Binomial","AdultBodyMass_g","MSW05_Family","MSW05_Order")]
names(masses) <- c("binomial","mass","family","order") 
masses$mass[masses$binomial=="Chlorocebus_pygerythrus"] <- masses$mass[masses$binomial=="Chlorocebus_aethiops"]

dat <- left_join(dat,masses)
dat$binomial[is.na(dat$mass)]

# require(ape)
# tree <- read.tree("../../Data/mammals.tre")

# unique(dat$binomial[!dat$binomial%in%tree$tip.label])
# dat$binomial[dat$binomial=="Equus_quagga"] <- "Equus_burchellii"
# dat$binomial[!dat$binomial%in%tree$tip.label]#0

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
     
ggsave("total_individuals_detected_barplot_large.pdf", plot = total, device = "pdf",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)

ggsave("total_individuals_detected_barplot_large.png", plot = total, device = "png",
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

ggsave("total_individuals_barplot_large.pdf", plot = total_bw, device = "pdf",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)

ggsave("total_individuals_barplot_large.png", plot = total_bw, device = "png",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)


# Removing species what we do not have barcode sequences for
dat <- dat[!dat$binomial%in%c("Mungos_mungo","Paraxerus_cepapi"),]
sort(unique(dat$binomial))
sort(unique(dat$sample_num))


# with(dat, plot(log(n_waterhole) ~ log(mass), pch=(detected*18)+1,
#       ylab="log (Visitation per sample per species)",
#       xlab="log (Species mass (g))"))
# legend(x=8.25,y=7.5, legend=c("Missed by eDNA","Detected by eDNA"),  pch=c(1,19))
# abline(v=min(log(dat$mass[dat$detected==1]))-0.1, lty=2)
# text(10,8.2, "<53 kg")
# text(11.5,8.2, ">53 kg")



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
ggsave("visitation_vs_mass_detected_large.pdf", plot = r1, device = "pdf",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)

ggsave("visitation_vs_mass_detected_large.png", plot = r1, device = "png",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)


# with(dat, plot((last_visit) , log(n_waterhole), pch=(detected*18)+1,
#       ylab="log( Visitation par site par espèce )",
#       xlab="Dernière observation d'une espèce avant l'échantillonnage (heures)"))
# legend(x=120,y=8, legend=c("Manqué","Détecté"),  pch=c(1,19))
# abline(v=max((dat$lastseen[dat$detected==1]))+1.5, lty=2)
# text(55,8, ">36 heures")

r2 <- ggplot(dat, aes(last_visit, n_waterhole)) + 
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
r2

ggsave("visitation_vs_last_visit_detected_large.pdf", plot = r2, device = "pdf",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)

ggsave("visitation_vs_last_visit_detected_large.png", plot = r2, device = "png",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)





## Creating stan.data
data <- dat

# Normalizing continuous predictors to mean = 0, sd = 0.5
scale_half <- function(x){
    return( (x-mean(x)) * 0.5/sd(x))
}

# with(dat, hist(log(n_waterhole)))
data$n_waterhole <- scale_half(log(data$n_waterhole))
# with(dat, hist(log(last_visit)))
# with(dat, hist(sqrt(last_visit)))
data$last_visit <- scale_half(sqrt(data$last_visit))
# with(dat, hist(log(mass)))
data$mass <- scale_half(log(data$mass))
# with(dat, hist((temp)))
# with(dat, hist((ms)))
# with(dat, hist((do_p)))
with(dat, hist((ph)))
names(dat)
pairs(dat[,c(4:5,8:14)])

stan.data <- with(data, list(
    N=nrow(data),
    M=length(unique(binomial)), 
    species=as.integer(as.factor(binomial)),
    F=length(unique(family)), 
    family=as.integer(as.factor(unique(cbind(family,binomial))[,1])), 
    O=length(unique(order)), 
    order=as.integer(as.factor(unique(cbind(order,family))[,1])), 
    mass=as.numeric(unique(cbind(mass,binomial))[,1]),
    N_ID=length(unique(sample_num)), 
    ID=as.integer(as.factor(sample_num)),
    temp=as.numeric(unique(cbind(temp,sample_num))[,1]),
    mS=as.numeric(unique(cbind(ms,sample_num))[,1]),
    DO=as.numeric(unique(cbind(do_ch,sample_num))[,1]),
    pH=as.numeric(unique(cbind(ph,sample_num))[,1]),
    S=length(unique(site)), 
    site=as.integer(as.factor(unique(cbind(site,sample_num))[,1])), 
    n_waterhole=n_waterhole,
    last_visit=last_visit,
    detected=detected))
str(stan.data)

## MODEL 
model <- stan_model(file="edna_tax_enviro.stan")
fit <- sampling(model, data=stan.data, iter=5000, warmup=2500, chains=4, cores=2, thin=1,
                control = list(adapt_delta = 0.90, max_treedepth = 12))
# saveRDS(fit, "edna_tax_enviro_5k_ppcExact.rds")
fit <- readRDS("edna_tax_enviro_5k_ppcExact.rds")
# print(fit)


## Model Diagnostics
# Extracting posterior
posterior <- as.array(fit)

params <- c(names(fit)[grep("beta", names(fit))],
            names(fit)[grep("sigma", names(fit))])

outcome <- summary(fit, pars=params)

to_print <- outcome$summary[,c("mean","se_mean","sd","2.5%","97.5%","n_eff","Rhat")]

# kable(to_print, digits=2)

pairs(fit, pars=params)

# *Parameter Estimates*

require(bayesplot)

mcmc_areas(
  posterior, 
  pars = params,
  prob = 0.95, # 95% intervals
  prob_outer = 1.0,
  point_est = "mean"
)

# mcmc_combo(posterior, pars = params, combo = c("areas", "trace"))

# *Posterior predictive plots*
# Extracting detected_pred
detected_pred <- as.matrix(fit, pars = "detected_pred_exact")
y <- stan.data$detected
y_rep <- detected_pred

p_fit <- as.matrix(fit, pars = "p_fit")

ppc_scatter_avg(y = y, yrep = y_rep)
ppc_hist(y, y_rep[1:8, ], binwidth = 0.2) # boring
ppc_dens_overlay(y, y_rep[1:300, ]) #nice
ppc_stat_2d(y, y_rep, stat = c("mean", "sd")) # boring
ppc_error_hist(y, y_rep[1:12, ], binwidth = 0.01) + xlim(-1, 1) # boring

bw <- 0.075
ppc_stat_grouped(y,y_rep[1:200,], group=data$binomial, binwidth=bw)
ppc_stat_grouped(y,y_rep[1:200,], group=data$family, binwidth=bw)
ppc_stat_grouped(y,y_rep[1:200,], group=data$order, binwidth=bw)
ppc_stat_grouped(y,y_rep[1:200,], group=data$site, binwidth=bw)
ppc_stat_grouped(y,y_rep[1:200,], group=data$sample_num, binwidth=bw)

p_fit_means <- apply(p_fit, 2, mean)
p_fit_95_CI <- apply(p_fit, 2, function(x) quantile(x, c(0.05, 0.975)))
p_fit_50_CI <- apply(p_fit, 2, function(x) quantile(x, c(0.25, 0.75)))

max(p_fit)
max(p_fit_95_CI)

plot(p_fit_means, data$detected)
plot(p_fit_means ~ data$last_visit)
plot(p_fit_means ~ data$mass)

# p_fit verus last_visit
plot(p_fit_means ~ data$last_visit, pch=20, col="red")#, xaxt='n')
for (i in 1:100) {
    points(p_fit[i,] ~ data$last_visit, pch=19, col=rgb(0, 0, 0, 0.05))
    
}
points(p_fit_means ~ data$last_visit, pch=20, col="red")
last_visit_mean <- mean(dat$last_visit)
last_visit_sd <- sd(dat$last_visit)

# scale_half_last_visit <- function(x,mean,sd){
#         (x-mean(x)) * 0.5/sd(x)  
# }
# at_lv <- scale_half_last_visit(sqrt(c(1,2,6,12,24,36,48,72,48*2,48*3)), last_visit_mean, last_visit_sd)

# axis(1, at=at_lv, labels=c(1,2,6,12,24,36,48,72,48*2,48*3))

blank_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), 
                     axis.title=element_text(size=14))

blank_theme_large <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=18), 
                     axis.title=element_text(size=18))

d <- ggplot(dat, aes(last_visit, p_fit_means)) + 
    blank_theme_large + geom_point(alpha = 1, col="black", cex=1.8) + 
    scale_x_continuous(expand = c(0, 0.01), trans = "sqrt", breaks = c(0,1,2,6,12,24,36,48,72,48*2,48*3)) +
    scale_y_continuous(expand = c(0, 0.006), breaks = c(0,0.2,0.4,0.6,0.8,1.0)) + 
    geom_errorbar(aes(ymin=p_fit_50_CI[1,], ymax=p_fit_50_CI[2,]), width=0, col="black", alpha=0.25)+
    xlab("\nTime since last visit (hours)") + ylab("Probability of detection\n")
d

# inv_output <- plogis(outcome$summary)

# # Interpretation

# # for Lastseen (On probability scale)
# plogis(outcome$summary[3,1])


# # On odds scale)
# exp(outcome$summary[3,1])

# p_loss_raw <- 1-exp(outcome$summary[3,1])
# # ~ 97% decrease in detection probability per 0.5 standard deviation of square root of last_seen

# last_visit_mean <- mean(dat$last_visit)

# last_visit_halfsd <- sd(dat$last_visit)/2


# p_loss2 <- 1- 
# (exp(outcome$summary[3,1]) / last_visit_halfsd )


# TRYING TO AD ABLINE>>>>

# # Normalizing continuous predictors to mean = 0, sd = 0.5
# scale_half <- function(x){
#     return( (x-mean(x)) * 0.5/sd(x))
# }

# last_visit_mean <- mean(dat$last_visit)
# last_visit_sd <- sd(dat$last_visit)

# dat$last_visit
# sd(data$last_visit)
# (data$last_visit^2* last_visit_sd*2) + last_visit_mean^2

# scale_half(dat$last_visit)
# scale_half(sqrt(dat$last_visit))
# data$last_visit

# last_visit_mean_sqrt <- mean(dat$last_visit^2)
# last_visit_sd_sqrt <- sd(dat$last_visit^2)

# sd(dat$last_visit^2)
# sd(data$last_visit * last_visit_sd_sqrt*2)

# sd(data$last_visit * last_visit_sd_sqrt*2 + last_visit_mean_sqrt)
# dat$last_visit

# # 

# (data$last_visit* last_visit_sd*2) + last_visit_mean


# d + geom_abline(aes(slope=inv_output[3,1]^2,intercept=inv_output[1,1],color="blue"))

# This below can't be right? right?
# d + geom_smooth(method = "glm", 
#     method.args = list(family = "binomial"), 
#     se = F)


ggsave("p_detect_vs_last_visit_large.pdf", plot = d, device = "pdf",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)

ggsave("p_detect_vs_last_visit_large.png", plot = d, device = "png",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)


# p_fit verus n_waterhole
plot(p_fit_means ~ log(dat$n_waterhole), pch=20, col="red")
for (i in 1:100) {
    points(p_fit[i,] ~ log(dat$n_waterhole), pch=19, col=rgb(0, 0, 0, 0.05))
    
}
points(p_fit_means ~ log(dat$n_waterhole), pch=20, col="red")

dim(p_fit)


# pdetect by species
names(fit)
betas_species <- summary(fit, pars="species_eff")$summary
rownames(betas_species) <- unique(data$binomial)
betas_species


p_by_species <- data.frame(species=dat$binomial, 
                            p_fit_mean = p_fit_means, 
                            CI_50_low = p_fit_50_CI[1,],
                            CI_50_hi = p_fit_50_CI[2,])%>%
                    group_by(species)%>%
                    summarize(
                        mean_p = mean(p_fit_mean),
                        CI_50_low = mean(CI_50_low),
                        CI_50_hi = mean(CI_50_hi))%>%
                    arrange(-mean_p)

with(p_by_species, plot(mean_p))

species_labs <- gsub("_"," ", p_by_species$species)

d2 <- ggplot(p_by_species, aes(1:length(species), mean_p)) + 
    blank_theme_large + geom_point(alpha = 1, col="black", cex=1.8) + 
    scale_x_continuous(expand= c(0,1), breaks = 1:length(species_labs), labels=species_labs) +
    scale_y_continuous(limits=c(0,1), expand = c(0, 0.006), breaks = c(0,0.2,0.4,0.6,0.8,1.0)) + 
    geom_errorbar(aes(ymin=CI_50_low, ymax=CI_50_hi), width=0, col="black", alpha=0.25)+
    xlab("") + ylab("Probability of detection\n") + 
    theme(axis.text.x = element_text(angle = 65, hjust = 1))
d2

ggsave("p_detect_by_species_large.pdf", plot = d2, device = "pdf",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)


ggsave("p_detect_by_species_large.png", plot = d2, device = "png",
  scale = 1, width = 10, height = 8, units = "in",
  dpi = 300, limitsize = TRUE)

# try coloring by order? or last_seen? or n_waterhole?




# Last visit by species...

last_visit_by_species <- data.frame(species=dat$binomial,
                            last_visit=dat$last_visit, 
                            p_fit_mean = p_fit_means)%>%
                            group_by(species)%>%
                            summarize(
                            mean_p = mean(p_fit_mean),
                            mean_last_visit = mean(last_visit),
                            min_last_visit = min(last_visit),
                            max_last_visit = max(last_visit))%>%
                            arrange(mean_last_visit)



species_labs <- gsub("_"," ", last_visit_by_species$species)

d3 <- ggplot(last_visit_by_species, aes(1:length(species), mean_last_visit)) + 
    blank_theme + geom_point(alpha = 1, col="black", cex=1) + 
    scale_x_continuous(expand= c(0,1), breaks = 1:length(species_labs), labels=species_labs) +
    scale_y_continuous(limits=c(0,95), expand = c(0, 0.006), breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
    geom_errorbar(aes(ymin=min_last_visit, ymax=max_last_visit), width=0, col="black", alpha=0.25)+
    xlab("") + ylab("Mean time since last visit (hours)\n") + 
    theme(axis.text.x = element_text(angle = 65, hjust = 1))
d3


# estimating number of elephants in last x hours


# str(trap)

# ele1 <- trap %>% 
#         group_by(sample_num)%>%
#         mutate(lastphoto = max(timestamp))%>%
#         filter(timestamp>=(lastphoto-36*60*60))%>%
#         summarize(elephants_36h=sum(n_waterhole[binomial=="Loxodonta_africana"]))
# ele1



