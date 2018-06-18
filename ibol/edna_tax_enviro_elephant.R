# Animal detection via waterhole eDNA
# Basic model - pooling AB replicates
# Includes taxonomy (family + order)
# Includes environmental data (water quality per sample)
# Includes observations of elephants in 24hours prior to sampling
# MJ Farrell

# Last updated November 8th 2017

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


ele1 <- trap %>% 
        group_by(sample_num)%>%
        mutate(lastphoto = max(timestamp))%>%
        filter(timestamp>=(lastphoto-48*60*60))%>%
        summarize(elephants_36h=sum(n_waterhole[binomial=="Loxodonta_africana"]))
ele1


dat <- full_join(dat1,dat2) %>%
        # arrange(site, sample_num, binomial, AB)
        arrange(site, sample_num, binomial)

dat <- left_join(dat, dat3)

dat <- left_join(dat, ele1)

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

# Removing species what we do not have barcode sequences for
dat <- dat[!dat$binomial%in%c("Mungos_mungo","Paraxerus_cepapi"),]
sort(unique(dat$binomial))
sort(unique(dat$sample_num))

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
data$temp <- scale_half(data$temp)
# with(dat, hist((ms)))
data$ms <- scale_half(data$ms)
# with(dat, hist((do_ch)))
data$do_ch <- scale_half(data$do_ch)
# with(dat, hist((ph)))
data$ph <- scale_half(data$ph)
# names(dat)
# pairs(dat[,c(4:5,8:15)])
# with(dat, hist(log(elephants_36h),breaks=20))
data$elephants_36h <- scale_half(log(data$elephants_36h))

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
    ele36=as.numeric(unique(cbind(elephants_36h,sample_num))[,1]),
    S=length(unique(site)), 
    site=as.integer(as.factor(unique(cbind(site,sample_num))[,1])), 
    n_waterhole=n_waterhole,
    last_visit=last_visit,
    detected=detected))
str(stan.data)

## MODEL 
model <- stan_model(file="edna_tax_enviro_elephant.stan")
fit <- sampling(model, data=stan.data, iter=5000, warmup=2500, chains=4, cores=2, thin=1,
                control = list(adapt_delta = 0.99, max_treedepth = 12))
# saveRDS(fit, "edna_tax_enviro_elephant_5k_ppcExact.rds")
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

d <- ggplot(dat, aes(last_visit, p_fit_means)) + 
    blank_theme + geom_point(alpha = 1, col="black", cex=1) + 
    scale_x_continuous(expand = c(0, 0.01), trans = "sqrt", breaks = c(0,1,2,6,12,24,36,48,72,48*2,48*3)) +
    scale_y_continuous(expand = c(0, 0.006), breaks = c(0,0.2,0.4,0.6,0.8,1.0)) + 
    geom_errorbar(aes(ymin=p_fit_50_CI[1,], ymax=p_fit_50_CI[2,]), width=0, col="black", alpha=0.25)+
    xlab("\nTime since last visit (hours)") + ylab("Probability of detection\n")
d

ggsave("p_detect_vs_last_visit.pdf", plot = d, device = "pdf",
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
    blank_theme + geom_point(alpha = 1, col="black", cex=1) + 
    scale_x_continuous(expand= c(0,1), breaks = 1:length(species_labs), labels=species_labs) +
    scale_y_continuous(limits=c(0,1), expand = c(0, 0.006), breaks = c(0,0.2,0.4,0.6,0.8,1.0)) + 
    geom_errorbar(aes(ymin=CI_50_low, ymax=CI_50_hi), width=0, col="black", alpha=0.25)+
    xlab("") + ylab("Probability of detection\n") + 
    theme(axis.text.x = element_text(angle = 65, hjust = 1))
d2

ggsave("p_detect_by_species.pdf", plot = d2, device = "pdf",
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


str(trap)

ele1 <- trap %>% 
        group_by(sample_num)%>%
        mutate(lastphoto = max(timestamp))%>%
        filter(timestamp>=(lastphoto-36*60*60))%>%
        summarize(elephants_36h=sum(n_waterhole[binomial=="Loxodonta_africana"]))
ele1




# p_fit verus mass
plot(p_fit_means ~ log(dat$mass), pch=20, col="red")
for (i in 1:100) {
    points(p_fit[i,] ~ log(dat$mass), pch=19, col=rgb(0, 0, 0, 0.05))
    
}
points(p_fit_means ~ log(dat$mass), pch=20, col="red")


