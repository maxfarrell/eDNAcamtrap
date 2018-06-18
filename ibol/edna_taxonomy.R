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

dat <- full_join(dat1,dat2) %>%
        # arrange(site, sample_num, binomial, AB)
        arrange(site, sample_num, binomial)
head(dat)

dat$detected[is.na(dat$detected)] <- 0
dat <- dat[dat$binomial!="unknown",]
dat <- dat[!is.na(dat$binomial),]
# one instance of seeing an animal but not detecting it -
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
    S=length(unique(site)), 
    site=as.integer(as.factor(unique(cbind(site,sample_num))[,1])), 
    n_waterhole=n_waterhole,
    last_visit=last_visit,
    detected=detected))
str(stan.data)

with(stan.data, order[family[species]])

## MODEL 
model <- stan_model(file="edna_taxonomy.stan")
fit <- sampling(model, data=stan.data, iter=5000, warmup=2500, chains=4, cores=2, thin=1,
                control = list(adapt_delta = 0.99))
saveRDS(fit, "edna_taxonomy_5k_ppcExact.rds")
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

ppc_stat_grouped(y,y_rep[1:200,], group=data$binomial)
ppc_stat_grouped(y,y_rep[1:200,], group=data$family)
ppc_stat_grouped(y,y_rep[1:200,], group=data$order)
ppc_stat_grouped(y,y_rep[1:200,], group=data$site)
ppc_stat_grouped(y,y_rep[1:200,], group=data$sample_num)

dim(y_rep)

dim(p_fit)
p_fit_means <- apply(p_fit, 2, mean)

plot(p_fit_means, data$detected)
plot(p_fit_means ~ data$last_visit)

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
    

d <- ggplot(dat, aes(last_visit, p_fit_means))
d + geom_point(alpha = 1, col="red") + 
    scale_x_continuous(trans = "sqrt", breaks = c(1,2,6,12,24,36,48,72,48*2,48*3)) +
for (i in 1:100) {
    points(p_fit[i,] ~ data$last_visit, pch=19, col=rgb(0, 0, 0, 0.05), ylim=c(-1,1))
}


# p_fit verus n_waterhole
plot(p_fit_means ~ log(dat$n_waterhole), pch=20, col="red")
for (i in 1:100) {
    points(p_fit[i,] ~ log(dat$n_waterhole), pch=19, col=rgb(0, 0, 0, 0.05))
    
}
points(p_fit_means ~ log(dat$n_waterhole), pch=20, col="red")


