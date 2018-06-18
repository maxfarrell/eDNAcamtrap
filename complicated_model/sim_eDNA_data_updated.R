# simulating eDNA and camera data
# Updated using J Mihaljevic's framework
# https://jrmihalj.github.io/hierarchical-data-management-and-visualization/

rm(list=ls())

# source("merging_eDNA_camptrap_data.R")
# Doesn't work

set.seed(666)

# functions
logit <- function(p){ 
	return(log(p/(1 - p)))
}

inv_logit <- function(x){
	return(exp(x)/(1+exp(x)))
}


## Sample sizes:
sites_N <- 10 # Number of sites sampled
minutes_N <- 7*24*60 # 7 days in minutes

obs_N_perSite <- minutes_N
obs_N <- sum(obs_N_perSite)

species_N <- 20 # Number of species detected

dat <- data.frame(time=rep(minutes_N:1,sites_N), site=factor(rep(1:sites_N)), each=minutes_N)
# str(dat)
# head(dat)

# Simulating photographs for each species
dat$species <- factor(sample(c(1:species_N), obs_N, replace=TRUE, prob=round(rexp(species_N, 0.01))))

# Simulating number of individuals per species - zero inflated 
n_waterhole <- ifelse(rbinom(dat$time, size = 1, prob = 0.9) > 0, 0, rpois(dat$time, lambda = 10))
# # hist(n_waterhole)
# # sum(n_waterhole==0)/length(n_waterhole) # ~90 % zeros

dat$n_waterhole <- n_waterhole


# require(tidyr)
# require(ggplot2)
# ggplot(dat, aes(x=n_waterhole))+
#   geom_bar()+
#   facet_wrap(~site, ncol=5, scales = "free")

# ggplot(dat, aes(x=n_waterhole))+
#   geom_bar()+
#   facet_wrap(~species, ncol=5, scales = "free")


# Exponential weights
# w <- k * exp(-timeElapsed/tau)
# where w are the weights, k is scaling constant, tau is time constant for decay
k <-1
tau <- 24*60*1 # mean lifetime in minutes (1440)
dat$w <- k*exp(-dat$time/tau)
dat$edna <- dat$n_waterhole*dat$w
 
# Find some way to get detected (1/0)
require(dplyr)

edna_sums <- dat %>% 
	group_by(site, species) %>% 
	summarize(total_edna = sum(edna))

# inv_logit((edna_sums$total_edna))
edna_sums$p <- inv_logit((edna_sums$total_edna/100))
edna_sums$detected <- rbinom(nrow(edna_sums), 1, edna_sums$p)

# Some species are not found at some sites...
# sum(with(edna_sums, table(site , species)))

# sum(edna_sums$detected)/length(edna_sums$detected)# %60-70 site-species combinations detected.

# First attempt: Don't merge them, just set the variables as list to read in
# structure is a dataframe with photo data, and a dataframe with edna data

stan.data <- with(dat, list(N_photo=nrow(dat), SI=length(unique(site)), SP=length(unique(species)), 
	time=time, n_waterhole=n_waterhole, site_photo=site, species_photo=species,
	N_edna = length(edna_sums$detected), edna_detected=edna_sums$detected, 
	site_edna=edna_sums$site, species_edna=edna_sums$species))

# str(stan.data)

with(dat, plot(time,edna))
with(edna_sums, plot(time,edna))


rm(list=ls()[grep("stan.data",ls(), invert=TRUE)])

