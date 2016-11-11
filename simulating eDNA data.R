# simulating eDNA and camera data

rm(list=ls())

# functions
logit <- function(p){ 
	return(log(p/(1 - p)))
}

inv_logit <- function(x){
	return(exp(x)/(1+exp(x)))
}

# setup
n_sites <- 8
n_spec <- 1
n_minutes <- 7*24*60 # 7days in minutes

dat <- data.frame(spec=1:n_spec, time=rep(n_minutes:1,n_sites), site=rep(1:n_sites, each=n_minutes))
# dim(dat)#80640

# zero-inflated poisson for waterhole visitation
dat$n_waterhole <- ifelse(rbinom(dat$time, size = 1, prob = 0.9) > 0, 0, rpois(dat$time, lambda = 5))

# hist(dat$n_waterhole)
# hist(dat$n_waterhole[dat$n_waterhole>0])
# mean(dat$n_waterhole[dat$n_waterhole>0])

# aggregate(n_waterhole ~ site, dat, sum)

# creating variablitiy across sites
# thin two sites, so every even entry is removed
dat$n_waterhole[dat$site%in%c(1:2) & c(FALSE,TRUE)] <- 0

# for two sites remove visitation for second half of timeperiod 
dat$n_waterhole[dat$site%in%c(3:4) & dat$time<(max(dat$time)/2)] <- 0

# for two sites remove visitation for last two thirds of timeperiod 
dat$n_waterhole[dat$site%in%c(5:6) & dat$time<2*(max(dat$time)/3)] <- 0

# for one site remove visitation for last three quarters of timeperiod 
dat$n_waterhole[dat$site==7 & dat$time<3*(max(dat$time)/4)] <- 0

# aggregate(n_waterhole ~ site, dat, sum)

# Exponential weights
# w <- k * exp(-timeElapsed/tau)
# where w are the weights, k is scaling constant, tau is time constant for decay

k <-0.1
tau <- 24*60*1 # mean lifetime in minutes
dat$w <- k*exp(-dat$time/tau)
dat$edna <- dat$n_waterhole*dat$w
# hist(dat$edna)
# hist(dat$edna[dat$edna>0])
# mean(dat$edna[dat$edna>0])

# with(dat, plot(edna ~ time))
# with(dat, plot(w ~ time))
# with(dat, plot(edna ~ w))


# p detection is based on total 
# log(sum(dat$edna))

# inverse logit
# inv_logit(log10(aggregate(edna ~ site, dat, sum)$edna))

# Find some way to get detected (1/0)
edna_sums <- aggregate(edna ~ site, dat, sum)
names(edna_sums)[2] <- "total_edna"
# inv_logit(log10(edna_sums$total_edna))
edna_sums$p <- inv_logit(log10(edna_sums$total_edna))
edna_sums$detected <- rbinom(nrow(edna_sums), 1, edna_sums$p)
edna_sums$detected

# Looks okay after making some sites have low eDNA

# Merge them so that detected is repeated for each time observation
require(plyr)
dat <- join(dat,edna_sums)
str(dat)
stan.data <- with(dat, list(N=nrow(dat), time=time, n_waterhole=n_waterhole, 
	site=site, S=length(unique(site)), detected=detected))

# Don't merge them, just set the variables as list to read in
# stan.data <- with(dat, list(N=nrow(dat), time=time, n_waterhole=n_waterhole, 
# 	site_long=site, S=length(unique(site)), site_name=as.integer(unique(site)), detected=edna_sums$detected))

rm(list=ls()[grep("stan.data",ls(), invert=TRUE)])

# str(stan.data)

# SCRAPS

# threshold <- 1
# dat$detected <- 0
# dat$detected[dat$edna>=threshold] <- 1
# sum(dat$edna)

# x <- sum(dat$detected)
# x <- sum(dat$edna)
# x <- log(x)
# exp(x)/(1+exp(x))


# with(dat, plot(edna ~ time))

# hist(dat$detected)

