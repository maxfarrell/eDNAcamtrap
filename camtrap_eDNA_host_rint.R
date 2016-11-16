# eDNA detection vs camtrap
# stan model with n_waterhole and exponential decay rate + site as random intercept

require(rstan)

# sim data
source("simulating eDNA data.R")
# detected_site <- unique(subset(sim_data, select=c(detected,site)))[,1]
# stan.data <- with(sim_data, list(N=nrow(sim_data), S=length(unique(site)), 
# 	site=as.integer(site), time=time, n_waterhole=n_waterhole, detected=detected_site))

str(stan.data)

## STAN MODEL
# for dave (multi-core support):
rstan_options(auto_write=TRUE)
options(mc.cores = parallel;;detectCores())
# model <- stan_model(file="camtrap_eDNA_detected_merged.stan")
# fit <- sampling(model, data=stan.data, iter=4000, chains=2, cores=2, thin=1)
# # save.image("camtrap_eDNA_host_detected_merged_alpha_noK.RData")

# # load("camtrap_eDNA_host_detected_merged.RData")
# traceplot(fit)
# print(fit)

# with k, and tau - no alpha...
# k est at 97.56
# tau at 955.82
# alpha not estimated...



# maybe try without k or parameters?
# updated simulation to have tau at 0.5 of a day (720 mins)
# and alpha = 0 (obviously wouldn't expect any chance of detection when eDNA is zero) 
# and k=1 (no effect of scale parameter in simulation)
model <- stan_model(file="camtrap_eDNA_detected_merged_noK_noAlpha.stan")
fit <- sampling(model, data=stan.data, iter=4000, chains=2, cores=2, thin=1)
save.image("camtrap_eDNA_host_detected_merged_noK_noAlpha.RData")
# load("camtrap_eDNA_host_detected_merged_noK_noAlpha.RData")
# traceplot(fit)
print(fit)







# ########## STOPPED HERE ##############


# # load("virulence_sesPD_host_simData.RData")
# # alpha = -11.65
# # beta_sesPD = -2.01
# # Looks good!

# stan_diag(fit, information = c("sample","stepsize", "treedepth","divergence"),
# 			chain = 0)
# plot(fit, pars=c("alpha","beta_sesPD","int_host"))


# # NEED TO TRANSFORM TIME TO MINUTES / SECONDS BEFORE SAMPLING>>>>

# # real data
# data <- read.csv("../Final/Clean Data/Virulence_merged_data_oct6_2016.csv")
# names(data) <- tolower(names(data))
# names(data) <- gsub("\\.", "_", names(data))
# data$multi_strain[is.na(data$multi_strain)] <- 0

# # There are NAs in sesPD metrics for single host parasites... must remove them before fitting in stan
# data <- data[!is.na(data$ses_pd_z_euth),]

# # There are some entries that have cases==0
# # REMOVE IN CLEANING SCRIPT
# # doing here for now..
# data <- data[data$cases>0,]

# # real data
# stan.data <- with(data, list(N=nrow(data), deaths=deaths, 
# 	cases=cases, sesPD=ses_pd_z_euth, H=length(unique(host)), host=as.integer(host)))

# ## STAN MODEL
# model <- stan_model(file="virulence_sesPD_host.stan")
# fit <- sampling(model, data=stan.data, iter=4000, chains=4, cores=2, thin=1)
# # default iter is 2000, warmup (burnin) is default half of iter
# # didn't completely converge after 2000. Rhat ranges 1.01-1.08
# # save.image("virulence_sesPD_host_realData.RData")

# # load("virulence_sesPD_host_realData.RData")
# traceplot(fit)
# # looking at quantiles, n_eff, and Rhat
# print(fit)

# # extracting parameters
# beta_sesPD <- extract(fit, 'beta_sesPD')
# beta_sesPD <- unlist(beta_sesPD)
# hist(beta_sesPD) 


# print(fit, pars=c("alpha","beta_sesPD","int_host"))

# p_int <- exp(-2.58)/(1 + exp(-2.58)) # inverse logit
# p_slope <- exp(-0.33)/(1 + exp(-0.33)) # inverse logit
# p_int 
# p_slope

# str(data)

# pred_p <- exp(-2.58 + -0.33*data$ses_pd_z_euth)/(1+exp(-2.58 + -0.33*data$ses_pd_z_euth))
# mean(pred_p)

# # ACTUAL PLOT: mortality rate vs sesPD
# # plot(data$ses_pd_z_euth, data$deaths/data$cases)
# plot(data$ses_pd_z_euth, pred_p, ylim=c(0,1), col="red")

# pred_deaths <- rbinom(data$cases,data$cases,pred_p)
# hist(pred_deaths)
# hist(data$deaths/data$cases)
# hist(pred_deaths/data$cases)
# # still funky - actual death rate is more likely to be 0 or 1
# # predicted has more even distribution 