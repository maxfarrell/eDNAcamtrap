# eDNA detection vs camtrap
# stan model with pooled exponential decay rate
# predicting detected of species within sites
# only estimating tau

require(rstan)

# sim data
source("sim_eDNA_data_updated.R")
# detected_site <- unique(subset(sim_data, select=c(detected,site)))[,1]
# stan.data <- with(sim_data, list(N=nrow(sim_data), S=length(unique(site)), 
# 	site=as.integer(site), time=time, n_waterhole=n_waterhole, detected=detected_site))

str(stan.data)

## STAN MODEL
# for dave (multi-core support):
rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())

model <- stan_model(file="camtrap_eDNA_host_site_pooledTau.stan")
fit <- sampling(model, data=stan.data, iter=4000, chains=2, cores=2, thin=1)
# save.image("camtrap_eDNA_host_detected_merged_noK_noAlpha.RData")
# load("camtrap_eDNA_host_detected_merged_noK_noAlpha.RData")
# traceplot(fit)
print(fit)



