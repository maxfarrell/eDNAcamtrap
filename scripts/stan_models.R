# COI pipeline
# Stan models for detection

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

# loading data
load("../data/merged_eDNA_camtrap_data_nov20_2019.RData")

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

# including time of sample collection to calclate time since last visit for each species per sample
collected <- seqs %>% select(sample_num, collected) %>% unique()

dat1 <- trap %>% full_join(collected)%>%
        filter(n_waterhole!=0)%>%
        group_by(sample_num, binomial)%>%
        mutate(n_waterhole=sum(n_waterhole), 
            seen=((n_waterhole>0)*1),
            last_visit=as.numeric(unique(collected)-max(timestamp), unit="hours")) %>%
        # select(site,sample_num,AB,binomial,n_waterhole,last_visit,seen)%>%
        select(site,sample_num,binomial,n_waterhole,last_visit,seen)%>%
        unique()%>%
        arrange(site, sample_num, binomial)     

names(seqs) <- tolower(names(seqs))
seqs <- subset(seqs, select=-c(sample))

dat2 <- seqs %>%
        filter(is.na(s_xs))%>%
        # group_by(sample_num, AB, species)%>%
        group_by(sample_num, species)%>%
        mutate(detected=1, abundance=sum(abundance))%>%
        select(site, sample_num, species, detected, abundance) %>%
        unique() %>%
        arrange(site, sample_num, species)      

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
        arrange(site, sample_num, binomial)

dat <- left_join(dat, dat3)

dat$detected[is.na(dat$detected)] <- 0

sort(unique(dat$binomial))

dat[dat$binomial=="unknown",]
# remove these species
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
# Vulture and Guinea fowl

# View(dat)
# There are two birds, and two mammals with positive detections...


####################
##  brms analyses ##
####################
require(brms)
require(kableExtra)
require(xtable)

# model with site, species, family, order hierarchical effects

# Removing rows where animals were not seen
dat <- dat[!is.na(dat$seen),]

# Prior predictive distribution
m1 <- brm(detected ~ log(n_waterhole) + sqrt(last_visit) + log(mass) + 
                            temp + ms + do_ch + ph + 
                            (1|site) + 
                            (1|binomial),
    # custom priors based on McElreath
    prior = c(prior(normal(0, 1), class = Intercept),
                prior(normal(0, 1.5), class = b),
                prior(normal(0, 1), class = sd)),
  data=dat, family=bernoulli(), 
  iter=1000,
  control=list(adapt_delta=0.95, max_treedepth=10),
  sample_prior="only")

pp <- brms::posterior_predict(m1)
ppc_dens_overlay(y = dat$detected, yrep=pp[1:50,])


# Fitting model


if (!file.exists("../results/m1.rds")) {

m1 <- brm(detected ~ log(n_waterhole) + sqrt(last_visit) + log(mass) + 
                            temp + ms + do_ch + ph + 
                            (1|site) + 
                            (1|binomial),
    # custom priors based on McElreath
    prior = c(prior(normal(0, 1), class = Intercept),
                prior(normal(0, 1.5), class = b),
                prior(normal(0, 1), class = sd)),
  data=dat, family=bernoulli(), 
  iter=10000,
  control=list(adapt_delta=0.95, max_treedepth=10))

  saveRDS(m1, "../results/m1.rds")

} else { m1 <- readRDS("../results/m1.rds")}

summary(m1)

outcome <- summary(m1$fit)
to_print <- outcome$summary[grep("^[sd,b]_*",rownames(outcome$summary)),]
to_print <- to_print[-1,]
names1 <- c("log (N waterhole)","sqrt (last visit)","log (mass)",
            "Temperature","mS","DO","pH","sigma species","sigma site")

## nicer mcmc plot
color_scheme_set("darkgray")
fp_m1 <- mcmc_plot(m1, pars=rownames(to_print), type="intervals", 
                    prob=0.80, prob_outer=0.95, point_est="mean") + scale_y_discrete(labels=rev(names1), limits=rev(rownames(to_print)))

fp_m1 <- fp_m1 + theme_classic() 

ggsave("../plots_tables/fp_m1.pdf",fp_m1, width=5, height=4)

# fp_m1F

## Betas plot
output <- summary(m1$fit)$summary
betas <- output[grep("b_", rownames(output)),]
sigmas <- output[grep("sd_", rownames(output)),]
r_biom <- output[grep("r_binom", rownames(output)),]
r_site <- output[grep("r_site", rownames(output)),]

# plot effects
col4fig <- c("mean","sd","2.5%","25%","50%","75%","97.5%","n_eff","Rhat")
beta_params <- rownames(betas)
betas.table <- output[beta_params,col4fig]
betas.table <- betas.table[-1,]

rownames(betas.table) <- c("log (N waterhole)",
                          "sqrt (time since last visit)",
                          "log (mass)",
                          "Temperature  (C)",
                          "Conductivity (mS)",
                          "Dissolved  Oxygen",
                          "pH")


pdf("../plots_tables/forestPlot_m1.pdf", height=5, width=7)
par(mar = c(2.2, 11.5, 0.5, 0.5))
plot(seq(-1.8, 2.5, length.out = nrow(betas.table)), 
       1:nrow(betas.table),
       type="n",
       xlab = "",
       ylab = "",
       yaxt = "n",
       xaxt="n")
  
axis(2, at = nrow(betas.table):1, labels = rownames(betas.table), las = 1, cex.axis = 1.1)
axis(1, cex.axis = 1.1)

arrows(betas.table[,"2.5%"], nrow(betas.table):1, betas.table[,"97.5%"], nrow(betas.table):1,
       len = 0, col = "dodgerblue", lwd=4, lty=1)
arrows(betas.table[,"75%"], nrow(betas.table):1, betas.table[,"25%"], nrow(betas.table):1,
       len = 0, col = "#1564b2", lwd=8)
points(betas.table[,'mean'],
       nrow(betas.table):1,
       pch = 19,
       col = "#0c3966", cex=1.5, fill="white")
abline(v = 0, lty = 2, col = "darkgray")
dev.off()

# pdf("../plots_tables/fp_m1.pdf", width=5, height=4)
# fp_m1 + ggtitle("")
# dev.off()



# Table
rownames(to_print) <- names1
to_print <- to_print[,c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")]
print(xtable(to_print, type = "latex",caption="Esimtated model parameters for model 1",
                        digits=2), file = "../plots_tables/m1.tex")


##Posterior predictive checks
pp <- brms::posterior_predict(m1)

# pdf("../plots_tables/pp_m1.pdf", width=5, height=4)
# ppc_dens_overlay(y = dat$detected, yrep=pp[1:50,])
# dev.off()

prop_zero <- function(x) mean(x == 0)
prop_zero(dat$detected) # check proportion of zeros in y
# 0.76
ppc_stat(dat$detected, pp[1:2000,], stat = "prop_zero", binwidth = 0.005)

prop_one <- function(x) mean(x == 1)
prop_one(dat$detected) # check proportion of zeros in y
# 0.24
ppc_stat(dat$detected, pp[1:2000,], stat = "prop_one", binwidth = 0.005)

# these look okay

bayes_R2(m1)#0.29



###############################
## Summary plots for model 1 ##
###############################


# marginal effects for last visit
lv_plot <- conditional_effects(m1, "last_visit", method="fitted")
lv1 <- plot(lv_plot, theme=theme_classic(), line_args=list(lty=2,col="#0c3966"))$last_visit

lv1 <- lv1 + labs(y="Probability of detection",x="Time since last visit (hours)") + 
      scale_x_continuous(expand= c(0,1),breaks=seq(0,24*7, by=12)) + 
      scale_y_continuous(expand= c(0.02,0.02),breaks=seq(0,1, by=0.1)) +
      geom_ribbon(fill="dodgerblue", alpha=0.2)

lv1
ggsave("../plots_tables/last_visit_marginal.pdf",lv1, height=4.5, width=6)


### Posterior predictions by species

########
# pp_check(m1, nsamples = 100, type = "stat_grouped", group = "binomial")

pp <- brms::posterior_predict(m1)

d <- data.frame(species=factor(dat$binomial), yrep=t(pp))
colnames(d) <- gsub(".", "_", colnames(d), fixed = TRUE)
molten_d <- reshape2::melt(d, id.vars = "species")
molten_d %<>% group_by(species) %>% 
              summarise(mean=mean(value), 
                var=var(value),
                low=mean-var, high=mean+var) %>%
              arrange(-mean)
molten_d

species_labs <- gsub("_"," ", molten_d$species)

blank_theme_large <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=18), 
                     axis.title=element_text(size=18))


d2 <- ggplot(molten_d, aes(1:length(species), mean)) + 
    blank_theme_large +  
    scale_x_continuous(expand= c(0,1), breaks = 1:length(species_labs), labels=species_labs) +
    scale_y_continuous(limits=c(0,0.8), expand = c(0, 0.02), breaks = c(0,0.2,0.4,0.6,0.8)) + 
    geom_errorbar(aes(ymin=low, ymax=high), width=0, col="dodgerblue", alpha=1, lty=2) +
    geom_point(alpha = 1, col="#0c3966", cex=3) +
    xlab("") + ylab("Probability of detection\n") + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
d2

ggsave("../plots_tables/pp_by_species.pdf",height=9, width=11)


# Plot of estimated species effects
m1
outcome <- as.data.frame(summary(m1$fit)$summary)
str(outcome)

sp_effs <- outcome[grep("r_binomial",rownames(outcome)),]
names(sp_effs) <- c("mean","se_mean","sd","CI_2_5","CI_25","CI_50","CI_75","CI_97_5","n_eff","Rhat")
sp_effs$species <- sort(unique(dat$binomial))

sp_effs <- arrange(sp_effs, -mean)
species_labs <- gsub("_"," ", sp_effs$species)


d3 <- ggplot(sp_effs, aes(1:length(species), mean)) + 
    blank_theme_large +  
    scale_x_continuous(expand= c(0,1), breaks = 1:length(species_labs), labels=species_labs) +
    # scale_y_continuous(limits=c(0,0.8), expand = c(0, 0.02), breaks = c(0,0.2,0.4,0.6,0.8)) + 
    geom_errorbar(aes(ymin=CI_25, ymax=CI_75), width=0, col="dodgerblue", alpha=1, lty=2) +
    geom_point(alpha = 1, col="#0c3966", cex=3) +
    xlab("") + ylab("Species Effects\n") + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
d3

ggsave("../plots_tables/species_effects_50CI.pdf",height=9, width=11)



########################
## ALTERNATIVE MODELS ##
########################


# Same model with hierarchical effect per sample

if (!file.exists("../results/m1_1.rds")) {

m1.1 <- brm(detected ~ log(n_waterhole) + sqrt(last_visit) + log(mass) + 
                            temp + ms + do_ch + ph + 
                            (1|site) + 
                            (1|binomial) +
                            (1|sample_num),
    # custom priors based on McElreath
    prior = c(prior(normal(0, 1), class = Intercept),
                prior(normal(0, 1.5), class = b),
                prior(normal(0, 1), class = sd)),
  data=dat, family=bernoulli(), 
  iter=10000,
  control=list(adapt_delta=0.95, max_treedepth=10))

  saveRDS(m1.1, "../results/m1_1.rds")

} else { m1.1 <- readRDS("../results/m1_1.rds")}

summary(m1.1)

outcome <- summary(m1.1$fit)
to_print <- outcome$summary[grep("^[sd,b]_*",rownames(outcome$summary)),]
to_print <- to_print[-1,]
names1 <- c("log (N waterhole)","sqrt (last visit)","log (mass)",
            "Temperature","mS","DO","pH","sigma species","sigma sample","sigma site")

fp_m1.1 <- mcmc_plot(m1.1, pars=rownames(to_print), type="intervals", 
                    prob=0.8, prob_outer=0.95, point_est="mean") + scale_y_discrete(labels=rev(names1), limits=rev(rownames(to_print)))

fp_m1.1

# pdf("../plots_tables/fp_m1.1.pdf", width=5, height=4)
# fp_m1.1 + ggtitle("")
# dev.off()

# Table
rownames(to_print) <- names1
to_print <- to_print[,c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")]
print(xtable(to_print, type = "latex",caption="Esimtated model parameters for model 1",
                        digits=2), file = "../plots_tables/m1.1.tex")


##Posterior predictive checks
pp <- brms::posterior_predict(m1.1)
ppc_dens_overlay(y = dat$detected, yrep=pp[1:50,])

# pdf("../plots_tables/pp_m1.1.pdf", width=5, height=4)
# ppc_dens_overlay(y = dat$detected, yrep=pp[1:50,])
# dev.off()

prop_zero <- function(x) mean(x == 0)
prop_zero(dat$detected) # check proportion of zeros in y
# 0.76
ppc_stat(dat$detected, pp[1:2000,], stat = "prop_zero", binwidth = 0.005)

prop_one <- function(x) mean(x == 1)
prop_one(dat$detected) # check proportion of zeros in y
# 0.24
ppc_stat(dat$detected, pp[1:2000,], stat = "prop_one", binwidth = 0.005)

# these look okay

bayes_R2(m1.1)#0.30



#####################################
## Adding Family and Order effects ##
#####################################


# Prior predictive distribution
m2 <- brm(detected ~ log(n_waterhole) + sqrt(last_visit) + log(mass) + 
                            temp + ms + do_ch + ph + 
                            (1|site) + 
                            (1|binomial) +
                            (1|family) +
                            (1|order),
    # custom priors based on McElreath
    prior = c(prior(normal(0, 1), class = Intercept),
                prior(normal(0, 1.5), class = b),
                prior(normal(0, 1), class = sd)),
  data=dat, family=bernoulli(), 
  iter=1000,
  control=list(adapt_delta=0.95, max_treedepth=10),
  sample_prior="only")

pp <- posterior_predict(m2, draws=50)
ppc_dens_overlay(y = dat$detected, yrep=pp)


# Fitting model
if (!file.exists("../results/m2.rds")) {

  m2 <- brm(detected ~ log(n_waterhole) + sqrt(last_visit) + log(mass) + 
                              temp + ms + do_ch + ph + 
                              (1|site) + 
                              (1|binomial) +
                              (1|family) +
                              (1|order),
      # custom priors based on McElreath
      prior = c(prior(normal(0, 1), class = Intercept),
                  prior(normal(0, 1.5), class = b),
                  prior(normal(0, 1), class = sd)),
    data=dat, family=bernoulli(), 
    iter=10000,
    control=list(adapt_delta=0.95, max_treedepth=10))
  
    saveRDS(m2, "../results/m2.rds")

} else { m2 <- readRDS("../results/m2.rds")}

summary(m2)

## Model Diagnostics

outcome <- summary(m2$fit)
to_print <- outcome$summary[grep("^[sd,b]_*",rownames(outcome$summary)),]
to_print <- to_print[-1,]
names1 <- c("log (N waterhole)","sqrt (last visit)","log (mass)",
            "Temperature","mS","DO","pH","sigma species",
            "sigma_family", "sigma_order","sigma site",)

fp_m2 <- mcmc_plot(m2, pars=rownames(to_print), type="intervals", 
                    prob=0.8, prob_outer=0.95, point_est="mean") + scale_y_discrete(labels=rev(names1), limits=rev(rownames(to_print)))

fp_m2

pdf("../plots_tables/fp_m2.pdf", width=5, height=4)
fp_m2 + ggtitle("")
dev.off()

# Table
rownames(to_print) <- names1
to_print <- to_print[,c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")]
print(xtable(to_print, type = "latex",caption="Esimtated model parameters for model 2",
                        digits=2), file = "../plots_tables/m2.tex")


##Posterior predictive checks
pp <- brms::posterior_predict(m2)
ppc_dens_overlay(y = dat$detected, yrep=pp[1:50,])

# pdf("../plots_tables/pp_m2.pdf", width=5, height=4)
# ppc_dens_overlay(y = dat$detected, yrep=pp[1:500,])
# dev.off()


prop_zero <- function(x) mean(x == 0)
prop_zero(dat$detected) # check proportion of zeros in y
# 0.76
ppc_stat(dat$detected, pp[1:2000,], stat = "prop_zero", binwidth = 0.005)

prop_one <- function(x) mean(x == 1)
prop_one(dat$detected) # check proportion of zeros in y
# 0.24
ppc_stat(dat$detected, pp[1:2000,], stat = "prop_one", binwidth = 0.005)

# these look quite good though...

bayes_R2(m2)#0.3159


### Adding sample level hierarchical effect
# Fitting model

if (!file.exists("../results/m2_2.rds")) {

  m2.2 <- brm(detected ~ log(n_waterhole) + sqrt(last_visit) + log(mass) + 
                              temp + ms + do_ch + ph + 
                              (1|site) + 
                              (1|binomial) +
                              (1|family) +
                              (1|order) +
                              (1|sample_num),

      # custom priors based on McElreath
      prior = c(prior(normal(0, 1), class = Intercept),
                  prior(normal(0, 1.5), class = b),
                  prior(normal(0, 1), class = sd)),
    data=dat, family=bernoulli(), 
    iter=10000,
    control=list(adapt_delta=0.95, max_treedepth=10))
  
    saveRDS(m2.2, "../results/m2_2.rds")

} else { m2.2 <- readRDS("../results/m2_2.rds")}

summary(m2.2)

## Model Diagnostics

outcome <- summary(m2.2$fit)
to_print <- outcome$summary[grep("^[sd,b]_*",rownames(outcome$summary)),]
to_print <- to_print[-1,]
names1 <- c("log (N waterhole)","sqrt (last visit)","log (mass)",
            "Temperature","mS","DO","pH","sigma species",
            "sigma_family", "sigma_order","sigma sample","sigma site")

fp_m2.2 <- mcmc_plot(m2.2, pars=rownames(to_print), type="intervals", 
                    prob=0.8, prob_outer=0.95, point_est="mean") + scale_y_discrete(labels=rev(names1), limits=rev(rownames(to_print)))

fp_m2.2

# pdf("../plots_tables/fp_m2.2.pdf", width=5, height=4)
# fp_m2.2 + ggtitle("")
# dev.off()

# Table
rownames(to_print) <- names1
to_print <- to_print[,c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")]
print(xtable(to_print, type = "latex",caption="Esimtated model parameters for model 2",
                        digits=2), file = "../plots_tables/m2.2.tex")


##Posterior predictive checks
pp <- brms::posterior_predict(m2.2)
ppc_dens_overlay(y = dat$detected, yrep=pp[1:50,])

# pdf("../plots_tables/pp_m2.2.pdf", width=5, height=4)
# ppc_dens_overlay(y = dat$detected, yrep=pp[1:500,])
# dev.off()


prop_zero <- function(x) mean(x == 0)
prop_zero(dat$detected) # check proportion of zeros in y
# 0.76
ppc_stat(dat$detected, pp[1:2000,], stat = "prop_zero", binwidth = 0.005)

prop_one <- function(x) mean(x == 1)
prop_one(dat$detected) # check proportion of zeros in y
# 0.24
ppc_stat(dat$detected, pp[1:2000,], stat = "prop_one", binwidth = 0.005)

# these look quite good though...

bayes_R2(m2.2)#0.324





#############################

comp <- loo_compare(m1.loo, m1.1.loo, m2.loo, m2.2.loo)
print(comp, simplify=FALSE)

# looks like m1 is mildly favored



#### Dropping species which were never detected...

d_summ <- dat %>% group_by(binomial) %>% summarise(detections=sum(detected),seens=sum(seen), rate = detections/seens)
never_seen <- d_summ$binomial[d_summ$detections==0]

dat_detected <- dat[!dat$binomial%in%never_seen,] 

# Fitting model
if (!file.exists("../results/m2_det.rds")) {

  m2.det <- brm(detected ~ log(n_waterhole) + sqrt(last_visit) + log(mass) + 
                              temp + ms + do_ch + ph + 
                              (1|site) + 
                              (1|binomial) +
                              (1|family) +
                              (1|order)  +
                              (1|sample_num),
      # custom priors based on McElreath
      prior = c(prior(normal(0, 1), class = Intercept),
                  prior(normal(0, 1.5), class = b),
                  prior(normal(0, 1), class = sd)),
    data=dat_detected, family=bernoulli(), 
    iter=10000,
    control=list(adapt_delta=0.95, max_treedepth=10))
  
    saveRDS(m2.det, "../results/m2_det.rds")

} else { m2.det <- readRDS("../results/m2_det.rds")}

summary(m2.det)

## Model Diagnostics

outcome <- summary(m2.det$fit)
to_print <- outcome$summary[grep("^[sd,b]_*",rownames(outcome$summary)),]
to_print <- to_print[-1,]
names1 <- c("log (N waterhole)","sqrt (last visit)","log (mass)",
            "Temperature","mS","DO","pH","sigma species",
            "sigma family", "sigma order","sigma sample","sigma site")

fp_m2.det <- mcmc_plot(m2.det, pars=rownames(to_print), type="intervals", 
                    prob=0.8, prob_outer=0.95, point_est="mean") + scale_y_discrete(labels=rev(names1), limits=rev(rownames(to_print)))

fp_m2.det

# pdf("../plots_tables/fp_m2.det.pdf", width=5, height=4)
# fp_m2.det + ggtitle("")
# dev.off()

# Table
rownames(to_print) <- names1
to_print <- to_print[,c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")]
print(xtable(to_print, type = "latex",caption="Esimtated model parameters for model 2",
                        digits=2), file = "../plots_tables/m2.det.tex")


##Posterior predictive checks
pp <- brms::posterior_predict(m2.det)
ppc_dens_overlay(y = dat_detected$detected, yrep=pp[1:50,])

# pdf("../plots_tables/pp_m2.det.pdf", width=5, height=4)
# ppc_dens_overlay(y = dat_detected$detected, yrep=pp[1:500,])
# dev.off()


prop_zero <- function(x) mean(x == 0)
prop_zero(dat_detected$detected) # check proportion of zeros in y
# 0.67
ppc_stat(dat_detected$detected, pp[1:2000,], stat = "prop_zero", binwidth = 0.005)

prop_one <- function(x) mean(x == 1)
prop_one(dat_detected$detected) # check proportion of zeros in y
# 0.33
ppc_stat(dat_detected$detected, pp[1:2000,], stat = "prop_one", binwidth = 0.005)

# these look quite good though...

bayes_R2(m2.det)#0.296
