// Animal detection via waterhole eDNA
// MJ Farrell

data {
    int<lower=0> N;                 // # observations
    
    // Factors
    //int<lower=0> H;                 // number of host species
    //int<lower=1> host[N];           // map observations to hosts (group)
    
    int<lower=0> S;                 // number of sites
    int<lower=1> site[S];           // map observations to sites (group)
    
    // Reponse
    int<lower=0> detected[S];       // detected in eDNA (binary 0-1)
    
    // Time variable
    vector[N] time;                 // time coded as integer (minute) before sampling

    // Waterhole Visitation
    vector[N] n_waterhole;          // number of individuals seen 
                                    // at the waterhole at a given time

}

parameters {

    // parameters for exponential weighting of n_waterhole by time
    real<lower=0> k;                // scale parameter
    // maybe put upper bound on k to help limit nan return...?

    real<lower=0.000001> tau;       // mean lifetime (in units of time)

    // to modify
    real alpha;                     // global intercept
    real beta_n_waterhole;          // slope for sesPD
    
    // Random intercepts for site
    real int_site[S];
    // SD of random intercepts
    real<lower=0.000001> sigma_site;

    // Random intercepts for host
    // real int_host[H];
    // SD of random intercepts
    //real<lower=0.000001> sigma_host;


}

transformed parameters {

    vector[N] weights;
    vector[N] eDNA_N;
    vector[S] eDNA;
    real<lower=0, upper=1> p[S];

    weights = k * exp(time/log10(tau)); // log10 tau and sampler doesn't give nan warning... but not actually running...
    
    eDNA_N = weights .* n_waterhole;    
    
    for (s in 1:S) {
        eDNA[s] = sum(eDNA_N[site]);    
    }

    for (i in 1:S) {
        p[i] = inv_logit(alpha + beta_n_waterhole * eDNA[site[i]] + int_site[site[i]]);
    }

}

model {
    alpha ~ normal(0.0,1.0E3);      // 
    beta_n_waterhole ~ normal(0.0,1.0E3);      // 
    int_site ~ normal(0,sigma_site);//

    for (i in 1:S) {
        detected[i] ~ bernoulli(p[i]);
    }
    
}

