// Animal detection via waterhole eDNA
// detected per site merged with data frame
// simple model - no site intercept
// MJ Farrell

data {
    int<lower=0> N;                 // # observations
    
    // Factors
    //int<lower=0> H;                 // number of host species
    //int<lower=1> host[N];           // map observations to hosts (group)
    
    # int<lower=0> S;                 // number of sites
    int<lower=1> site[N];           // map observations to sites (group)
    
    // Reponse
    int<lower=0, upper=1> detected[N];       // detected in eDNA (binary 0-1)
    
    // Time variable
    real time[N];                 // time coded as integer (minute) before sampling

    // Waterhole Visitation
    real n_waterhole[N];          // number of individuals seen 
                                    // at the waterhole at a given time

}

parameters {

    // parameters for exponential weighting of n_waterhole by time
    // real k;                // scale parameter
    // maybe put upper bound on k to help limit nan return...?

    real<lower=0.000001> tau;       // mean lifetime (in units of time)

    // to modify
    // real alpha;                     // global intercept
    // real beta_n_waterhole;          // slope for n_waterhole
    
    // Random intercepts for site
    // real int_site[S];
    // SD of random intercepts
    // real<lower=0.000001> sigma_site;

    // Random intercepts for host
    // real int_host[H];
    // SD of random intercepts
    //real<lower=0.000001> sigma_host;


}

transformed parameters {

    // vector[N] weights;
    // vector[N] eDNA_N;
    // vector[S] eDNA;
    // real<lower=0, upper=1> p[N];

    //weights = k * exp((time/10080)/tau); // time is divided by 10080 (max of time)
    //eDNA_N = weights .* n_waterhole;    
    
    //for (s in 1:S) {
    //    eDNA[s] = sum(eDNA_N[site]);    
    //}

   /* for (i in 1:N) {
        p[i] = inv_logit(alpha + k * exp((time[i]/10080)/tau) * n_waterhole[i]);
    }
*/

}

model {
    // alpha ~ normal(0.0,1.0E3);      // 
    // beta_n_waterhole ~ normal(0.0,1.0E3);      // 
    // int_site ~ normal(0,sigma_site);//
    // k ~ normal(0.0,100);
    tau ~ normal(5000,1000);

     for (i in 1:N) {
        detected[i] ~ bernoulli_logit(exp(-(time[i])/tau) * n_waterhole[i]);
        //detected[i] ~ bernoulli_logit(k * exp(-(time[i])/tau) * n_waterhole[i]);
        // detected[i] ~ bernoulli_logit(alpha + k * exp(-(time[i])/tau) * n_waterhole[i]);
    }

/*    for (i in 1:N) {
        detected[i] ~ bernoulli(p[i]);
    }
*/    

}

