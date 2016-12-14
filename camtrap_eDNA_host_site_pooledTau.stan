// Animal detection via waterhole eDNA
// detected per site and photo data  as independent data frames
// multiple sites and species, but no random effects, just pooled decay rate

// MJ Farrell

data {
    
    // Groups
    int<lower=0> SP;                    // number of species
    int<lower=0> SI;                    // number of sites

    // Photo data
    int<lower=0> N_photo;               // # observations
    int<lower=1> species_photo[N_photo];// map photo observations to species (group)
    int<lower=1> site_photo[N_photo];   // map photo observations to site (group)
    
    // Time variable
    real time[N_photo];                 // time coded as integer (minute) before sampling

    // Waterhole Visitation
    real n_waterhole[N_photo];          // number of individuals seen for each photo

    // eDNA data        
    int<lower=0> N_edna;                // # eDNA detection observations (site x species)
    int<lower=1> site_edna[N_edna];     // map edna detection to site (group)
    int<lower=1> species_edna[N_edna];  // map edna detection to species (group)
    
    // Reponse
    int<lower=0, upper=1> detected[N_edna]; // detected in eDNA (binary 0-1)
    
}

parameters {

    // parameters for exponential weighting of n_waterhole by time
    // real k;                // scale parameter
    // maybe put upper bound on k to help limit nan return...?

    real<lower=0.0001> tau;           // mean lifetime (in units of time)

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

   vector[N_photo] total_edna;          // amount of edna expected to be leftover at sampling
   for (q in 1:N_photo) {               //  based on decay rate and n_waterhole (per photo)
        total_edna[q] = (exp(-(time[q])/tau) * n_waterhole[q]);
    }


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

    // This is accepted...
    // real sum_edna;
    // sum_edna = sum(total_edna);

    // This is accepted...
    // real sum_edna;
    // sum_edna = sum(total_edna[1:20]);
 
    // This is accepted...
    //  vector[SP] sum_edna;
    //  for (i in 1:SP){
    //      sum_edna[i] = total_edna[species_photo[i]];
    //  }
    
    // This is accepted...
    vector[SP] sum_edna;
    for (i in 1:SP){
        sum_edna[i] = total_edna[species_photo[i]];
    }
    
    // BUT TAKING sum(total_edna[species_photo[i]]) fails...


    // HOW DO YOU DO THIS INDEXING?
    // WANT TO GET SUM OF "edna" per species per site then use this to model detected...     
    for (i in 1:N_edna) {

        detected[i] ~ bernoulli_logit(sum_edna);
        //detected[i] ~ bernoulli_logit(k * exp(-(time[i])/tau) * n_waterhole[i]);
        // detected[i] ~ bernoulli_logit(alpha + k * exp(-(time[i])/tau) * n_waterhole[i]);
    }

    tau ~ normal(5000,1000);


}

