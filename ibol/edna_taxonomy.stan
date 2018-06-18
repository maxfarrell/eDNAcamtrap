// Animal detection via waterhole eDNA
// Basic model - pooling AB replicates
// MJ Farrell

// Last updated November 6th 2017

data {
    int<lower=1> N;                     // # observations

    // Response
    int<lower=0, upper=1> detected[N];  // detected in eDNA (binary 0-1)

    // Continuous predictors
    vector[N] n_waterhole;              // total number of individuals seen across all photos in preceeding week
    vector[N] last_visit;               // time since last visit to site (hours)

    // Hierarchical Predictors
    int<lower=1> M;                     // number of mamal species
    int<lower=1, upper=M> species[N];   // map observations to species (group)
    vector[M] mass;                     // mean adult female body mass (g)

    int<lower=1> F;                     // number of mamal families
    int<lower=1, upper=F> family[M];    // map species to family (group)

    int<lower=1> O;                     // number of mamal families
    int<lower=1, upper=O> order[F];     // map family to order (group)
    
    int<lower=1> N_ID;                  // number of samples (site-times)    
    int<lower=1, upper=N_ID> ID[N];     // map observations to sample ID (group)
 
    int<lower=1> S;                     // number of sites
    int<lower=1, upper=S> site[N_ID];   // map sample ID to sites (group)
    

    // env data to be added ...
/*    vector[N_ID] temp;                  // mean temperature at time of sampling
    vector[N_ID] pH;                    // mean pH at time of sampling
    vector[N_ID] mS;                    // mean conductivity (mS) at time of sampling
    vector[N_ID] DO;                    // mean dissolved oxygen (DO) at time of sampling
*/    
}

parameters {

    real beta0;                         // global intercept
    real beta_n_waterhole;              // n_waterhole coef
    real beta_last_visit;               // last_visit coef    

    // Host
    vector[M] species_raw;              // intercepts for mammal species
    real<lower=0> sigma_species;        // variance for species intercepts
    real beta_mass;                     // coef for mass

    vector[F] family_raw;              // intercepts for mammal families
    real<lower=0> sigma_family;        // variance for family intercepts

    vector[O] order_raw;              // intercepts for mammal orders
    real<lower=0> sigma_order;        // variance for order intercepts

    // Site
    vector[S] site_raw;                 // intercepts for sites
    real<lower=0> sigma_site;           // variance for site intercepts

    // Sample
    vector[N_ID] ID_raw;                // intercepts for samples
    real<lower=0> sigma_ID;             // variance for sample intercepts

}

transformed parameters {

    vector[N_ID] ID_eff;
    vector[M] species_eff; 

    ID_eff = sigma_site * site_raw[site]
            + sigma_ID * ID_raw; 
   
    species_eff = beta_mass * mass
                + sigma_species * species_raw
                + sigma_family * family_raw[family]
                + sigma_order * order_raw[order[family]];

}

model {

    vector[N] p;

    beta0 ~ cauchy(0,10);           
    
    beta_n_waterhole ~ cauchy(0,2.5);     
    beta_last_visit ~ cauchy(0,2.5);     

    site_raw ~ normal(0, 1); 
    sigma_site ~ cauchy(0,2.5); 

    ID_raw ~ normal(0, 1); 
    sigma_ID ~ cauchy(0,2.5);  

    species_raw ~ normal(0, 1); 
    sigma_species ~ cauchy(0,2.5);  
    beta_mass ~ cauchy(0,2.5);     

    family_raw ~ normal(0, 1); 
    sigma_family ~ cauchy(0,2.5);  

    order_raw ~ normal(0, 1); 
    sigma_order ~ cauchy(0,2.5);  
    
    p = beta0 + beta_n_waterhole * n_waterhole 
                + beta_last_visit * last_visit
                + ID_eff[ID] + species_eff[species];

    detected ~ bernoulli_logit(p);
    
}

generated quantities{

    //vector[N] log_lik;
    int detected_pred_exact[N];
    real p_fit[N];

    for (n in 1:N) {

    // Using fit parameter values to estimate p, log-lik & detected_pred_exact
        p_fit[n] = inv_logit(beta0 + beta_n_waterhole * n_waterhole[n] 
                            + beta_last_visit * last_visit[n]
                            + ID_eff[ID[n]] + species_eff[species[n]]);

        // log_lik[n] = bernoulli_lpmf(detected[n] | p_fit[n]);
 
        // Generating predicted effects for "exact" grouping factors
        detected_pred_exact[n] = bernoulli_rng(p_fit[n]);

        }

}
