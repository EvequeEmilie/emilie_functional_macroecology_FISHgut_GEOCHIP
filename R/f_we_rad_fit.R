# This function fits SAd models to abundance data using the weecology tools library
# it used in 04_rank_abundance_distribution.R

we_rad_fit <- function(abund) {
  
  abund <- abund[abund != 0]
  
  # Calculate Akaike weight of species abundance models:
  # Parameter k is the number of fitted parameters
  k1 <- 1
  k2 <- 2            
  
  # Calculate log-likelihoods of species abundance models and calculate AICc values:
  
  S <- length(abund)
  N <- sum(abund)
  
  # Logseries
  p_untruncated <- logser_solver(abund)
  L_logser_untruncated <- logser_ll(abund, p_untruncated) # Log-likelihood of untruncated logseries
  AICc_logser_untruncated <- AICc(k1, L_logser_untruncated, S) # AICc logseries untruncated
  relative_ll_logser_untruncated <- AICc_logser_untruncated# Relative likelihood untruncated logseries
  
  # Poisson lognormal
  p_pln <- unlist(pln_solver(abund))
  L_pln <- pln_ll(abund, p_pln[1], p_pln[2]) # Log-likelihood of Poisson lognormal
  AICc_pln <- AICc(k2, L_pln, S) # AICc Poisson lognormal
  relative_ll_pln <- AICc(k1, L_pln, S) #Relative likelihood, Poisson lognormal
  
  # Negative binomial
  p_negbin <- unlist(nbinom_lower_trunc_solver(abund))
  L_negbin <- nbinom_lower_trunc_ll(abund, p_negbin[1], p_negbin[2]) # Log-likelihood of negative binomial
  AICc_negbin <- AICc(k2, L_negbin, S)# AICc negative binomial
  relative_ll_negbin <- AICc(k1, L_negbin, S) # Relative log-likelihood of negative binomial
  
  # Zipf distribution
  p_zipf <- zipf_solver(abund)
  L_zipf <- zipf_ll(abund, p_zipf) #Log-likelihood of Zipf distribution
  AICc_zipf <- AICc(k1, L_zipf, S)
  relative_ll_zipf <- AICc_zipf
  
  #Making lists
  
  p_list <- list(logseris  = p_untruncated,
                 lognormal = p_pln,
                 binomial  = p_negbin,
                 zipf      = p_zipf)
  
  AICc_list <- list(logseris  = AICc_logser_untruncated,
                    lognormal = AICc_pln,
                    binomial  = AICc_negbin,
                    zipf      = AICc_zipf)
  
  likelihood_list <- list(logseris  = L_logser_untruncated,
                          lognormal = L_pln,
                          binomial  = L_negbin,
                          zipf      = L_zipf)
  
  relative_likelihood_list <- list(logseris  = relative_ll_logser_untruncated,
                                   lognormal = relative_ll_pln,
                                   binomial  = relative_ll_negbin,
                                   zipf      = relative_ll_zipf)
  
  res <- list(parameters = p_list,
              AICc = AICc_list,
              likelihood = likelihood_list,
              relative_likelihood = relative_likelihood_list)
  
  res
  
}#eo we_rad_fit



# from : https://cran.r-project.org/web/packages/blmeco/index.html
aic_comp <- function(AIC_values){
  # AIC_values: vector with AIC values
  delta_aic <- AIC_values - min(AIC_values)
  expdelta_aic <- exp(-0.5*delta_aic)
  aic_weights <- expdelta_aic/(sum(expdelta_aic))
  best_mod <- ifelse(delta_aic == 0, 1, 0)
  return(data.frame(aic = AIC_values, delta_aic = delta_aic, aic_weights = aic_weights, best_mod = best_mod))
}



