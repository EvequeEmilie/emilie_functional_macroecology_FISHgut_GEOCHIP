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
  
  p_list <- list(p_untruncated,
                 p_pln,
                 p_negbin,
                 p_zipf)
  
  AICc_list <- list(AICc_logser_untruncated,
                    AICc_pln,
                    AICc_negbin,
                    AICc_zipf)
  
  likelihood_list <- list(L_logser_untruncated,
                          L_pln,
                          L_negbin,
                          L_zipf)
  
  relative_likelihood_list <- list(relative_ll_logser_untruncated,
                                   relative_ll_pln,
                                   relative_ll_negbin,
                                   relative_ll_zipf)
  
  res <- list(parameters = p_list,
              AICc = AICc_list,
              likelihood = likelihood_list,
              relative_likelihood = relative_likelihood_list)
  
  res
  
}#eo we_rad_fit
