################################################################################
# ANALYZE THE DITRIBUTION OF FUNCTIONAL GENES AND FUNCTIONAL POTENTIAL
#  ACROSS SOIL MICROBIAL COMMUNITIES
#
# ------------------------------------------------------------------------------
# arthur.escalas@gmail.com 
# ------------------------------------------------------------------------------
#
#	SCRIPT TO :
# - characterize the frequency-occupancy distribution of genes across communities
#
################################################################################


################################################################################
#
#  LOAD DATA  
#
################################################################################

# Load the common label matrix (CLM)
CLM <- read.csv(paste0(dir_res_01, "table_common_label_matrix.csv"), row.names = 1)

# Load the raw data list
data_list <- readRDS(paste0(dir_res_01, "data_list_raw.rds"))

# Load metadata
metadata <- read.csv(paste0(dir_res_01, "metadata.csv"), row.names = 1)

# the directory to save the results

dir_save <- dir_res_05


################################################################################
#
#	  CHARACTERIZE THE FREQUENCY-OCCUPANCY DISTRIBUTIONS
#
################################################################################

# ==============================================================================
#  Load data
# ==============================================================================

occupancy_data <- readRDS(paste0(dir_res_03, "ls_matrices_occupancy_data.rds"))
mat_occupancy  <- read.csv(paste0(dir_res_03, "matrice_gene_occupancy.csv"), 
                           row.names = 1)

# ==============================================================================
#   Test bimodality using tokeshi test (vegan package)
#   as done in Lindh 2016 Env Micr
#    Metapopulation theory identifies biogeographical patterns among 
#     core and satellite marine bacteria scaling from tens to thousands of kilometers
# ==============================================================================

# USING QUANTILES --------------------------------------------------------------

# Test the relationship between gene frequency and rank across all sites 
#  and using different numbers of quantiles

res_tokeshi_global <- lapply(occupancy_data, function(X) { MOStest(X$rank, X$funct_gene) })

table_tokeshi_global <- lapply(res_tokeshi_global, function(x) x$coefficients) %>% 
  reformat_as_df(new_var_name = "num_of_quantile")

write.csv(table_tokeshi_global, file = paste0(dir_save, 'res_tokeshi_test_global.csv'))

# extract coefficients for plotting
# (Intercept) / x / I(x^2)

coeff_tokeshi_test <-  do.call(rbind, lapply(res_tokeshi_global, function(x) {
  attributes(profile(x))$summary$coefficients[,1]
}))
colnames(coeff_tokeshi_test) <- letters[1:3]

write.csv(coeff_tokeshi_test, file = paste0(dir_save, 'coeff_tokeshi_test_global.csv'))

# Test the relationship in each site using 10 occupancy quantiles

res_tokeshi_site <- apply(mat_occupancy, 2, function(X) {
  h <- hist(X[X!=0], 10)
  x <- h$count[h$counts != 0]
  return(list(test = MOStest(1:length(x), x),
              data = x))
})

table_tokeshi_site <- do.call(rbind, lapply(res_tokeshi_site, function(X) {
  X$test$coefficients
})) %>% rownames_to_column(var = "rownames") %>% separate(rownames, 
                                                          into = c("site", "min/max"),
                                                          sep= "\\.")
write.csv(table_tokeshi_site, 
          file = paste0(dir_save, 'table_tokeshi_test_ecosystem_10_quantiles.csv'))


# USING QUANTILES --------------------------------------------------------------

X <- as.vector(mat_occupancy)
x <- X[X!=0]
res <- MOStest(1:length(table(x)), table(x))
res$coefficients

# Test the relationship in each site 

res_tokeshi_site_q <- apply(mat_occupancy, 2, function(X) {
  x <- X[X!=0]
  res <- MOStest(1:length(x), x)
  res <- MOStest(1:length(table(x)), table(x))
  return(res)
})

table_tokeshi_site <- do.call(rbind, lapply(res_tokeshi_site_q, function(X) {
  X$coefficients
})) %>% rownames_to_column(var = "rownames") %>% separate(rownames, 
                                                          into = c("site", "min/max"),
                                                          sep= "\\.")

write.csv(table_tokeshi_site, 
          file = paste0(dir_save, 'table_tokeshi_test_ecosystem_all_ranks.csv'))

# ==============================================================================
#   Plot the global relationship
# ==============================================================================

num_q <- "q10"
dat <- occupancy_data[[num_q]]

png(paste0(dir_save, "plot_occupancy_frequency_global_model.png"),
    height = 12, width = 15, unit = "cm", res = 600)
par(mar = c(4, 4, 2, 2), oma=c(0,1,0,0))

ylim    <- setrge(c(dat$funct_gene), 10)
plotTokeshi(res_tokeshi_global[[num_q]], colPch = "grey30", 
            colLine = 'black', xlim = c(0.5, 10 + 0.5), ylim = c(0, max(ylim)), 
            xlab = "", ylab = "", cex=0.5)
axis(1, 1:10, labels =  1:10)
mtext(2, text="Gene frequency", line = 3, font = 2, 
      cex = 1, outer = F, las = 0)
mtext(1, text="Occupancy quantile", line = 2, font = 2, cex = 1, outer = F)

graphics.off()



png(paste0(dir_save, 'plot_occupancy_frequency_within_ecosystems.png'), 
    height = 20, width = 15, unit = 'cm', res = 400)

par(mfrow=c(5,2),oma=c(3,4,1,1), mar=c(2,2,1,1), las=1, mgp=c(3,0.3,0), tcl=-0.2) 

for (nm_site in names_sites) {
  
  restest <- res_tokeshi_site[[nm_site]]$test
  coeffs  <- restest$coefficients
  tab     <- res_tokeshi_site[[nm_site]]$data
  ylim    <- setrge(c(tab), 20)
  plotTokeshi(restest, main_text = names_sites_cplte[[nm_site]],
              colPch = "grey30", 
              colLine = 'black', xlim = c(0.5, length(tab) + 0.5), 
              ylim = c(0, max(ylim)), 
              xlab = "", ylab = "", cex=1)
  axis(1, 1:10, labels =  1:10, cex.lab = 0.1)
}
par(las=0)
mtext(2, text = "Gene frequency", line=1.5, font=2, cex=1, outer=T)
mtext(1, text = "Occupancy quantile", line=1, font=2, cex=1, outer=T)

graphics.off()



