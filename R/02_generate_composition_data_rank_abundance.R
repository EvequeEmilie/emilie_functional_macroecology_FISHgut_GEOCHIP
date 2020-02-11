################################################################################
# ANALYZE THE DITRIBUTION OF FUNCTIONAL GENES AND FUNCTIONAL POTENTIAL         #
#  OF FISH MICROBIAL COMMUNITIES                                               #
#                                                                              #
# ------------------------------------------------------------------------------
# arthur.escalas@gmail.com 
# ------------------------------------------------------------------------------
#
#	- classify genes into abundance quantiles, n = c(5, 10, 15, 20)              #
# - generate functional composition table for every abundance quantile         #
# - this is done for observed data and matrices with randomized abundance      #
#                                                                              #
################################################################################


################################################################################
#                                                                              #
#  LOAD DATA                                                                   #
#                                                                              #
################################################################################

# Load the common label matrix

geochip <- read.csv(paste0(dir_res_01, "labels_geochip.csv"), row.names = 1)

# Load the table matching different functional levels

table_match_funct_levels <- read.csv(paste0(dir_res_01, 'table_match_funct_levels.csv'))

# Load the abundance data 

table_data_abundance_log <- read.csv(paste0(dir_res_01, "table_data_abundance_log.csv"),
                                     row.names = 1)

# the directory to save the results
dir_save <- dir_res_02


################################################################################
#                                                                              #
# GENERATE COMPOSITION TABLES FOR ABUNDANCE QUANTILES                          #
#                                                                              #
################################################################################

# ==============================================================================
# Determine the distribution of genes into abundance quantiles                 #
# ==============================================================================

# Loop on various number of quantiles

lapply(c(5, 10, 15, 20), function(num_q) {
  
  mat_gene_per_quantiles <- getAbQuantId(table_data_abundance_log,
                                 method = "trueQuant",
                                 nullModel = FALSE,
                                 num_of_quantiles = num_q)
  
  saveRDS(mat_gene_per_quantiles, 
          file = paste0(dir_save, "matrices_gene_per_quantile_q", num_q,".rds"))
})


# ==============================================================================
# Associate function with genes for each functional level                      #
# ==============================================================================

lapply(c(5, 10, 15, 20), function(num_q) {

  mat <- readRDS(paste0(dir_save, "matrices_gene_per_quantile_q", num_q, ".rds"))
    
  tables_funct_compo <- lapply(names_funct_levels, function(funct_lev) {
  getCompoFromQuantilTable(mat, abundance, CLM = geochip, funcLevel = funct_lev, 
                                    abWeight = FALSE, scaled = TRUE)

  })
  
  saveRDS(tables_funct_compo, file = paste0(dir_save, "tables_funct_compo_q", 
                                            num_q, ".rds"))
})



################################################################################
#                             END OF SCRIPT                                    #
################################################################################
