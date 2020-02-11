################################################################################
# ANALYZE THE DITRIBUTION OF FUNCTIONAL GENES AND FUNCTIONAL POTENTIAL         #
#  OF FISH MICROBIAL COMMUNITIES                                               #
#                                                                              #
# ------------------------------------------------------------------------------
# arthur.escalas@gmail.com 
# ------------------------------------------------------------------------------
#                                                                              #
#	- classify genes into occupancy quantiles                                    #
# - generate functional composition table for every occupancy quantile         #
# - this is done for observed and randomized data                              #
#                                                                              #
################################################################################


################################################################################
#                                                                              # 
#  LOAD DATA                                                                   #
#                                                                              #
################################################################################

# Load the common label matrix

geochip <- read.csv(paste0(dir_res_01, "labels_geochip.csv"))

# Load the table matching different functional levels
table_match_funct_levels <- read.csv(paste0(dir_res_01, 'table_match_funct_levels.csv'))

# Load the presence-absence data abundance
table_abundance_presabs <- read.csv(paste0(dir_res_01, "table_data_abundance_presabs.csv"),
                                    row.names = 1)


# the directory to save the results
dir_save <- dir_res_03


################################################################################
#                                                                              #
#	  DETERMINE OCCUPANCY-FREQUENCY DISTRIBUTION OF GENES                        #
#                                                                              #
################################################################################


# ==============================================================================
# Determine the distribution of genes into occupancy quantiles for each site
#
# Apply the getOccupancyTables() function that compute:
#  tot_num_of_genes
#  gene_occurrence: for each gene the number of times it is detected
#  gene_occupancy: for each gene the proportion of samples in which it is detected
#  gene_per_quantile_matrix: presence-absence of genes into quantiles
#  gene_id_per_occurence_level: ids of genes in each level of occurrence 
#  results_occurrence: table of number and proportion of gene per occurrence level
#  results_occurrence_quantil: table of number and occupancy of gene per occurrence quantile
# ==============================================================================


data_occupancy <- getOccupancyTables(table_abundance_presabs, CLM = geochip, 
                                     num_quants = 12)
saveRDS(data_occupancy, 
        file = paste0(dir_save, "data_occupancy.rds"))


#	Summarize the results  -------------------------------------------------------

# tables of number of genes per occurrence level

res_occurrence <- data_occupancy$results_occurrence_quantil

# tables of proportion of genes per occurrence level

res_occurrence_pct <- data_occupancy$results_occurrence_quantil

# Extract and save the gene_per_quantile matrices

mat_gene_per_quantiles <- data_occupancy[["gene_X_quantile_matrix"]]


# save the objects summarizing the occurrence/occupancy results:  --------------

write.csv(res_occurrence, 
        file = paste0(dir_save, "res_number_of_genes_per_quantiles.csv"))
write.csv(res_occurrence_pct, 
        file = paste0(dir_save, "res_proportion_of_genes_per_quantiles.csv"))
write.csv(mat_gene_per_quantiles, 
        file = paste0(dir_save, "res_matrices_genes_per_quantiles.csv"))


################################################################################
#                                                                              #
# GENERATE COMPOSITION TABLES FOR OCCUPANCY QUANTILES                          #
#                                                                              #
################################################################################

# ==============================================================================
# Associate function with genes for each functional level                      #
# ==============================================================================


for (funct_lev in names_funct_levels) {
    
  res <- do.call(cbind, getCompoFromQuantilTable(mat_gene_per_quantiles, 
                                                     NULL, 
                                                     geochip, 
                                                     funcLevel = funct_lev, 
                                                     abWeight  = FALSE, 
                                                     scaled    = TRUE))
write.csv(res, file = paste0(dir_save, "table_function_composition_",
                           funct_lev, ".csv"))
}


