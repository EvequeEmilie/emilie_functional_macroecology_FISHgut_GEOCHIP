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

geochip = read.csv(paste0(dir_res_01, "labels_geochip.csv"), row.names = 1)

# Load the table matching different functional levels
table_match_funct_levels <- read.csv(paste0(dir_res_01, 'table_match_funct_levels.csv'), row.names = 1)

# Load the presence-absence data abundance
table_abundance_presabs <- read.csv(paste0(dir_res_01, "table_data_abundance_presabs.cvs"))


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


lapply(c(5, 10, 12, 15, 20), function(num_q) {
  
  data_occupancy <- getOccupancyTables(table_abundance_presabs, CLM = geochip, num_quants = num_q)
  saveRDS(data_occupancy, 
          file = paste0(dir_save, "data_occupancy_q", num_q,".rds"))
}) 



#[data_occupancy <- sapply(c(5,10,12,15,20), quantile)]

#	Summarize the results  -------------------------------------------------------

ls_res_occurrence <- list()
ls_res_occurrence_pct <- list()
ls_gene_per_quantiles_mat <- list()


for( num_q in c(5, 10, 12, 15, 20)) {

   data_occupancy <- readRDS(paste0(dir_save, "data_occupancy_q", num_q,".rds"))
   nm <- paste0("q", num_q)

   # tables of number of genes per occurrence level

   ls_res_occurrence[[nm]] <- do.call(rbind, lapply(data_occupancy, function(X) {
     X$results_occurrence_quantil[1,]
   }))

  # tables of proportion of genes per occurrence level

   ls_res_occurrence_pct[[nm]] <- do.call(rbind, lapply(data_occupancy, function(X) {
     X$results_occurrence_quantil[2,]
   }))

  # Extract and save the gene_per_quantile matrices

   ls_gene_per_quantiles_mat[[nm]] <- lapply(data_occupancy, function(X) {
     X[["gene_X_quantile_matrix"]]
   }) %>% setNames(names_sites)

}


# save the objects summarizing the occurrence/occupancy results:  --------------

saveRDS(ls_res_occurrence, 
        file = paste0(dir_save, "res_number_of_genes_per_quantiles.rds"))
saveRDS(ls_res_occurrence_pct, 
        file = paste0(dir_save, "res_proportion_of_genes_per_quantiles.rds"))
saveRDS(ls_gene_per_quantiles_mat, 
        file = paste0(dir_save, "res_matrices_genes_per_quantiles.rds"))

# Save a table of gene occupancy X sites ---------------------------------------

mat_gene_occupancy <- do.call(cbind, lapply(data_occupancy, function(X) {
  X$gene_occupancy
}))

write.csv(mat_gene_occupancy, file = paste0(dir_save, "matrice_gene_occupancy.csv"))



################################################################################
#                                                                              #
# GENERATE COMPOSITION TABLES FOR OCCUPANCY QUANTILES                          #
#                                                                              #
################################################################################

# ==============================================================================
# Associate function with genes for each functional level                      #
# ==============================================================================

ls_tables_funct_compo <- list()

for (num_q  in names(ls_gene_per_quantiles_mat)) {

  for (funct_lev in names_funct_levels) {
    
    tables_funct_compo[[num_q]][[funct_lev]] <- list()
    
    for (name_site in names_sites) {
      
      gene_per_quant_matrix <- ls_gene_per_quantiles_mat[[num_q]][[name_site]]
      res <- do.call(cbind, getCompoFromQuantilTable(gene_per_quant_matrix, 
                                                     NULL, 
                                                     geochip, 
                                                     funcLevel = funct_lev, 
                                                     abWeight  = FALSE, 
                                                     scaled    = TRUE))
      ls_tables_funct_compo[[num_q]][[funct_lev]][[name_site]] <- res
    }
  }
}


saveRDS(ls_tables_funct_compo, 
        file = paste0(dir_save, "ls_tables_function_composition.rds"))


