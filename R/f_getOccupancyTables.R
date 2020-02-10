################################################################################
# This function associates each gene/OTU from a community to an occupancy quantile
#
# input: 
#  - a presence/absence matrix with genes/OTUs in row and communities in columns
#  - the common label  matrix matching genes and functions
# output: 
#  a list containing:
#   tot_num_of_genes
#   gene_occurrence: for each gene the number of times it is detected
#   gene_occupancy: for each gene the proportino of samples it is detected
#   gene_X_quantile_matrix: presence-absence of genes into quantiles
#   gene_id_per_occurence_level: ids of genes in each level of occurrence 
#   results_occurrence: table of number and proportion of gene per occurrence level
#   results_occurrence_quantil: table of number and occupancy of gene per occurrence level
#
################################################################################

getOccupancyTables <- function(dataset, CLM, num_quants = 10){
  
  # Get the number of genes in the dataset
  tot_num_of_genes <- sum(dataset)
  
  # get the number of samples in dataset
  num_of_sples <- ncol(dataset)

  # Get the occurence of each gene in the dataset
  gene_occurrence <- apply(dataset, 1, function(x) sum(x != 0))  

  # Get the occupancy for each gene
  gene_occupancy <- gene_occurrence / num_of_sples

  # Get the indexes of genes for each occurence level
  #  this is used later to match genes with their functions
  gene_id_per_occurence_level <- lapply(0:ncol(dataset), function(oc) {  
    as.vector(names(which(gene_occurrence == oc))) 
    })
  names(gene_id_per_occurence_level) <- paste0("occ_", 0:ncol(dataset))

  # Calculate occurrence statistics
  # ----------------------------------------------------------------------------
  
  # Number of genes in each occurence rank
  vect_to_fill <- rep(0, ncol(dataset) + 1)
  names(vect_to_fill) <- 0:ncol(dataset)
  vect_to_fill[names(table(gene_occurrence))] <- table(gene_occurrence)
  num_of_genes <- as.integer(vect_to_fill)
  
  #	Proportion of genes in each occurence rank
  prop_of_genes <- round(vect_to_fill / sum(vect_to_fill) * 100, 2)
  
  # Aggregate results into quantiles
  # ----------------------------------------------------------
  
  # test whether there is more than 10 samples
  if (num_of_sples < num_quants) {
    num_quants <- num_of_sples
  }

  tplte <- matrix(0, nrow(dataset), num_quants)
  quantiles_limits <- seq(0, 1, length.out = num_quants + 1)					

  for (i in 1:(length(quantiles_limits)-1))  {
    
    lims <- c(quantiles_limits[i], quantiles_limits[i + 1])
    for (j in 1:length(gene_occupancy)) {
    
      val <- gene_occupancy[j]
      if (lims[1] < val && val <= lims[2]) {
        tplte[j, i] <- 1
      } else {
        tplte[j, i] <- 0
      }
    }
  } 
  colnames(tplte)  <- paste0("Q", 1:num_quants)
  row.names(tplte) <- 1:nrow(tplte)
  gene_X_quantile_matrix <- tplte

  out <- list(tot_num_of_genes = tot_num_of_genes, 
              gene_occurrence  = gene_occurrence,
              gene_occupancy   = gene_occupancy, 
              gene_X_quantile_matrix      = gene_X_quantile_matrix, 
              gene_id_per_occurence_level = gene_id_per_occurence_level, 
              results_occurrence          = rbind(num_of_genes, prop_of_genes), 
              results_occurrence_quantil  = 
                rbind(num_of_genes  = colSums(gene_X_quantile_matrix),
                      prop_of_genes = round(colSums(gene_X_quantile_matrix) /
                                               length(gene_occurrence[gene_occurrence != 0]) * 100, 2)
                      )
              )
 return(out)
  
}
##eo function










