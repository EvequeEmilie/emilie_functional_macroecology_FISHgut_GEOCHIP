################################################################################
# ANALYZE THE DITRIBUTION OF FUNCTIONAL GENES AND FUNCTIONAL POTENTIAL
#  WITHIN AND ACROSS SOIL MICROBIAL COMMUNITIES
#
# ------------------------------------------------------------------------------
# arthur.escalas@gmail.com 
# ------------------------------------------------------------------------------
#
#	- estimate functional dissimilarity across quantiles (occupancy or abundance)
# - make an ordination of quantiles
# - this is done for observed and randomized data
#
################################################################################


################################################################################
#
#  LOAD DATA  
#
################################################################################

# Load the common label matrix
CLM <- read.csv(paste0(dir_res_01, "table_common_label_matrix.csv"), row.names = 1)

# Load the table matching different functional levels
table_match_funct_levels <- read.csv(paste0(dir_res_01, 'table_match_funct_levels.csv'))

# load composition data 
num_q <- "q10"

data_compo_ab <- readRDS(paste0(dir_res_02, "tables_funct_compo_", num_q, ".rds"))
data_compo_oc <- readRDS(paste0(dir_res_03, "ls_tables_function_composition.rds"))
data_compo_oc <- data_compo_oc[[num_q]]

# the directory to save the results
dir_save <- dir_res_07

################################################################################
#
#	   ANALYZE DISSIMILARITY BETWEEN ABUNDANCE QUANTILES IN EACH ECOSYSTEM
#   
################################################################################

# dissimilarity method used
diss_meth <- "jaccard"


# ==============================================================================
# Estimate dissimilarity within ecosystem
# ==============================================================================

diss_within_site_abundance <- lapply(data_compo_ab, function(funct_lev) {  
  lapply(funct_lev, function(site) {                
    vegdist(t(do.call(cbind, site)  ), diss_meth)
  })
})

saveRDS(diss_within_site_abundance, 
        file = paste0(dir_save, "diss_within_site_abundance_", diss_meth, ".rds"))


# ==============================================================================
# do the DCA
# ==============================================================================

dca_abundance <- lapply(data_compo_ab, function(funct_lev) {  
  lapply(funct_lev, function(site) {                
    site <- do.call(cbind, site)  
    decorana(t(site))
  })
})

saveRDS(dca_abundance, file = paste0(dir_save, "dca_abundance_", diss_meth, ".rds"))


# ==============================================================================
#     Plot the ordination of abundance quantiles within ecosystems
# ==============================================================================

for (X in 1:3) {
  
  funct_lev      <- names_funct_levels[[X]]
  funct_lev_name <- names_funct_levels_plots[X]
  
  png(paste0(dir_save, "plot_DCA_within_ecos_abundance_", diss_meth, "_", funct_lev,".png"), 
      height = 20, width = 15, unit = "cm", res = 200)
  
  par(mfrow = c(5,2), mar = c(2,1,.5,1.5), oma = c(3,4,3,1), mgp = c(3,.3,0), 
      las = 1, xaxs = 'i', yaxs = "i", tcl = -0.2)
  
  lapply(names_sites, function(nm_site) {
    toPlot <- summary(dca_abundance[[funct_lev]][[nm_site]])$site.scores[, c(1, 2)]
    cex <- as.vector(rep(2.5 * seq(2, 15, length.out=10) / 10, nrow(toPlot) / 10))
    xlim <- round(setrge(toPlot[,1]),2)
    ylim <- round(setrge(toPlot[,2]),2)
    plot(toPlot, pch = 21, bg = paste0(colors_sites[nm_site], 20),
          col = paste0(colors_sites[nm_site], 60), cex = cex, main = '', 
         xlab = '', ylab = '', xlim = xlim, ylim = ylim)
    text(apply(toPlot, 2, function(x) tapply(x, cex, mean)), 
         labels = 1:10, cex = 1, font = 2, adj = 0.5, col = "black")
    mtext(nm_site, 3, font = 2, cex = 0.8)
  })
  mtext(funct_lev_name,3, font = 4, cex = 1, out = T, line = 1.5)
  mtext('DCA axis 1', 1, font = 2, cex = 1, out = T)
  mtext('DCA axis 2', 2, font = 2, las = 0, cex = 1, out = T, line = 1.5)
  
  graphics.off()
}

# ==============================================================================
# Test effect of quantile and sample within ecosystem
# ==============================================================================

res_permanova_abundance_within <- lapply(names_funct_levels, function(funct_lev) {  
  lapply(names_sites, function(nm_site) {
    X <- diss_within_site_abundance[[funct_lev]][[nm_site]]
    fact_quant <- factor(rep(1:10, nrow(as.matrix(X)) / 10))
    fact_sple  <- factor(sapply(1:10, rep, vect_num_sple_per_site[nm_site]))
    adonis(X ~ factor_quant + fact_sple, permutations = 999)$aov.tab
  })
})

saveRDS(res_permanova_abundance_within, 
        file = paste0(dir_save, "res_permanova_abundance_within_", diss_meth, ".rds"))

lapply(names_funct_levels, function(funct_lev) {
  res_permanova_table <- do.call(rbind, 
                                 lapply(res_permanova_abundance_within[[funct_lev]], 
                                        function(X) X[1:3,]))
  write.csv(res_permanova_table, 
            file = paste0(dir_save, "res_permanova_within_abundance_", funct_lev, 
                          diss_meth, ".csv"))
})



################################################################################
#
#	  ANALYZE DISSIMILARITY BETWEEN OCCUPANCY QUANTILES IN EACH ECOSYSTEM
#   
################################################################################

# create factors for analysis

factor_samples <- factor( do.call(c, sapply(1:10, function(x) {
  rep(names_sites[x], vect_num_sple_per_site[x])
})), levels = names_sites)
factor_sites <- factor(sapply(names_sites, rep, 10), levels = names_sites)
factor_quantiles <- factor(rep(paste0("Q", 1:10), 10), levels = paste0("Q", 1:10))


# ==============================================================================
# Estimate dissimilarity between occupancy quantiles within ecosystem
#	  - for all pairs of quantiles
#   - for each quantile to the Q10 (widespread gene pool)
# ==============================================================================

# ------------------------------------------------------------------------------
#   Merge all the quantiles into one table (10quant x 10sites)
# ------------------------------------------------------------------------------

# rename the columns
data_compo_oc <- lapply(data_compo_oc, function(X) { 
  lapply(X, function(x) { 
    colnames(x) <- paste0("Q", 1:10)
    x
  })
})

data_compo_oc_all_sples <- lapply(data_compo_oc, function(X) {
  out <- do.call(cbind, X)
  colnames(out) <- paste0(factor_sites, "_", factor_quantiles)
  return(out)
})

# Add the functional genes level

gene_quantile_mat <- readRDS(paste0(dir_res_03, 'res_matrices_genes_per_quantiles.rds'))
tmp <- do.call(cbind, gene_quantile_mat[[num_q]])
colnames(tmp) <- paste0(factor_sites, "_", factor_quantiles)
data_compo_oc_all_sples$funct_gene <- tmp

# ------------------------------------------------------------------------------
# Estimate dissimilarity between occupancy quantiles within ecosystems
# ------------------------------------------------------------------------------

diss_within_site<- list()
diss_within_site$full_matrix <- list()
diss_within_site$diss_to_q10 <- list()

for (funct_lev in names_funct_levels) {
  
  res <- lapply(names_sites, function(nm_site) {
    vegdist(t(data_compo_oc[[funct_lev]][[nm_site]]), diss_meth) 
  }) %>% setNames(names_sites)
  diss_within_site$full_matrix[[funct_lev]] <- res
  tmp <- do.call(c, lapply(res, function(x) as.matrix(x)[10,-10] ))
  tmp2 <- factor(sapply(names(tmp), function(x) strsplit(x, "\\.")[[1]][1]),
                 levels = levels(factor_sites))
  tmp3 <- factor(sapply(names(tmp), function(x) strsplit(x, "\\.")[[1]][2]),
                 levels = levels(factor_quantiles))
  diss_within_site$diss_to_q10[[funct_lev]] <- data.frame(diss_q10 = tmp,
                                                   site = tmp2,
                                                   quantile = tmp3,
                                                   rank = as.numeric(tmp3),
                                                   row.names = 1:length(tmp))
}

saveRDS(diss_within_site, 
         file = paste0(dir_save, 'diss_occupancy_quantiles_within_site.rds'))

# ------------------------------------------------------------------------------
# Pairwise test between quantiles
# In pairwise we test only the quantile effect as testing the site effect makes 
#  no sense because quantiles from a given site do not have any gene in common
#  by definition 
# The pairwiseAdonis() function originates from:
# https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
# ------------------------------------------------------------------------------

res_permanova_pairwise_occupancy <- lapply(data_compo_oc_all_sples, function(X) {
  tmp <- pairwiseAdonis(t(X), factor_quantiles, stratum = NULL, 
                            sim.function = 'vegdist', sim.method = 'bray', 
                            p.adjust.m = 'fdr', reduce = NULL, perm = 999)
})

lapply(names_funct_levels_all, function(x) {
  write.csv(res_permanova_pairwise_occupancy[[x]], 
            file = paste0(dir_save, 'res_permanova_pairwise_occupancy_', x,'.csv'))
})

# ==============================================================================
#  Compare the strength of the difference between quantiles
# ==============================================================================

# ------------------------------------------------------------------------------
# Compare the F values for each pairwise comparison to Q10
# ------------------------------------------------------------------------------

fval_tab <- do.call(cbind, lapply(res_permanova_pairwise_occupancy, function(X) {
  X[grep(10, X$pairs), 'F.Model']
}))

png(paste0(dir_save, "plot_Fvalues_pairwise_permanova_occupancy.png"), 
    height = 20, width = 20, unit = "cm", res = 200)

par(mar = c(6, 5, 1, 1))
plot(NA, type = 'n', xlim = c(.75, 9.25), ylim = c(5, 40), xaxt = 'n',
     ylab = '', xlab = '', las = 1)
lapply(1:4, function(x) {
  lines(1:9, fval_tab[,x], lty = x,lwd = 2)
  points(1:9, fval_tab[,x], pch = 21, bg = 'grey', cex = 1.5)
})
tmp <- res_permanova_pairwise_occupancy[[1]]
text(1:9, rep(2,9), labels = tmp[grep(10, tmp$pairs),'pairs'], srt = 45, 
     cex = 0.8, xpd = NA, adj = c(1, 0.5))
axis(1,1:9, labels = rep('', 9))
mtext('F value of pairwise PERMANOVA', 2, line = 3, font = 2)
mtext('Pair of occupancy quantiles', 1, line = 4, font = 2)
legend('topright', legend = names_funct_levels_all_plots, 
       lty = 1:4, lwd = 2, cex = 0.8)

graphics.off()
