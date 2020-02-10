################################################################################
# ANALYZE THE DITRIBUTION OF FUNCTIONAL GENES AND FUNCTIONAL POTENTIAL
#  WITHIN AND ACROSS SOIL MICROBIAL COMMUNITIES
#
# ------------------------------------------------------------------------------
# arthur.escalas@gmail.com 
# ------------------------------------------------------------------------------
#
#	- fit the global abundance-occupancy relationship
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

# Load the logged data list
data_list <- readRDS(paste0(dir_res_01, "data_list_logged.rds"))

# the directory to save the results

dir_save <- dir_res_08

# ------------------------------------------------------------------------------
#   Prepare abundance data
# ------------------------------------------------------------------------------

# Estimate summed gene abundance within ecosystems

table_sum_gene_abundance <- do.call(cbind, lapply(data_list, function(X) { 
  apply(X, 1, sum) 
}))

# Estimate average gene abundance within ecosystems

table_avg_gene_abundance <- do.call(cbind, lapply(data_list, function(X) { 
  apply(X, 1, mean) 
}))

# Estimate average gene abundance across ecosystems
ab_across  <- apply(table_sum_gene_abundance, 1, mean)

# ------------------------------------------------------------------------------
#   Prepare occupancy data
# ------------------------------------------------------------------------------

# Load gene occupancy within ecosystems
mat_gene_occupancy <- read.csv(paste0(dir_res_03, "matrice_gene_occupancy.csv"),
                               row.names = 1)

# Estimate gene occupancy across ecosystems
avg_occ_across <- apply(mat_gene_occupancy, 1, mean)

# ------------------------------------------------------------------------------
#  Merge occupancy and abundance data
# ------------------------------------------------------------------------------

# For each ecosystem -----------------------------------------------------------
# Actual occupancy and average abundance

list_ab_occ_within <- lapply(names_sites, function(nm_site) {
  out <- data.frame(occupancy = mat_gene_occupancy[, nm_site],
                    abundance = table_avg_gene_abundance[, nm_site])
  out[rowSums(out) != 0,]
})

#  Across ecosystems -----------------------------------------------------------
# Actual occupancy and average abundance

tab_ab_occ_across <- data.frame(occupancy = avg_occ_across, abundance = ab_across)
tab_ab_occ_across <- tab_ab_occ_across[rowSums(tab_ab_occ_across) != 0,]


# ==============================================================================
# Classify genes into quantiles of groups
# ==============================================================================

num_q <- 10

# Using quantiles

list_ab_occ_within <- lapply(list_ab_occ_within, function(X) {
  
  # Associate occupancy and abundance values with a quantile and create a new 
  #  grouping corresponding to the sum of quantiles
  
  out <- getAbQuantId(X, method= "trueQuant", nullModel=FALSE, 
                      num_of_quantiles = num_q) %>% 
    as.data.frame() %>% 
    dplyr::rename(occupancy_q = occupancy, abundance_q = abundance) %>%
    mutate(group = occupancy_q + abundance_q) %>% 
    mutate(ab_occ_q = sapply(group, function(x) { 
      if (x %in% seq(2, 20, 2)) { x / 2 } else { (x + 1) / 2 } 
    }))

  cbind(X, out)
})

tmp <- getAbQuantId(tab_ab_occ_across, method= "trueQuant", nullModel=FALSE, 
                    num_of_quantiles = num_q) %>% 
  as.data.frame() %>% 
  dplyr::rename(occupancy_q = occupancy, abundance_q = abundance) %>%
  mutate(group = occupancy_q + abundance_q) %>% 
  mutate(ab_occ_q = sapply(group, function(x) { 
    if (x %in% seq(2, 20, 2)) { x / 2 } else { (x + 1) / 2 } 
}))

tab_ab_occ_across <- cbind(tab_ab_occ_across, tmp)

saveRDS(list_ab_occ_within, 
        file = paste0(dir_save, "data_abundance_occurrence_within_sites.rds"))
saveRDS(tab_ab_occ_across, 
        file = paste0(dir_save, "table_abundance_occurrence_across_sites.rds"))


# ==============================================================================
# Fit models to the abundance-occupancy relationship 
# ==============================================================================

# ------------------------------------------------------------------------------
# For each ecosystem
# ------------------------------------------------------------------------------

list_model_fits <- lapply(list_ab_occ_within, function(X) {
  lm(abundance ~ occupancy, data = X)
})

saveRDS(list_model_fits, file = paste0(dir_save, "results_model_fits_within_sites.rds"))

# ------------------------------------------------------------------------------
#  Across ecosystems
# ------------------------------------------------------------------------------

model_fits <- lm(abundance ~ occupancy, data = tab_ab_occ_across)

saveRDS(model_fits, file = paste0(dir_save, "results_model_fits_across_sites.rds"))


# ==============================================================================
#   Plot the abundance-occurrence relationship 
# ==============================================================================

# ------------------------------------------------------------------------------
#  Within ecosystems
# ------------------------------------------------------------------------------

png(paste0(dir_save, 'plot_abundance_occupancy_within_ecosystems.png'), 
    height = 20, width = 15, unit = 'cm', res = 400)

par(mfrow=c(5,2),oma=c(3,3,1,1), mar=c(2,2,1,1), las=1, mgp=c(3,0.3,0), tcl=-0.2) 

for (nm_site in names_sites) {
  
  tab <- list_ab_occ_within[[nm_site]]
  mod <- list_model_fits[[nm_site]]
  ylim <- c(0, max(setrge(tab$abundance)))
  plot(tab$abundance  ~ tab$occupancy, pch = 21, 
       main = names_sites_cplte[nm_site],
       bg  = "#75757530", #paste0(colors_sites[nm_site], '20') , 
       col = "#75757530", #paste0(colors_sites[nm_site], '20'), 
       ylim = ylim, cex = 0.2, xlim = c(0,1), 
       xaxt = "n")  
  axis(1, seq(0,1,0.1))
  coeff <- mod$coefficients
  addsegment(a = coeff[1], b = coeff[2], limX = c(0,1), 
             limY = setrge(tab$abundance), 
             col_seg = "black", lwd_seg = 1.5)
  text_lm <- equa_lm(mod, nmy = "abundance", nmx = "occupancy", prec = 1)
  text(0.4, max(ylim) - 0.2 * diff(ylim), labels = text_lm, cex=.7)
       
}
par(las=0)
mtext(2, text = "Average gene abundance within site", line=.8, font=2, cex=1, outer=T)
mtext(1, text = "Gene occupancy within site", line=1, font=2, cex=1, outer=T)

graphics.off()


# to visualize the groups

cols <- c("#EF61D1", "#C25BF0", "#7056F1", "#5088F2", "#4BDDF3", "#45F4B0", 
          "#40F550", "#8AF63A", "#EDF734", "#F89A2E", "#F92C28")

png(paste0(dir_save, 'plot_abundance_occupancy_within_ecosystems_groups.png'), 
    height = 20, width = 15, unit = 'cm', res = 400)

par(mfrow=c(5,2),oma=c(3,3,1,1), mar=c(2,2,1,1), las=1, mgp=c(3,0.3,0), tcl=-0.2) 

for (nm_site in names_sites) {
  
  tab <- list_ab_occ_within[[nm_site]]
  groups <- tab$ab_occ_q
  ylim <- c(0, max(setrge(tab$abundance)))
  plot(tab$abundance ~ tab$occupancy, ylim = ylim, cex = 0.2, 
       xlim = c(0,1), main = names_sites_cplte[nm_site], pch = 21,
       bg = cols[groups], col = cols[groups])
}
par(las=0)
mtext(2, text = "Average gene abundance within site", line=.8, font=2, cex=1, outer=T)
mtext(1, text = "Gene occupancy within site", line=1, font=2, cex=1, outer=T)

graphics.off()



# ------------------------------------------------------------------------------
#  Across ecosystems
# ------------------------------------------------------------------------------

df <- tab_ab_occ_across

png(paste0(dir_save, 'plot_abundance_occupancy_across_ecosystems.png'), 
    height = 15, width = 15, unit = 'cm', res = 400)

par(las=1, xaxs = "i", yaxs = "i", mgp = c(2.5,0.5,0), oma = c(1,1,1,1),
    mar = c(3.5,3.5,1,1), tcl = -0.2)
ylim <- setrge(df$abundance)
plot(df, col = "#75757530", 
     bg = "#75757530", pch = 21, 
     xlim = c(-0.05,1.05), ylim = ylim, 
     xlab = "Average occupancy across sites", 
     ylab = "Average abundance across sites", 
     xaxt = 'n', font.lab = 2 , cex.lab = 1.2, cex = .2)
axis(1, seq(0,1,0.1))

coeff <- model_fits$coefficients
addsegment(a = coeff[1], b = coeff[2], limX = c(0,1), 
           limY = setrge(df$abundance), 
           col_seg = "black", lwd_seg = 1.5)
text_lm <- equa_lm(mod, nmy = "abundance", nmx = "occupancy", prec = 1)
text(0.4, max(ylim) - 0.2 * diff(ylim), labels = text_lm, cex=.7)
graphics.off()


