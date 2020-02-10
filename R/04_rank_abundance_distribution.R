################################################################################
# ANALYZE THE DITRIBUTION OF FUNCTIONAL GENES AND FUNCTIONAL POTENTIAL
# OF FISH MICROBIAL COMMUNITIES
#
# ------------------------------------------------------------------------------
# arthur.escalas@gmail.com 
# ------------------------------------------------------------------------------
#
#	SCRIPT TO :
#	- characterize the rank-abundance distribution of genes into local communities
#
################################################################################


################################################################################
#
#  LOAD DATA  
#
################################################################################

# Load the common label matrix (geochip)

geochip = read.csv(paste0(dir_res_01, "labels_geochip.csv"), row.names = 1)

# Load the raw data list

abundance <- read.csv(paste0(dir_res_01, "rawdata_geochip.csv"))

# Load metadata

metadata <- read.csv(paste0(dir_res_01, "metadata.csv"), row.names = 1)

# the directory to save the results

dir_save <- dir_res_04


################################################################################
#
# GENERATE DATA FOR RANK-ABUNDANCE ANALYSIS
#
################################################################################

table_rank_abundance_mean <- do.call(cbind, lapply(abundance, function(X) { 
  apply(x, 1, function(x) mean(x))
}))

table_rank_abundance_sum <- do.call(cbind, lapply(abundance, function(X) { 
  apply(x, 1, sum) 
}))


# Create a table with the mean of colonne data abundance 

table_rank_abundance_mean <- abundance %>% summarise_all(funs(mean))

#Create table with the sum of colonne data abundance 

table_rank_abundance_sum <- abundance %<% summarise(nb=n())

# Save data 
saveRDS(table_rank_abundance_mean,
        file = paste0(dir_save, "table_rank_abundance_mean.rds"))
saveRDS(table_rank_abundance_sum, 
        file = paste0(dir_save, "table_rank_abundance_sum.rds"))

# ==============================================================================
#  Fit several RAD models using macroecotools

"https://github.com/weecology/macroecotools/archive/master.zip"
#from https://github.com/weecology/sad-comparison/blob/master/sad-comparisons.py

# SAD models and packages used:
#Logseries (macroecotools/macroecodistributions)
#Poisson lognormal (macroecotools/macroecodistributions)
#Negative binomial (macroecotools/macroecodistributions)
#Zipf (macroecotools/macroecodistributions)
# ==============================================================================

# ------------------------------------------------------------------------------
# for each sample
# ------------------------------------------------------------------------------

# create directories to store the results

dir_radfit <- paste0(dir_save, "radfit/")
dir.create(dir_radfit)

# Loop on each sample and fit RAD models

for (i in row.names(rads_tab)) {
  
  x <- rads_tab[i, ]
  cat("fitting RAD models on community", i, "\n")
  fit <- radfit(x, model = c("Null", "Preemption", "Lognormal", "Zipf"))
  saveRDS(fit, file = paste0(dir_radfit, i, ".rds"))
}


# ------------------------------------------------------------------------------
# for each ecosystem
# ------------------------------------------------------------------------------

fit_ecosystem <- lapply(names_sites, function(nm_site) {
  radfit(round(table_sum_gene_abundance[, nm_site]), model = c("Null", "Preemption", "Lognormal", "Zipf"))
})

save(fit_ecosystem, file = paste0(dir_save, "res_radfit_ecosystem_sum"))


# ==============================================================================
# Look at the results
# ==============================================================================

# ------------------------------------------------------------------------------
# for each sample
# ------------------------------------------------------------------------------

# Load the fitted models into a list

files <- rownames(rads_tab)

rad_list <- lapply(files, function(x) {
  load(paste0(dir_radfit, x))
  fit
}) %>% setNames(files)

# save(rad_list, file = paste0(dir_save, "rad_list.R"))

# Get a table of models aic for each community

tab_modaic_rad <- do.call(rbind, lapply(rad_list, function(fit) {
  do.call(c, lapply(fit$models, function(xx) xx$aic))
}))

# Get the best model for each community

names_best_models <- apply(tab_modaic_rad, 1, function(x) {
  names(x)[which(x == min(x, na.rm = TRUE))]
})

# Count the number of times each model is selected in each site

counts_mod <- tapply(names_best_models, factor_samples, function(x) {
  table(factor(x, levels = unique(names_best_models)))
}) %>% bind_rows() %>% t() %>% as.data.frame()%>% setNames(c("Zipf","Lognormal"))

write.csv(counts_mod, file = paste0(res_dir_02, "table_counts_of_best_models.csv"))

# Make a table with all the models coefficients

tab_mod_coeff <- do.call(rbind, lapply(rad_list, function(X) {
  do.call(c, lapply(X$models, coefficients))
})) %>% as.data.frame()

write.csv(tab_mod_coeff, file = paste0(dir_save, 
                                       "tab_mod_coeff_per_community.csv"))

# Analyse effect of metadata on the model parameters ---------------------------

# effect of categorical variables

vars_factor <- names(metadata)[c(1,2,3,8,9)]

res_kw_categvar <- lapply(vars_factor, function(x) {
  kruskal.test(tab_mod_coeff$Lognormal.log.sigma ~ metadata[, x])[1:3]
}) %>% bind_rows() %>%
  mutate(variable = vars_factor) %>% as.data.frame()

# pairwise test between groups

res_dunn_categvar <- lapply(lapply(vars_factor[1:3], function(x) {
  dunn.test(tab_mod_coeff$Lognormal.log.sigma, metadata[, x], altp = TRUE)
}), function(X) {  do.call(cbind, X[c(5,2,3)])}) %>% setNames(vars_factor[1:3])


# effect of continous variables

vars_cont <- names(metadata)[c(7:15)]

fit_contvar <- lapply(vars_cont, function(x) {
  lm(tab_mod_coeff$Lognormal.log.sigma ~ metadata[, x])
}) %>% setNames(vars_cont)

res_aov_contvar <- lapply(fit_contvar, function(x) { summary(aov(x))})


# plot model slope as function of metadata

par(mfrow = c(3,3), mar = c(4,4,1,1), oma = c(1,1,1,1), las = 1, mgp = c(2,0.5,0))

lapply(vars_factor, function(x) {
  boxplot(tab_mod_coeff$Lognormal.log.sigma ~ metadata[, x], 
          ylab = "Model slope", xlab = "")
})
lapply(vars_cont, function(x) {
  plot(tab_mod_coeff$Lognormal.log.sigma ~ metadata[, x],
       ylab = "Model slope", xlab = x)
})

# Maybe keep only the 3 first categorical variables and remove CN

# ------------------------------------------------------------------------------
# for each ecosystem
# ------------------------------------------------------------------------------

# Get a table of models aic for each community

tab_modaic_rad <- do.call(rbind, lapply(fit_ecosystem, function(fit) {
  do.call(c, lapply(fit$models, function(xx) xx$aic))
}))

# Get the best model for each community

names_best_models <- apply(tab_modaic_rad, 1, function(x) {
  names(x)[which(x == min(x, na.rm = TRUE))]
})

# Make a table with all the models coefficients

tab_mod_coeff <- do.call(rbind, lapply(fit_ecosystem, function(X) {
  do.call(c, lapply(X$models, coefficients))
})) %>% as.data.frame()

write.csv(tab_mod_coeff, file = paste0(dir_save, 
                                       "tab_mod_coeff_per_site.csv"))


# ------------------------------------------------------------------------------
# represent the models
# ------------------------------------------------------------------------------

png(paste0(dir_save, 'plot_rank_abundance_within_sites.png'), 
height = 20, width = 15, unit = 'cm', res = 400)

par(mfrow=c(5,2),oma=c(3,3,1,1), mar=c(2,2,1,1), las=1, mgp=c(3,0.3,0), tcl=-0.2) 

for (nm_site in names_sites) {

  tab <- round(table_sum_gene_abundance[, nm_site])
  tab <- sort(tab[tab != 0], decreasing = TRUE)
  bm <- names_best_models[nm_site]
  mod <- fit_ecosystem[[nm_site]]$models[[bm]] 
  ylim <- setrge(log(c(tab, mod$fitted.values)))
  plot(1:length(tab), log(tab), xlim = c(0, 40000), ylim = ylim, 
       main = names_sites_cplte[nm_site], cex = 1, col = "grey")
  lines(log(mod$fitted.values), col = colors_sites[nm_site], lty = 1, lwd = 2)
  legend(25000, max(ylim) - 0.01*diff(ylim), legend = bm, bty = "n", text.font = 3)
  
}
par(las=0)
mtext(2, text = "log(Summed gene abundance in the site)", line=.8, font=2, cex=1, outer=T)
mtext(1, text = "Gene rank", line=1, font=2, cex=1, outer=T)

graphics.off()

png(paste0(dir_save, 'plot_rank_abundance_within_communities.png'), 
    height = 20, width = 15, unit = 'cm', res = 400)

par(mfrow=c(5,2),oma=c(3,3,1,1), mar=c(2,2,1,1), las=1, mgp=c(3,0.3,0), tcl=-0.2) 

for (nm_site in names_sites) {
  
  tab <- round(table_avg_gene_abundance[, nm_site])
  tab <- sort(tab[tab != 0], decreasing = TRUE)
  ylim <- c(3, 14)  #round(setrge(log(tab)))
  plot(NA, xlim = c(0, 34000), ylim = ylim, main = names_sites_cplte[nm_site],
       cex = 1, col = "grey")
  # loop on each community
  lapply(rad_list[grep(nm_site, names(rad_list))], function(mod) {
    dat <- mod$models[[bm]]$fitted.values
    lines(log(dat), col = "#00000050", lty = 1, lwd = 2)
  })
  #legend("topright", legend = bm, bty = "n", text.font = 3)
}
par(las=0)
mtext(2, text = "log(Gene abundance within the community)", line=.8, font=2, cex=1, outer=T)
mtext(1, text = "Gene rank", line=1, font=2, cex=1, outer=T)

graphics.off()
