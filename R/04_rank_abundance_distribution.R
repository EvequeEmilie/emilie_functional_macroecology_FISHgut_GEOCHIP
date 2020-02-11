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

metadata <- read.csv(paste0(dir_res_01, "metadata.csv"))

# the directory to save the results

dir_save <- dir_res_04


################################################################################
#
# GENERATE DATA FOR RANK-ABUNDANCE ANALYSIS
#
################################################################################

table_rank_abundance_sum <- apply(abundance, 1, sum)

# # Create a table with the mean of colonne data abundance 
# 
# table_rank_abundance_mean <- abundance %>% summarise_all(funs(mean))
# 
# #Create table with the sum of colonne data abundance 
# 
# table_rank_abundance_sum <- abundance %<% summarise(nb=n())

# sum the abundance of genes in each habitat

table_gene_abundance_habitat <- do.call(cbind, lapply(levels(metadata$Type), function(x) {
  apply(abundance[, metadata$Type == x], 1, sum)
})) 
colnames(table_gene_abundance_habitat) <- names_habitats


# Save data 
saveRDS(table_rank_abundance_sum, 
        file = paste0(dir_save, "table_rank_abundance_sum.rds"))
saveRDS(table_gene_abundance_habitat, 
        file = paste0(dir_save, "table_rank_abundance_sum_habitats.rds"))

# ==============================================================================
#  Fit several RAD models using macroecotools

# "https://github.com/weecology/macroecotools/archive/master.zip"
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

dir.create(paste0(dir_save, "radfit/", "samples/" ))

# transpose data

rad_tab <- t(abundance)
row.names(rad_tab) <- names(abundance)

# Loop on each sample and fit RAD models

mclapply(row.names(rad_tab), function(i) {
  x <- rad_tab[i, ]
  cat("fitting RAD models on community", i, "\n")
  fit <- we_rad_fit(x)
  save(fit, file = file.path(dir_radfit, "samples", i))
}, mc.cores = 7)

# ------------------------------------------------------------------------------
# for each sample type
# ------------------------------------------------------------------------------

dir.create(paste0(dir_res_04, "radfit/", "habitats/" ))

fit_habitat <- lapply(names_habitats, function(x) {
  we_rad_fit(table_gene_abundance_habitat[, x])
})

lapply(names(fit_habitat), function(f){
  obj <- fit_habitat[[f]]
  save(obj, file = file.path(dir_radfit, "habitats", f))
})

# ------------------------------------------------------------------------------
# for the whole ecosystem
# ------------------------------------------------------------------------------

fit <- we_rad_fit(table_rank_abundance_sum)
save(fit, file = file.path(dir_radfit, "ecosystem"))


# ==============================================================================
# Analysing the results
# ==============================================================================

# ------------------------------------------------------------------------------
# for each sample
# ------------------------------------------------------------------------------

# Load the fitted models into a list

s_files <- list.files(paste0(dir_res_04, "radfit/", "samples/" ))

rad_list <- lapply(s_files, function(x) {
  load(paste0(dir_res_04, "radfit/", "samples/", x))
  fit
}) %>% setNames(s_files)

# save(rad_list, file = paste0(dir_res_04, "rad_list.R"))

# Get a table of models aic, delta_aics and aic_weights for each community

list_aic_samples <- lapply(rad_list, function(fit) {
  
  AIC_values <- unlist(fit$AICc)
  names(AIC_values) <- rad_models
  
  aic_comp(AIC_values)
  
})

# Get the best model for each community

names_best_models <- do.call(c, lapply(list_aic_samples, function(x) {
  rownames(x)[which(x$best_mod == 1)]
}))

# Count the number of times each model is selected in each site

counts_mod <- tapply(names_best_models, metadata$Site, function(x) {
  table(factor(x, levels = unique(names_best_models)))
}) %>% bind_rows() %>% t() %>% as.data.frame() #%>% setNames(c("Zipf","Lognormal"))

write.csv(counts_mod, file = paste0(dir_res_04, "table_counts_of_best_models.csv"))

# Make a table with all the models coefficients

tab_mod_coeff <- do.call(rbind, lapply(rad_list, function(X) {
  X$parameters[[2]] %>% setNames(c("mu","sigma"))
})) %>% as.data.frame() %>% rownames_to_column(var = "Sample") %>% 
  left_join(metadata, by = "Sample")

write.csv(tab_mod_coeff, file = paste0(dir_res_04, 
                                       "table_model_coeff_per_community.csv"))


# ------------------------------------------------------------------------------
# for each habitat
# ------------------------------------------------------------------------------

h_files <- list.files(paste0(dir_res_04, "radfit/", "habitats/" ))

rad_list <- lapply(h_files, function(x) {
  load(paste0(dir_res_04, "radfit/", "habitats/", x))
  obj
}) %>% setNames(h_files)

# Get a table of models aic for each community

list_aic_habitat <- lapply(rad_list, function(fit) {
  
  AIC_values <- unlist(fit$AICc)
  names(AIC_values) <- rad_models
  
  aic_comp(AIC_values)
  
})

# Get the best model for each community

names_best_models <- do.call(c, lapply(list_aic_habitat, function(x) {
  rownames(x)[which(x$best_mod == 1)]
}))

# Make a table with all the models coefficients

tab_mod_coeff <- do.call(rbind, lapply(rad_list, function(X) {
  X$parameters[[2]] %>% setNames(c("mu","sigma"))
})) %>% as.data.frame() %>% rownames_to_column(var = "Site")

write.csv(tab_mod_coeff, file = paste0(dir_res_04, 
                                       "table_model_coeff_per_habitat.csv"))

# ------------------------------------------------------------------------------
# Ecosystem level
# ------------------------------------------------------------------------------

load(file.path(dir_radfit, "ecosystem")) # object name = "fit"








