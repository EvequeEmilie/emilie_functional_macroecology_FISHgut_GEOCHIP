################################################################################
# ANALYZE THE DITRIBUTION OF FUNCTIONS ALONG THE ABUNDANCE-OCCUPANCY GRAIDENT
#
# ------------------------------------------------------------------------------
# arthur.escalas@gmail.com 
# ------------------------------------------------------------------------------
#
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
table_match_funct_levels <- read.csv(paste0(dir_res_01, 'table_match_funct_levels.csv'),
                                     row.names = 1)

# Load abundance - occupancy data
list_ab_occ_within <- readRDS(paste0(dir_res_08,
                                     "data_abundance_occurrence_within_sites.rds"))

tab_ab_occ_across <- readRDS(paste0(dir_res_08, 
                                    "table_abundance_occurrence_across_sites.rds"))

# Load metadata
metadata <- read.csv(paste0(dir_res_01, "metadata.csv"), row.names = 1)
metadata_site <- read.csv(paste0(dir_res_01, "metadata_site.csv"), row.names = 1)

# dafine save directory
dir_save <- dir_res_09

# ==============================================================================
#   Compute the functional potential of core and satellites genes 
# ==============================================================================

# ------------------------------------------------------------------------------
#  Within ecosystems
# ------------------------------------------------------------------------------

list_compo_within <- lapply(names_funct_levels, function(funct_lev){
  
  lapply(names_sites, function(nm_site) {
    
    dat <- list_ab_occ_within[[nm_site]]
    groups <- dat$ab_occ_q %>% as.vector()
    names(groups) <- rownames(dat)
    
    # make a gene x quantile matrix
    tplte <- matrix(0, nrow(CLM), 10)
    rownames(tplte) <- rownames(CLM)
    
    for (i in 1:10) {
      ids <- names(which(groups == i))
      tplte[ids, i] <- groups[ids]
    }
    tplte[tplte != 0] <- 1
    colnames(tplte) <- paste0("gp_", 1:10)
    
    # get the composition
    compo_tab <- getCompoFromQuantilTable(tplte, originalTable = NA, CLM,
                                         funcLevel = funct_lev, 
                                         abWeight = FALSE, scaled = TRUE)
    out <- do.call(cbind, compo_tab)
    colnames(out) <- colnames(tplte)
    return(out)
  })
})


# Reformat the data and add Rank, Site and Function name

num_q <- 10

data_within_site <- lapply(names_funct_levels, function(funct_lev) {
  
  names_functions <- sort(unique(CLM[, funct_lev]))
  
  # loop on functions
  tmp <- do.call(rbind, lapply(names_functions, function(f) { 
    
    # loop on sites
    funct_tab <- do.call(rbind, lapply(names_sites, function(nm_site) {
      
      dat <- list_compo_within[[funct_lev]]
      a   <- do.call(rbind, lapply(dat, function(tab) { tab[f, ] }))
      gather(as.data.frame(a)) %>%
        mutate(Site   = rep(names_sites, num_q),
               Rank   = as.numeric(sapply(key, function(x) {
                 strsplit(x, "gp_")[[1]][2]
               })),
               Funct = rep(f, num_q)
        )
    }))
  }))
  # tmp[, "value"] <- tmp[, "value"] - 1 
  val <- tmp[, "value"]
  val[val == 0] <- NA
  tmp[, "value"] <- val
  tmp[, c("Site", "Rank", "Funct", "value")]
}) %>% setNames(names_funct_levels)

# data_within_site <- lapply(data_within_site, function(X) {
#   X$Site <- as.factor(do.call(c, X$Site))
#   left_join(X, metadata_site, "Site")
# })
# 
data_within_site$funct_process <- data_within_site$funct_process %>%
  mutate(funct_process = Funct) %>%
  left_join(table_match_funct_levels, "funct_process")

data_within_site$funct_trait <- data_within_site$funct_trait %>%
  mutate(funct_trait = Funct) %>%
  left_join(table_match_funct_levels, "funct_trait")


saveRDS(data_within_site, file = paste0(dir_save, "data_within_site.rds"))

# ==============================================================================
# Fit models to the distribution abundance
# ==============================================================================

model_list <-  list(mod1 = "value ~ Rank + (1 | Site)",
                    mod2 = "value ~ Rank + I(Rank^2) + (1 | Site)",
                    mod3 = "log(value) ~ Rank + (1 | Site)"
                    )

model_list <-  list(mod1 = "value ~ Rank",
                    mod2 = "value ~ Rank + I(Rank^2)",
                    mod3 = "log(value) ~ Rank"
                    )
# ------------------------------------------------------------------------------
# Using data across sites
# ------------------------------------------------------------------------------

list_model_fits <- lapply(data_within_site[3], function(X) {
  lapply(split(X, X$Funct), function(x) {
    x$Site <- factor(x$Site, levels = names_sites)
    fit_lm(x, model_list, nm_var_Y = NULL)
  })
})

saveRDS(list_model_fits, file = paste0(dir_save, "results_model_fits_within.rds"))


# Compute estimates of values in Q1 to Q10 using evey models

list_predicted_values <- lapply(list_model_fits, function(X) {
  lapply(X, function(x) {
    out <- lapply(x$mods, function(xx) { predict(xx, list(Rank = 1:num_q))})
    out$mod3 <- exp(out$mod3)
    out
  })
})

saveRDS(list_model_fits, file = paste0(dir_save, "results_model_fits_within.rds"))


# Extract the results of model fit for all models ------------------------------
# name of the best model
# parameters of the best model
# p value of the parameters
# F and p values of the ANOVA on main effects
# difference in AIC between the lm and the lmm
# differences in AIC, R2m and R2c between both lmm models (mod1 vs mod2)
# AIC, R3m and R2c of the best model

list_best_models <- lapply(list_model_fits, function(X) {
  out <- do.call(rbind, lapply(X, extract_lm_results, force_linear = FALSE)) %>% 
    as.data.frame()
  out
})

# add informations on other classification levels

table_best_models <- list_best_models$funct_trait %>% as.data.frame() %>% 
  rownames_to_column("funct_trait")

table_best_models <- table_best_models %>%
  left_join(table_match_funct_levels, "funct_trait") %>% 
  select(names_funct_levels, everything())


# Save the results

saveRDS(table_best_models, file = paste0(dir_save, "results_best_model_within.rds"))


# ==============================================================================
# Analyse results of the models
# ==============================================================================

# Determine the shape of the relationship
# ------------------------------------------------------------------------------

model_shape <- do.call(rbind, lapply(table_best_models$funct_trait, function(func) {
  determine_model_shape(table_best_models %>% filter(funct_trait == func), 
                        x_pred = list_predicted_values$funct_trait[[func]],
                        pval_thd = 0.05, intercept_thd = 1)
                  })) %>% as.data.frame()

table_best_models$Type  <- factor(model_shape$V1, levels = c("No relation", "Linear", 
                                              "Quadratic", "Exponential"))
table_best_models$Shape <- factor(model_shape$V2, levels = c("No relation", "U shaped", "Bell shaped",
                                              "Negative", "Positive"))

data_model <- table_best_models %>% arrange(Shape, funct_category, mod1_Rank)

# Analyse the results

table(data_model$Type)
table(data_model$Type) / 194 *100

table(data_model$Shape)
table(data_model$Shape) / 194 *100


res_shape <- do.call(cbind, lapply(split(data_model, data_model$Shape), function(x) { 
  table(x$funct_category)
}))

res_shape_pct <- apply(res_shape, 2, function(x) round(x / sum(x) * 100, 1))

write.csv(res_shape, file = file.path(dir_save, "res_model_shape.csv"))
write.csv(res_shape_pct, file = file.path(dir_save, "res_model_shape_pct.csv"))
write.csv(data_model, file = file.path(dir_save, "res_table_models.csv"))



res_shape_pct <- apply(res_shape, 1, function(x) round(x / sum(x) * 100, 1))
barplot(res_shape_pct)


# Make panel plot for each type of relationship --------------------------------


png(paste0(dir_save, 'plot_models_positive.png'), 
    height = 30, width = 40, unit = "cm", res = 200)

model_res <- data_model %>% filter(Shape %in% c("Positive")) %>% 
  arrange(Shape, funct_category, mod1_Rank)
raw_data <- data_within_site$funct_trait
predictions <- list_predicted_values$funct_trait

plot_pars <- list(mfrow = c(8,10), mar = c(1,1,2,1), las = 1, oma = c(3,3,2,1), xpd = NA,
    mgp = c(1.5,0.5,0))

fill_panel_plot(model_res, raw_data, predictions, plot_pars)

graphics.off()

png(paste0(dir_save, 'plot_models_no_relation_U_and_Bell_shape.png'), 
    height = 30, width = 40, unit = "cm", res = 200)

model_res <- data_model %>% filter(Shape %in% c("No relation", "U shaped", "Bell shaped")) %>% 
  arrange(Shape, funct_category, mod1_Rank)
raw_data <- data_within_site$funct_trait
plot_pars <- list(mfrow = c(8,10), mar = c(1,1,2,1), las = 1, oma = c(3,3,2,1), xpd = NA,
                  mgp = c(1.5,0.5,0))

fill_panel_plot(model_res, raw_data, predictions, plot_pars)

graphics.off()


png(paste0(dir_save, 'plot_models_negative.png'), 
    height = 40, width = 50, unit = "cm", res = 200)

model_res <- data_model %>% filter(Shape %in% c("Negative")) %>% 
  arrange(Shape, funct_category, mod1_Rank)
raw_data <- data_within_site$funct_trait
plot_pars <- list(mfrow = c(9,10), mar = c(1,1,2,1), las = 1, oma = c(3,3,2,1), xpd = NA,
                  mgp = c(1.5,0.5,0))

fill_panel_plot(model_res, raw_data, predictions, plot_pars)

graphics.off()







