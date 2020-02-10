################################################################################
# ANALYZE THE DITRIBUTION OF FUNCTIONAL GENES AND FUNCTIONAL POTENTIAL
#  WITHIN AND ACROSS SOIL MICROBIAL COMMUNITIES
#
# ------------------------------------------------------------------------------
# arthur.escalas@gmail.com 
# ------------------------------------------------------------------------------
#
#	- model functional composition across quantiles (occupancy or abundance)
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

# Load metadata
metadata <- read.csv(paste0(dir_res_01, "metadata.csv"), row.names = 1)

# the directory to save the results

dir_save <- dir_res_06

################################################################################

#   COMPUTE THE ABUNDANCE OF FUNCTION IN EACH ABUNDANCE QUANTILE

################################################################################

# For observed data ------------------------------------------------------------

for (num_q in c(5, 10, 15, 20)) { # loop on the number of quantiles

# load composition data
tables_funct_compo <- readRDS(paste0(dir_res_02, "tables_funct_compo_q", num_q, ".rds"))

data_function_quantile <- lapply(names_funct_levels, function(funct_lev) {
  
  names_functions <- sort(unique(CLM[, funct_lev]))
  
  # loop on functions
  tmp <- do.call(rbind, lapply(names_functions, function(f) { 
    
    # loop on sites
    funct_tab <- do.call(rbind, lapply(names_sites, function(nm_site) {
      
      dat <- tables_funct_compo[[funct_lev]][[nm_site]]
      a   <- do.call(rbind, lapply(dat, function(tab) { tab[f, ] }))
      gather(as.data.frame(a)) %>%
        mutate(Sample = rep(row.names(a), num_q),
               Site   = rep(nm_site, num_q * length(dat)),
               Rank   = as.numeric(sapply(key, function(x) {
                 strsplit(x, "Q")[[1]][2]
               })),
               Funct = rep(f, num_q * length(dat))
        )
    }))
  }))
  tmp[, c("Site", "Sample", "Rank", "Funct", "value")]
}) %>% setNames(names_funct_levels)

saveRDS(data_function_quantile, 
        file = paste0(dir_save, "data_function_quantile_q", num_q, ".rds"))
}

# For permuted data ------------------------------------------------------------

# dir.create(paste0(dir_save, "data/"))
# 
# # loop on functional levels
# for (funct_lev in names_funct_levels) {
#   
#   names_functions <- sort(unique(CLM[, funct_lev]))
#   
#   # loop on functions
#   mclapply(names_functions, mc.cores = 8, function(f) {
#     
#     # loop on sites
#     data_function_quantile_permuted <- do.call(rbind, lapply(names_sites, function(nm_site) {
#       
#       load(paste0(nm_compo_dir, "tables_funct_compo_permuted_", nm_site, ".R"))
#       
#       # loop on permutations
#       do.call(rbind, lapply(1:100, function(rdm) {
#         
#         dat <- tables_funct_compo_permuted[[funct_lev]][[rdm]]
#         a   <- do.call(rbind, lapply(dat, function(tab) { tab[f, ] }))
#         gather(as.data.frame(a)) %>%
#           mutate(Sample = rep(row.names(a), 10),
#                  site   = rep(nm_site, 10 * length(dat)),
#                  Rank   = as.numeric(sapply(key, function(x) {
#                    strsplit(x, "Q")[[1]][2]
#                  })),
#                  rdm    = rep(rdm, 10 * length(dat)),
#                  funct  = rep(f, 10 * length(dat))
#           )
#       }))
#     })) %>% select(site, Sample, Rank, rdm, funct, value)
#     
#     # export the table
#     save(data_function_quantile_permuted, 
#          file = paste0(dir_save, "data/data_function_quantile_abundance_permuted_", funct_lev, "_", f))
#   })
# }
# 

################################################################################

# FIT MODELS TO DESCRIBE VARIATION OF FUNCTION ABUNDANCE ACROSS RARITY GRADIENT

# Fit models to the relationship between rank and function abundance:
#
# 1. First we fit a simple linear model to describe the relationship between the
#    rank along the rarity-to-abundance gradient and the weight of the function
#   We use mixed models and the function lmer() from the package lme4
#   - site and sample are used as random effects
#   - Ranks is used as fixed effect
################################################################################

num_q <- 10

data_function_quantile <- readRDS(paste0(dir_save, "data_function_quantile_q", num_q, ".rds"))

# ==============================================================================
#  Model relation between abundance quantile and abundance for each function
# ==============================================================================

# ------------------------------------------------------------------------------
#  Fit liner model between rank and function abundance in each community
# ------------------------------------------------------------------------------

# fit the models

for (funct_lev in names_funct_levels) {
  
  dat <- data_function_quantile[[funct_lev]]
  res_lm_fit <- mclapply(split(dat, dat$Funct), mc.cores = 2, function(x_funct) {
    fit <- lapply(split(x_funct, x_funct$Sample), function(x_sple) {
      lm(value ~ Rank, data = x_sple)
    })
  })
  
  saveRDS(res_lm_fit, file = paste0(dir_save, "res_model_fits_q", num_q, "_", 
                                    funct_lev, "_abundance.rds"))
}


# reformat the results

for (funct_lev in names_funct_levels) {
  
  dat <- readRDS(paste0(dir_save, "res_model_fits_q", num_q, "_", 
                        funct_lev, "_abundance.rds"))
  
  res_lm_tab <- mclapply(dat, mc.cores = 2, function(x_funct) {
    
    res <- do.call(rbind, lapply(x_funct, function(fit) {
        out <- c(as.vector(summary(fit)$coefficients)[c(1,2,5,6,7,8)],
                 round(summary(fit)$r.squared, 4),
                 round(summary(fit)$adj.r.squared, 4),
                 round(summary(fit)$fstatistic,4))
        names(out) <- c("intercept", "slope", "tval_intercept", "tval_Rank",
                        "pval_intercept", "pval_Rank", "rsquared", 
                        "rsquared_adj", "F_value", "df_num", "df_den")
        out
      })) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "Sample") %>% 
      left_join(metadata, by = "Sample")
      res
    }) %>% reformat_as_df(new_var_name = funct_lev) %>% 
    left_join(table_match_funct_levels, by = funct_lev)
  
  saveRDS(res_lm_tab, file = paste0(dir_save, "res_model_table_", 
                                    funct_lev, "_abundance.rds"))
}


# ==============================================================================
#  Analyse the factors influencing the slope and intercept of linear models
#   describing the distribution of functions along the abundance gradient
# ==============================================================================


# ------------------------------------------------------------------------------
# Define models
#   - Mod 1: only the environment
#   - Mod 2: only functions

# ------------------------------------------------------------------------------

model_list <- list(mod1 = '~ Ecosystem_type + Country + Location + Sample + (1 | Country/Location/Sample)',
                   mod2 = '~ funct_category + funct_process + funct_trait + Sample + (1 | funct_category/funct_process/funct_trait)')




                   mod3 = '~ Rank + Site + (1 | Site/Sample)',
                   mod4 = '~ Rank * Site + (1 | Site/Sample)',
                   mod5 = '~ Rank + Site + I(Rank^2) + (1 | Site/Sample)',
                   mod6 = '~ Rank * Site + I(Rank^2) + (1 | Site/Sample)')


fit <- lm(Rank ~ funct_category + funct_process + funct_trait + Ecosystem_type + Continent + Country + Location + Site,
          data = res_lm_tab)

fit <- lm(Rank ~ funct_category * funct_process * funct_trait * Ecosystem_type * Continent + Country + Location + Site,
          data = res_lm_tab)


summary(fit)
summary(aov(fit))



# ------------------------------------------------------------------------------
#  Fit models and save the best one
# ------------------------------------------------------------------------------

for (funct_lev in names_funct_levels) { 
  
  # Fit models
  
  model_fits <- mclapply(res_lm_tab, mc.cores = 4, function(X) {
                           fit_mixed_models(input_data = X,
                                            model_list = model_list,
                                            nm_var_Y   = "Rank" # slope
                           )
                           
                         })
  saveRDS(model_fits, file = paste0(dir_save, "model_fits_per_function_", funct_lev, ".rds"))
  
  # Find the best model in each case
  
  best_models <-  lapply(model_fits, extract_best_mixmodel)
  saveRDS(best_models, file = paste0(dir_save, "best_models_per_function_", funct_lev, ".rds"))
  
  do.call(rbind, lapply(model_fits, function(x) {
    out <- c(as.vector(summary(x)$coefficients)[c(1,2,5,6,7,8)],
             round(summary(x)$r.squared, 4),
             round(summary(x)$adj.r.squared, 4),
             round(summary(x)$fstatistic,4))
    names(out) <- c("intercept", "Rank", "tval_intercept", "tval_Rank",
                    "pval_intercept", "pval_Rank", "rsquared", 
                    "rsquared_adj", "F_value", "df_num", "df_den")
    out
  }))
}













par(mfrow = c(9,2), mar = c(1,1,1,1), mgp = c(2,0.5,0))

lapply(levels(CLM$funct_category), function(func) {
  
  X <- res_lm_tab$funct_category[[func]]
  hist(X$intercept, 30, xlim = c(0.4,1.5), main = func)
  abline(v = 1)
  hist(X$Rank, 30, xlim = c(-0.15,0.15), main = func)
  abline(v = 0)
})

par(mfrow = c(2,1), mar = c(1,1,1,1), mgp = c(2,0.5,0))

boxplot(res_lm_tab$intercept ~ res_lm_tab$funct_category)
boxplot(res_lm_tab$Rank ~ res_lm_tab$funct_category)



# HERE TEST TO SEE THE EFFECT OF METADATA ON THE PARAMETERS OF THE MODEL RELATING
#  FUNCTION VALUE AND RANK

list_fits <- lapply(names(list_data_function$abundance$funct_category), function(funct) {

data <- list_data_function$abundance$funct_category[[funct]] %>% 
  dplyr::rename(Sple_code = sample)
tab <- left_join(data, metadata, by = "Sple_code") %>% 
  dplyr::rename(sample = Sple_code)
tab$site <- as.factor(tab$site)
tab$sample <- factor(tab$sample, levels = unique(tab$sample))

fits <- lapply(split(tab, tab$sample), function(x) {
  lm(value ~ rank, data = x)
})

mod_coeff <- do.call(rbind, lapply(fits, coefficients)) %>% as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  dplyr::rename(intercept = "(Intercept)", slope = "rank") %>% 
  left_join(tab, by = "sample")

mod_coeff
}) %>% setNames(names(list_data_function$abundance$funct_category))


vars_factor <- names(metadata)[c(1,2,3)]
par(mfrow = c(3,3), mar = c(1,2,1,1), oma = c(1,1,1,1), las = 1)

lapply(list_fits, function(X) {
  lapply(vars_factor[3], function(x) {
    boxplot(X$slope ~ X[,x])
  })
})

vars_cont <- names(metadata)[c(5,6,7,10)]
par(mfrow = c(3,3), mar = c(1,2,2,1), oma = c(1,1,1,1), las = 1)

lapply(names(list_fits), function(f) {
  X <- list_fits[[f]]
  lapply(vars_cont[3], function(x) {
    plot(X$slope ~ X[,x], main = f, col = X$site)
  })
})


################################################################################

# Plot the weight of functions across quantiles

################################################################################

# select the observed data

num_q <- 10
data_function_quantile <- readRDS(paste0(dir_save, "data_function_quantile_q", num_q, ".rds"))


# ==============================================================================
# FOR FUNCTIONAL CATEGORIES
# ==============================================================================

# Prepare data -----------------------------------------------------------------

# select the data to plot

funct_lev <- "funct_category"

# load the randomized data

# files_to_get <-  list.files(paste0(dir_save, "data/"))[grep(paste0(data_type, "_permuted_", funct_lev), list.files(paste0(dir_save, "data/")))]
# 
# ll_p <- lapply(files_to_get, function(f) {
#   load(paste0(dir_save, "data/", f))
#   data_function_quantile_permuted
# }) %>% setNames(unique(table_match_funct_levels[, funct_lev]))


data_obs <- data_function_quantile[[funct_lev]]
ll <- split(data_obs, data_obs$Funct)

# best models

model_fits <- readRDS( paste0(dir_save, "res_model_fits_q", num_q, "_", 
                              funct_lev, "_abundance.rds"))

# Make the plot for functional categories --------------------------------------

png(paste0(dir_save, 'plot_funct_', funct_lev, '.png'), 
height = 20, width = 20, unit = "cm", res = 200)

par(mfrow = c(3,3), mar = c(3,2,1,1), las = 1, oma = c(2,3,1,1))

lapply(names(ll)[c(2,5,6,8,3,1,4,7,9)], function(funct) {

  tab   <- ll[[funct]]
  # tab_p <- ll_p[[funct]]
  ylim  <- setrge(bind_rows(tab)[, "value"], 15)
  xlim  <- c(.5, num_q + 0.5)
  plot(tab$value ~ tab$Rank, type = "n", xaxt = 'n', main = funct, ylim = ylim, xlim = xlim)
  abline(h = 1, lty = 2)
  # boxplot(tab_p[, "value"] ~ tab_p$Rank, col = paste0("#757575", '20'), 
  #          border = paste0("#757575", '60'), main = "", ylim = ylim, 
  #          names = 1:10, add = T, outline = F)
  rk <- sapply(tab$Rank, function(x) x + sample(seq(-0.18, 0.18, 0.001), 1))
  pt_cex <- 0.5
  points(tab$value ~ rk, pch = 21, cex = pt_cex, 
         col = paste0(colors_sites[tab$Site], '40'),
         bg = paste0(colors_sites[tab$Site], '40'))
  # coefs  <- summary(best_models[[funct]])$coefficients[c("(Intercept)", "rank","I(rank^2)"), "Estimate"]
  # lines(1:10, coefs[1] + coefs[2] * (1:10) + coefs[3] * (1:10)^2, 
  #       lwd = 3, col = 'grey20')
  # legend('top', bty = 'n', cex = 1, text.font = 4,
  #        legend = paste('R2=', round_text(r.squaredGLMM(best_models[[funct]])[, "R2m"], 2)))
})
mtext(1, text = "Abundance quantile", outer = T, font = 2, line = 0)
mtext(2, text = "Function weight", outer = T, las = 0, font = 2, line = 1)

graphics.off()


# average function weight per site

ll_avg <- lapply(ll, function(X) { 
  X %>% dplyr::group_by(Site, Rank) %>% dplyr::summarize(value = mean(value))
})

png(paste0(dir_save, 'plot_funct_', funct_lev, '_average_per_site.png'), 
    height = 20, width = 20, unit = "cm", res = 200)

par(mfrow = c(3,3), mar = c(3,2,1,1), las = 1, oma = c(2,3,1,1))

lapply(names(ll_avg)[c(2,5,6,8,3,1,4,7,9)], function(funct) {
  
  tab   <- ll_avg[[funct]]
  # tab_p <- ll_p[[funct]]
  ylim  <- setrge(bind_rows(tab)[, "value"], 15)
  xlim  <- c(.5, num_q + 0.5)
  plot(tab$value ~ tab$Rank, type = "n", xaxt = 'n', main = funct, ylim = ylim, xlim = xlim)
  abline(h = 1, lty = 2)
  rk <- sapply(tab$Rank, function(x) x + sample(seq(-0.18, 0.18, 0.001), 1))
  pt_cex <- 1
  points(tab$value ~ rk, pch = 21, cex = pt_cex, 
         col = paste0(colors_sites[tab$Site], '40'),
         bg = paste0(colors_sites[tab$Site], '40'))
})
mtext(1, text = "Abundance quantile", outer = T, font = 2, line = 0)
mtext(2, text = "Function weight", outer = T, las = 0, font = 2, line = 1)

graphics.off()


# ------------------------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------------------------

# select the data to plot

funct_lev <- "funct_process"

# randomized data

# files_to_get <-  list.files(paste0(dir_save, "data/"))[grep(paste0(data_type, "_permuted_", funct_lev), list.files(paste0(dir_save, "data/")))]
# 
# ll_p <- lapply(files_to_get, function(f) {
#   load(paste0(dir_save, "data/", f))
#   data_function_quantile_permuted
# }) %>% setNames(sort(unique(table_match_funct_levels[, funct_lev])))
# #%>% setNames(sapply(files_to_get, function(x) strsplit(x, paste0(funct_lev, "_"))[[1]][2]))

data_obs <- data_function_quantile[[funct_lev]]
ll <- split(data_obs, data_obs$Funct)

# best models

model_fits <- readRDS( paste0(dir_save, "res_model_fits_q", num_q, "_", 
                              funct_lev, "_abundance.rds"))

# ------------------------------------------------------------------------------
# Make the plot 
# ------------------------------------------------------------------------------

png(paste0(dir_save, 'plot_funct_', funct_lev, '.png'), 
height = 30, width = 30, unit = "cm", res = 200)

par(mfrow = c(7,8), mar = c(3,2,1,1), las = 1, oma = c(2,3,1,1))

lapply(names(ll), function(funct) {

  tab   <- ll[[funct]]
  tab <- tab[tab$value != 0,]
  # tab_p <- ll_p[[funct]]
  ylim  <- setrge(bind_rows(tab)[, "value"], 15)
  xlim  <- c(.5, 10.5)
  plot(tab$value ~ tab$Rank, type = "n", xaxt = 'n', main = funct, ylim = ylim, xlim = c(0.5,10.5))
  abline(h = 1, lty = 2)
  # boxplot(tab_p[, "value"] ~ tab_p$Rank, col = paste0("#757575", '20'), 
  #          border = paste0("#757575", '60'), main = "", ylim = ylim, 
  #          names = 1:10, add = T, outline = F)
  rk <- sapply(tab$Rank, function(x) x + sample(seq(-0.18, 0.18, 0.001), 1))
  pt_cex <- 0.5
  points(tab$value ~ rk, pch = 21, cex = pt_cex, 
         col = paste0(colors_sites[tab$Site], '40'),
         bg = paste0(colors_sites[tab$Site], '40'))
  # coefs  <- summary(best_models[[funct]])$coefficients[c("(Intercept)", "rank","I(rank^2)"), "Estimate"]
  # lines(1:10, coefs[1] + coefs[2] * (1:10) + coefs[3] * (1:10)^2, 
  #       lwd = 3, col = 'grey20')
  # legend('top', bty = 'n', cex = 1, text.font = 4,
  #        legend = paste('R2=', round_text(r.squaredGLMM(best_models[[funct]])[, "R2m"], 2)))
})
mtext(1, text = "Abundance quantile", outer = T, font = 2, line = 0)
mtext(2, text = "Function weight", outer = T, las = 0, font = 2, line = 1)

graphics.off()


# average function weight per site

ll_avg <- lapply(ll, function(X) { 
  X %>% dplyr::group_by(Site, Rank) %>% dplyr::summarize(value = mean(value))
})

png(paste0(dir_save, 'plot_funct_', funct_lev,'_average_per_site.png'), 
    height = 30, width = 30, unit = "cm", res = 200)

par(mfrow = c(7,8), mar = c(3,2,1,1), las = 1, oma = c(2,3,1,1))

lapply(names(ll_avg), function(funct) {
  
  tab   <- ll_avg[[funct]]
  tab <- tab[tab$value != 0,]
  # tab_p <- ll_p[[funct]]
  ylim  <- setrge(bind_rows(tab)[, "value"], 15)
  xlim  <- c(.5, 10.5)
  plot(tab$value ~ tab$Rank, type = "n", xaxt = 'n', main = funct, ylim = ylim, xlim = c(0.5,10.5))
  abline(h = 1, lty = 2)
  # boxplot(tab_p[, "value"] ~ tab_p$Rank, col = paste0("#757575", '20'), 
  #          border = paste0("#757575", '60'), main = "", ylim = ylim, 
  #          names = 1:10, add = T, outline = F)
  rk <- sapply(tab$Rank, function(x) x + sample(seq(-0.18, 0.18, 0.001), 1))
  pt_cex <- 0.5
  points(tab$value ~ rk, pch = 21, cex = pt_cex, 
         col = paste0(colors_sites[tab$Site], '40'),
         bg = paste0(colors_sites[tab$Site], '40'))
  # coefs  <- summary(best_models[[funct]])$coefficients[c("(Intercept)", "rank","I(rank^2)"), "Estimate"]
  # lines(1:10, coefs[1] + coefs[2] * (1:10) + coefs[3] * (1:10)^2, 
  #       lwd = 3, col = 'grey20')
  # legend('top', bty = 'n', cex = 1, text.font = 4,
  #        legend = paste('R2=', round_text(r.squaredGLMM(best_models[[funct]])[, "R2m"], 2)))
})
mtext(1, text = "Abundance quantile", outer = T, font = 2, line = 0)
mtext(2, text = "Function weight", outer = T, las = 0, font = 2, line = 1)

graphics.off()















