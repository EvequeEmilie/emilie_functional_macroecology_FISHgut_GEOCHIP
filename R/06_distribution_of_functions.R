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

geochip <- read.csv(paste0(dir_res_01, "labels_geochip.csv"), row.names = 1)

# Load the table matching different functional levels

table_match_funct_levels <- read.csv(paste0(dir_res_01, 'table_match_funct_levels.csv'))

# Load the abundance data 

table_data_abundance_log <- read.csv(paste0(dir_res_01, "table_data_abundance_log.csv"),
                                     row.names = 1)

# Load metadata

metadata <- read.csv(paste0(dir_res_01, "metadata.csv"))

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
  
  names_functions <- sort(unique(geochip[, funct_lev]))
  
  # loop on functions
  tmp <- do.call(rbind, lapply(names_functions, function(f) { 
      dat <- tables_funct_compo[[funct_lev]]
      funct_tab   <- do.call(rbind, lapply(dat, function(tab) { tab[f, ] }))
      funct_tab <-  gather(as.data.frame(funct_tab)) %>%
        mutate(Sample = rep(row.names(funct_tab), num_q),
               Rank   = as.numeric(sapply(key, function(x) {
                 strsplit(x, "Q")[[1]][2]
               })),
               Funct = rep(f, num_q * length(dat))
        )
    }))
  tmp[, c("Sample", "Rank", "Funct", "value")]
}) %>% setNames(names_funct_levels)

saveRDS(data_function_quantile, 
        file = paste0(dir_save, "data_function_quantile_q", num_q, ".rds"))
}


################################################################################

# FIT MODELS TO DESCRIBE VARIATION OF FUNCTION ABUNDANCE ACROSS RARITY GRADIENT

# Fit models to the relationship between rank and function abundance:
#
# 1. First we fit a simple linear model to describe the relationship between the
#    rank along the rarity-to-abundance gradient and the weight of the function
################################################################################

num_q <- 5

data_function_quantile <- readRDS(paste0(dir_save, "data_function_quantile_q", num_q, ".rds"))

# ==============================================================================
#  Model relation between abundance quantile and abundance for each function
# ==============================================================================

# ------------------------------------------------------------------------------
#  Fit liner model between rank and function weight in each community
# ------------------------------------------------------------------------------

# fit the models for functional traits -----------------------------------------

funct_lev <- names_funct_levels[3]

dat <- data_function_quantile[[funct_lev]]

res_lm_fit <- mclapply(split(dat, dat$Funct), mc.cores = 2, function(x_funct) {
  fit <- lapply(split(x_funct, x_funct$Sample), function(x_sple) {
    tmp <- x_sple %>% filter(value != 0)
    if (nrow(tmp) <= 3) {
      # do not try to fit a model if there is 3 values or less 
      res <- NA
    } else {
      res <- lm(value ~ Rank, data = tmp)
    }
    res
  })
  fit[! is.na(fit)]  # remove the models that were not fitted
})

res_lm_fit <- res_lm_fit[lapply(res_lm_fit, length) != 0]

saveRDS(res_lm_fit, file = paste0(dir_save, "res_lm_fits_rank_vs_weight_q", num_q, "_", 
                                  funct_lev, ".rds"))


# Make some plots to visualize the relationships -------------------------------

# for each function make a plot of the model for a random subset of 64 samples

dir_mods <- paste0(dir_save, "plots_model/", funct_lev, "/")
dir_mods_plot <- paste0(dir_mods,num_q)
dir.create(dir_mods)
dir.create(dir_mods_plot)


for (X in 1:length(res_lm_fit)) {
  
  nm <- names(res_lm_fit)[X]
  png(paste0(dir_mods_plot, "/plot_model_weight_Rank_", nm, ".png" ),
      height = 30, width = 30, unit = "cm", res = 200)
  
  par(mfrow = c(6,4), mar = c(2,2,1,1), oma = c(1,1,1,1))
  for (x in 1:length(res_lm_fit[[X]])) {
    
    fit <- res_lm_fit[[X]][[x]]
    if (is.na(fit)) {
      plot(1,1, main = names(res_lm_fit[[X]])[x])
    } else {
      plot(value ~ Rank, data = fit$model, main = names(res_lm_fit[[X]])[x])
      abline(a = coefficients(fit)[1], b = coefficients(fit)[2], col = "red")
    }
  }
  mtext(text = nm, out = TRUE, side = 3)
  graphics.off()
}


# take a look at the models R2

dir_r2 <- paste0(dir_save, "plots_rsquare/", funct_lev, "/")
dir_r2_plot <- paste0(dir_r2,num_q)
dir.create(dir_r2)
dir.create(dir_r2_plot)


lapply(list(1:25, 26:50, 51:75, 76:100, 101:125, 
            126:150, 151:175, 176:194), function(X) {
              
              png(paste0(dir_r2_plot, "/plot_model_rsquare_", X[1], "_", X[length(X)], ".png" ),
                  height = 30, width = 50, unit = "cm", res = 200)
              
              par(mfrow = c(5,5), mar = c(2,2,1,1), oma = c(1,1,1,1))
              
              for (x in X) {
                
                nm <- names(res_lm_fit)[x]
                r2 <- do.call(c, lapply(res_lm_fit[[x]], function(fit) {
                  if (! is.na(fit)) { 
                    summary(fit)$r.squared 
                  } else {
                    NA
                  }
                }))
                if (sum(is.na(r2)) != length(r2)) {
                  r2 <- r2[!is.na(r2)]
                  hist(r2, 20, xlim = c(0,1), main = nm)
                  abline(v = mean(r2), col = "red")
                  abline(v = median(r2), col = "blue")
                }
              }
              graphics.off()
            })



# reformat the results ---------------------------------------------------------
# make a table containing the results of lm for each sample and each function

# res_lm_fit <- readRDS(paste0(dir_save, "res_model_fits_rank_vs_weight_q", num_q, "_", 
#                       funct_lev, ".rds"))

res_lm_tab <- lapply(res_lm_fit, function(x_funct) {
  
  res <- do.call(rbind, lapply(x_funct[!is.na(x_funct)], function(fit) {
    out <- c(as.vector(summary(fit)$coefficients),
             round(summary(fit)$r.squared, 4),
             round(summary(fit)$adj.r.squared, 4),
             round(summary(fit)$fstatistic,4))
    names(out) <- c("intercept", "slope", "sd_intercept", "sd_Rank", 
                    "tval_intercept", "tval_Rank",
                    "pval_intercept", "pval_Rank", "rsquare", 
                    "rsquare_adj", "F_value", "df_num", "df_den")
    out
  })) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Sample") %>% 
    left_join(metadata, by = "Sample")
  res
}) %>% reformat_as_df(new_var_name = funct_lev) %>% 
  left_join(table_match_funct_levels, by = funct_lev)

saveRDS(res_lm_tab, file = paste0(dir_save, "res_table_lm_rank_vs_weight_q",
                                  num_q, "_", funct_lev, ".rds"))

write.csv(res_lm_tab, file = paste0(dir_save, "res_table_lm_rank_vs_weight_q", 
                                    num_q, "_", funct_lev, ".csv"), row.names = FALSE)

# 
dim(res_lm_tab)

res_lm_tab %>% filter(pval_Rank < 0.05) %>% dim()

do.call(rbind, lapply(seq(0, 1, 0.1), function(x) {
  out <- res_lm_tab %>% filter( rsquare > x) %>% nrow()
  c(out, round(out / nrow(res_lm_tab), 1))
})) %>% as.data.frame %>% 
  setNames(c("num_models_with_r2_>_thd", "prop_models_with_r2_>_thd")) %>% 
  mutate(r2_thd = seq(0, 1, 0.1))



# ==============================================================================
#  Analyse the factors influencing the slope and intercept of linear models
#   describing the distribution of functions along the abundance gradient
#   (i.e. the relationship between quantiles and weight)
# ==============================================================================

dat <- res_lm_tab %>% filter(funct_process == "Halogenated aliphatic compounds")

boxplot(dat$slope ~ dat$Type)



