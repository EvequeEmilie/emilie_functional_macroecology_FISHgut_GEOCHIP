################################################################################
# ANALYZE THE DITRIBUTION OF FUNCTIONAL GENES AND FUNCTIONAL POTENTIAL         #
#  OF FISH MICROBIAL COMMUNITIES                                               #
#                                                                              #
# ------------------------------------------------------------------------------
# arthur.escalas@gmail.com                       
# ------------------------------------------------------------------------------
#                                                                              #
# Project setup script                                                         # 
#                                                                              #
################################################################################

############################# CLEANUP & CO #####################################
rm(list = ls())                                                                #
################################################################################



################################################################################
# PATHS                                                                        # 
################################################################################

# Define where is the project folder

dir_media <- getwd()

# for data, scripts and results

dir_data   <- paste0("data/")
dir_res    <- paste0("res/")
dir_script <- paste0("R/")

dir.create(dir_data)
dir.create(dir_res)

################################################################################
#	DEFINE WORKING DIRECTORIES                                                   #
################################################################################

# for each result folder

dir_list <- list( dir_res_01 = paste0(dir_res, "01_clean_and_normalize_data/"),
                  dir_res_02 = paste0(dir_res, "02_generate_composition_data_rank_abundance/"),
                  dir_res_03 = paste0(dir_res, "03_generate_composition_data_frequency_occupancy/"),
                  dir_res_04 = paste0(dir_res, "04_rank_abundance_distribution/"),
                  dir_res_05 = paste0(dir_res, "05_frequency_occupancy_distribution/"),
                  dir_res_06 = paste0(dir_res, "06_distribution_of_functions/"),
                  dir_res_07 = paste0(dir_res, "07_functional_dissimilarity/"),
                  dir_res_08 = paste0(dir_res, "08_abundance_occupancy_relationship/"),
                  dir_res_09 = paste0(dir_res, "09_function_distribution_abundance_occupancy/"),
                  dir_res_10 = paste0(dir_res, "10_/")
)

lapply(dir_list, dir.create, showWarnings = FALSE)
lapply(names(dir_list), function(x) {
  assign(x, dir_list[[x]], envir = .GlobalEnv)
})

################################################################################
# FUNCTIONS & LIBRAIRIES                                                       #
################################################################################

source(paste0(dir_script, "f_utils.R"))
source(paste0(dir_script, "f_L_S.R"))
source(paste0(dir_script, "f_getAbQuantId.R"))
source(paste0(dir_script, "f_getCompoFromQuantilTable.R"))
source(paste0(dir_script, "f_randomizeThisTable.R"))
source(paste0(dir_script, "f_we_rad_fit.R"))
source(paste0(dir_script, "f_getOccupancyTables.R"))
source(paste0(dir_script, "f_pairwise_permanova.R"))
source(paste0(dir_script, "f_lm_and_lmm.R"))


libs <- c("vegan", "cluster", "cclust", "ca", "dunn.test", "tidyverse", "labdsv",
          "nlstools", "MuMIn", "nlstools", "lme4", "plotrix", "gambin", "abind",
          "Rmisc", "sads", "RADanalysis", "lmerTest", "parallel", "epitools",
          "reticulate")


load_libraries(libs)

################################################################################
# Weecology macroecotools source files                                         #
################################################################################

## get files and extract

# mt_url <- "https://github.com/weecology/macroecotools/archive/master.zip"
# dir.create(file.path("libs"))
# dir.create(file.path("libs","macroecotools"))
# mt_file <- file.path("libs","macroecotools","master.zip")
# download.file(url = mt_url, destfile = mt_file)
# unzip(mt_file, exdir = file.path("libs","macroecotools"))
# 
# ## python modules
# py_install("matplotlib")
# py_install("pandas")
# py_install("scipy")

## source libraries

reticulate::source_python(file.path("libs", "macroecotools",
                                    "macroecotools-master",
                                    "macroecotools",
                                    "macroecotools.py"))

reticulate::source_python(file.path("libs", "macroecotools",
                                    "macroecotools-master",
                                    "macroeco_distributions",
                                    "macroeco_distributions.py"))


################################################################################
# GLOBAL VARIABLES                                                             #
################################################################################

# names of RAD models
rad_models <- c("Logseries untruncated", "Poisson lognormal", "Negative binomial", "Zipf distribution")
rad_models_short <- c("logser", "pln", "negbin", "zipf")
names(rad_models_short) <- rad_models

# names of different habitats

names_habitats <- c("fish", "sediment", "water")
names(names_habitats) <- c("fish", "sediment", "water")

# names of functional levels
names_funct_levels <- c("funct_category", "funct_process", "funct_trait")
names_funct_levels_plots <- c("Functional categories", 
                              "Functional processes", 
                              "Functional traits")
names(names_funct_levels_plots) <- names(names_funct_levels) <- names_funct_levels

names_funct_levels_all <- c("funct_category", "funct_process", 
                            "funct_trait", "funct_gene")
names_funct_levels_all_plots <- c("Functional categories", 
                                  "Functional processes", 
                                  "Functional traits",
                                  "Functional genes")
names(names_funct_levels_all_plots) <- names_funct_levels_all


# Define category colors
colors_category <- c("#7435C7","#E21712","#ED6615","#DAD40D","#43BB1B",
                     "#1C841D", "#06B0BC","#2948D1", "#14122B")
names(colors_category) <- c("Antibiotic resistance", "Carbon cycling", 
                            "Energy process", "Metal resistance", 
                            "Nitrogen cycling", "Phosphorus cycling",
                            "Stress",  "Sulphur cycling", "Virulence")


################################################################################
#                             END OF SCRIPT                                    #
################################################################################
