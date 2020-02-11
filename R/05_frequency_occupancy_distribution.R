################################################################################
# ANALYZE THE DITRIBUTION OF FUNCTIONAL GENES AND FUNCTIONAL POTENTIAL
#  ACROSS SOIL MICROBIAL COMMUNITIES
#
# ------------------------------------------------------------------------------
# arthur.escalas@gmail.com 
# ------------------------------------------------------------------------------
#
#	SCRIPT TO :
# - characterize the frequency-occupancy distribution of genes across communities
#
################################################################################


################################################################################
#
#  LOAD DATA  
#
################################################################################

# Load the common label matrix (CLM)
geochip <- read.csv(paste0(dir_res_01, "labels_geochip.csv"))

# Load metadata
metadata <- read.csv(paste0(dir_res_01, "metadata.csv"), row.names = 1)

# the directory to save the results

dir_save <- dir_res_05


################################################################################
#
#	  CHARACTERIZE THE FREQUENCY-OCCUPANCY DISTRIBUTIONS
#
################################################################################

# ==============================================================================
#  Load data
# ==============================================================================

occupancy_data <- readRDS(paste0(dir_res_03, "data_occupancy.rds"))


# ==============================================================================
#   Test bimodality using tokeshi test (vegan package)
#   as done in Lindh 2016 Env Micr
#    Metapopulation theory identifies biogeographical patterns among 
#     core and satellite marine bacteria scaling from tens to thousands of kilometers
# ==============================================================================

# Test the relationship between gene frequency and rank 

dat <- data_occupancy$results_occurrence[1,-1]

res_tokeshi <- MOStest(1:length(dat), dat)

dat <- data_occupancy$results_occurrence_quantil[1,]

res_tokeshi_q <- MOStest(1:length(dat), dat)


coeff_tokeshi <- list( attributes(profile(res_tokeshi))$summary$coefficients[,1],
                       attributes(profile(res_tokeshi_q))$summary$coefficients[,1]
)

# ==============================================================================
#   Plot the global relationship
# ==============================================================================

num_q <- "q12"
dat <- data_occupancy$results_occurrence_quantil

png(paste0(dir_save, "plot_occupancy_frequency_global_model.png"),
    height = 12, width = 15, unit = "cm", res = 600)
par(mar = c(4, 4, 2, 2), oma=c(0,1,0,0))

ylim    <- setrge(c(dat), 10)
plotTokeshi(res_tokeshi_q, colPch = "grey30", 
            colLine = 'black', xlim = c(0.5, 12 + 0.5), ylim = c(0, max(ylim)), 
            xlab = "", ylab = "", cex=0.5)
axis(1, 1:12, labels =  1:12)
mtext(2, text="Gene frequency", line = 3, font = 2, 
      cex = 1, outer = F, las = 0)
mtext(1, text="Occupancy quantile", line = 2, font = 2, cex = 1, outer = F)

graphics.off()


