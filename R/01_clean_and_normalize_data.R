################################################################################
# ANALYZE THE DITRIBUTION OF FUNCTIONAL GENES AND FUNCTIONAL POTENTIAL         # 
#  OF FISH MICROBIAL COMMUNITIES                                               # 
#                                                                              # 
# ------------------------------------------------------------------------------
# arthur.escalas@gmail.com
# ------------------------------------------------------------------------------
#                                                                              #
#	SCRIPT TO CLEAN AND NORMALIZE DATA FOR DOWNSTREAM ANALYSES                   #
#	                                                                             # 
################################################################################


################################################################################
#                                                                              #
#  LOAD DATA                                                                   #
#                                                                              #
################################################################################

# define the direcctory to save the results

dir_save <- dir_res_01

# Load the common label matrix (geochip)

geochip = read.csv(paste0 ("data/", "labels_geochip.csv"), row.names = 1)


################################################################################
#   CLEAN THE DATASETS                                                         #
# - Removal of genes with less than 10 probes                                  # 
# - Removal of one functional categories "Other category" #
################################################################################

# ==============================================================================
#   Create the mask to remove unwanted genes                                   #
# ==============================================================================


# Identify genes with less than 15 probes
nProb <- do.call(rbind, lapply(split(geochip, geochip$Gene.name), dim))
lessThan15 <- rownames(nProb[nProb[, 1] < 15, ])
maskLessThan15 <- ! geochip$Gene.name %in% lessThan15

# make a mask to remove 'Other category' categories
# maskCateg <- ! geochip$Gene.category %in%  c("Other category")

maskCateg <- ! geochip$Gene.category %in%  c("Other category")

# match the two masks

finalMask <- (maskLessThan15 + maskCateg) == 2 
toRemove <- unique(geochip[! finalMask, 'Gene.name'])



# ==============================================================================
#	GEOCHIP DESIGN                                                               #
# ==============================================================================

geochip <- geochip[finalMask,]
row.names(geochip) <- geochip$Genbank.ID

# Change the levels of each column in the geochip

for(i in names(geochip)[-1]) {
  geochip[,i] <- factor(geochip[,i], levels = sort(unique(geochip[,i])))
}

# change names of columns in geochip

names(geochip)[c(1, 2, 3, 4, 5)] <- c("gene_id", "funct_trait", "funct_category",
                                   "funct_process", "type_trait")

write.csv(geochip, file = paste0(dir_res_01, "labels_geochip.csv"), row.names = F)

          tmp <- data.frame(
  Ecological_process = do.call(c, lapply(split(geochip, geochip$funct_category), function(x) {
    length(unique(x$funct_process))
  })),
  Gene_family = do.call(c, lapply(split(geochip, geochip$funct_category), function(x) {
    length(unique(x$funct_trait))
  })),
  Genes = do.call(c, lapply(split(geochip, geochip$funct_category), function(x) nrow(x))))

#%>% = then

tmp <- tmp %>% rownames_to_column(var = "Broad_ecosystem_function")

write.csv(tmp, file = paste0(dir_save, "table_distribution_of_genes_into_categories.csv"), row.names = F)

write.csv(tmp, file = paste0(dir_save, "table_distribution_of_genes_into_categories.csv"), row.names = F)

# 	MAKE A TABLE MATCHING ALL FUNCTIONAL LEVELS (gene, subcategory, category)

tmp <- geochip %>% group_by(funct_category, funct_process, funct_trait) %>% tally()

write.csv(tmp, file = paste0(dir_save, "table_match_funct_levels.csv"), row.names = F)


#==================================================================================
#PREPARE ABUNDANCE DATA FISH : RAWDATA_GEOCHIP                                    #
#==================================================================================

# Load the matrix
abundance <- read.csv(paste0(dir_data, "rawdata_geochip.csv"), row.names = 1)


# Make a presence/absence data list
abundance2 <- apply(abundance, 2, function(x) {
x[which(x!=0)] <- 1
return(x)
})

# Transformation of abundance data: Log the data 
abundancelog <- log(abundance)
# Filter data - Inf 
abundancelog[abundancelog == "-Inf"] <- 0
logabund <- round(abundancelog, 3)

# Save new csv for abundance presence/absence and abundance log 
write.csv(abundance2, file = paste0(dir_save, "table_data_abundance_presabs.cvs"), row.names = F)
write.csv(logabund, file = paste0(dir_save, "table_data_abundance_log.csv"), row.names = F)
write.csv(abundance, file = paste0(dir_res_01, "rawdata_geochip.csv"), row.names = F)



################################################################################
#   CLEAN THE METADATA                                                         #
#     - complet metadata                                                       #
#     - Filtre data                                                            #
#     - select data with ecological interest                                   #
################################################################################

# Load the metadata of each site
# Metadata <- read.csv(paste0(dir_data, "metadata.csv"))

metadata <- read.csv(paste0(dir_data, "metadata.csv"))

# Save metadata

write.csv(metadata, file = paste0(dir_save, "metadata.csv"), row.names = F)


################################################################################
#                             END OF SCRIPT                                    #
################################################################################




# Determine the variables present in several datasets

#vars <- c("Continent", "Ecosystem_type", "Site_name", "Site", "Sple_Code", "Country", 
#          "Location", "Latitude", "Longitude", "Elevation", "pH", "C_Pct",
#          "C_Total","N_Pct", "N_Total", "CN", "Moisture", "Temperature")

# Keep only these variables and make a metadata table

#tmp <- lapply(list_meta, function(x){ 
#  x <- x[, names(x) %in% vars]
#  vars_NA <- vars[! vars %in% names(x)]
#  df <- data.frame(matrix(NA, nrow(x), length(vars_NA))) %>% setNames(vars_NA)
#  cbind(x, df) %>% select(vars)
#})

#metadata <- do.call(rbind, tmp) %>% as.data.frame(row.names = 1:818) %>% 
#  dplyr::rename(Sample = Sple_Code)

# Look at the metadata distribution

#vars_fact <- c("Continent", "Country", "Location", "Ecosystem_type", "Site")
               
#vars_num <- c("Elevation", "pH", "C_Pct", "C_Total","N_Pct", "N_Total", "CN", 
#              "Moisture", "Temperature")

#metadata[, vars_fact] <- apply(metadata[, vars_fact], 2, as.factor)
#metadata[, vars_num]  <- apply(metadata[, vars_num], 2, as.numeric)

#apply(metadata[, vars_fact], 2, table)
#apply(metadata[, vars_num], 2, summary)

#out <- rbind(do.call(rbind, lapply(split(metadata, metadata$Site), function(X) {
#  apply(X, 2, function(x) round(sum(! is.na(x)) / nrow(X) * 100, 0))
#})),
#Total = apply(metadata, 2, function(x) round(sum(! is.na(x)) / 818 * 100, 0))
#)

# Save the metadata table

#write.csv(metadata, file = paste0(dir_save, "metadata.csv"))
#write.csv(out, file = paste0(dir_save, "table_summary_metadata.csv"))


# make a table of average metadata per site

#vars_fact <- c("Continent", "Country", "Location", "Ecosystem_type", "Site")

#vars_num <- c("Elevation", "pH", "C_Pct", "C_Total","N_Pct", "N_Total", "CN", 
#              "Moisture", "Temperature")

#meta_site_fact <- do.call(rbind, lapply(split(metadata, metadata$Site_name), function(X) {
#  X[1, vars_fact]
#}))

#meta_site_num  <- do.call(rbind, lapply(split(metadata, metadata$Site_name), function(X) {
#  apply(X[, vars_num], 2, mean, na.rm = TRUE)
#}))






