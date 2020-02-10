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

# Load raw data : GEOCHIHP
geochip = read.csv(paste0 ("data/", "labels_geochip.csv"), row.names = 1)
# Load table_match_levels_data of script 01
tmp <- read.csv(fille = paste0(dir_res_01, "table_match_funct_levels.csv"), row.names = 1)

################################################################################


################################################################################
#                           EXEMPLE LOOP                                       #
#                                                                              #
################################################################################
# 1 at 10

ls_toto=list()

for(i in 1:3){
  
  x=c(10,15,20)[i]
  tab <- filter(tmp, n > x)
  dim(tab)
  ls_toto[[i]]=table(tab$funct_category        )
}

#
for (i in c( > 10)) print(i)
x = ifelse(x>10, , )


################################################################################
#                                                                              #
#                         FILTER DATA LESS PROBES                              #
#                                                                              #
################################################################################
# Removal of genes with less than 10 probes 
filter10 <- filter(tmp, n > 10)
grp10 <- group_by(filter10, funct_category, n) %>% summarise(calcul_mean = mean(n)) %>% count ()

#-------------------------------------------------------------------------------

# Removal of genes with less than 15 probes 
filter15 <- filter(tmp, n > 15)

#-------------------------------------------------------------------------------

# Removal of genes with less than 20 probes 
filter20 <- filter(tmp, n > 20)

#-------------------------------------------------------------------------------

# Removal of genes with less than 25 probes 
filter25 <- filter(tmp, n > 25)

#-------------------------------------------------------------------------------

# Removal of genes with less than 30 probes 
filter30 <- filter(tmp, n > 30)

#-------------------------------------------------------------------------------

# Removal of genes with less than 35 probes 
filter35 <- filter(tmp, n > 35)

#-------------------------------------------------------------------------------

# Removal of genes with less than 40 probes 
filter40 <- filter(tmp, n > 40)

###################################################################################




###################################################################################
#                                                                                 #
#                             RESUM_TABLE_FILTER                                  #
#                                                                                 #
###################################################################################
# Create a base-matrix
matrice <- matrix(data = 1:32, nrow=4, ncol=8, byrow = T)
matrice
# Create a vector
levels.category <- c(12,12,12,12,12,12,12,12)
levels.porcess <- c(121,101,95,94,90,88,84,81) 
levels.trait <- c(420,283,260,253,239,229,213,202)
sum.n <- c(67628,67019,66730,66611,66292,66019,65498,65080)

# Group all data in the same vector 
resultat_matrice_filtre <- matrix(c(levels.category, levels.porcess, levels.trait, sum.n), nrow = 4, ncol =8, byrow = T)

# Names rows and columns
rownames(resultat_matrice_filtre) <- c("levels.category", "levels.process", "levels.trait", "sum.n")
colnames(resultat_matrice_filtre) <- c("f0","f10","f15","f20","f25","f30","f35","f40")
# Show the result
resultat_matrice_filtre

####################################################################################


###################################################################################
#                                                                                 #
#                             TABLE_funct_category_FILTER                         #
#                                                                                 #
###################################################################################
# Count funct_category
table_funct_category <- count(tmp$funct_category)
# Count funct_category for the filter 10
table_funct_category10 <- count(filter10$funct_category)
# Count funct_category for the filter 15
table_funct_category15 <- count(filter15$funct_category)
# Count funct_category for the filter 20
table_funct_category20 <- count(filter20$funct_category)
# Count funct_category for the filter 25
table_funct_category25 <- count(filter25$funct_category)
# Count funct_category for the filter 30
table_funct_category30 <- count(filter30$funct_category)
# Count funct_category for the filter 35
table_funct_category35 <- count(filter35$funct_category)
# Count funct_category for the filter 40
table_funct_category40 <- count(filter40$funct_category)

# Join all table 
total_table <- bind_cols(table_funct_category,table_funct_category10, table_funct_category15, table_funct_category20, table_funct_category25, table_funct_category30, table_funct_category35, 
          table_funct_category40)

# Removal the False data 
Table_funct_category <- total_table %>% group_by(x, freq, freq1, freq2, freq3, freq4, freq5, freq6, freq7) %>% tally()

names(geochip)[c(1, 2, 3, 4, 5, 6, 7)] <- c("category", "filter0", "filter10",
                                      "filter15", "filter20", "filter25", "filter30")

