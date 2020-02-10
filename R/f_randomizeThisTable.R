################################################################################
# This function randomizes the abundance in each communities from a gene/OTU table
#
# input: an abundance matrix with genes/OTUs in row and communities in columns
# output: a randomized abundance matrix
################################################################################

randomizeThisTable <- function(abTab){

  out <- apply(abTab, 2, function(sple){
   
    spleNo0 <- sple[sple!=0]
    
    # Randomize data for null model
      rdz <- sample(sple[sple!=0], length(sple[sple!=0]), replace= FALSE)
      names(rdz) <- names(sple[sple!=0])
      sple[sple!=0] <- rdz

	return(sple)
  })
  
colnames(out) <- colnames(abTab)
return(out)
}


