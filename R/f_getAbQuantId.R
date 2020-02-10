################################################################################
# This function associates each gene/OTU from a community to an abundance quantile
#
# input: 
#  an abundance matrix with genes/OTUs in row and communities in columns
# output: 
#  a matrix with genes/OTUs in row and communities in columns in which
#  the values correspond to the quantile (from 1 to 10)
#
# - method = "trueQuant" or NULL
#  define the method used to determine the abundance cutoffs used to separate 
#   genes/OTUs into quantiles. trueQuant uses the quantiles() function while
#   any other value generates a sequence of length 11 between the minimal and
#   maximal values observed in the community
#
# - nullModel = TRUE ro FALSE
#  determine whether abundance are randomized prior to quantile identification
################################################################################


getAbQuantId <- function(abTab, method= "trueQuant", nullModel=FALSE,
                         num_of_quantiles = 10){

  out <- apply(abTab, 2, function(sple){
   
    spleNo0 <- sple[sple!=0]
    
    # Randomize data for null model
    if(nullModel == TRUE){
      rdz <- sample(spleNo0, length(sple[sple!=0]), replace= FALSE)
      names(rdz) <- names(spleNo0)
      sple[sple!=0] <- rdz
    }
    
	  # Define the quantiles of abundance
	  if( method =='trueQuant'){
	    quantilesCuts <- stats::quantile(spleNo0, probs = seq(0, 1, 
	                                                   by = 1 / num_of_quantiles))
	  }else{
	    quantilesCuts <- seq(min(spleNo0), max(spleNo0), 
	                         length.out = num_of_quantiles + 1)
	  }
	  
    # Make a probeXquantile matrix
    tplte <- rep(0, length(sple))
    names(tplte) <- 1:length(sple)
    
    for(i in 1:(length(quantilesCuts)-1))  {
      lims <- c(quantilesCuts[i], quantilesCuts[i+1]) 
      
      mask <- which(lims[1] <= sple & sple < lims[2] )
      if(i == 10){
        mask <- which(lims[1] <= sple & sple <= lims[2] )
      }
      tplte[mask] <- i 
    }

	return(tplte)
  })
  
colnames(out) <- colnames(abTab)
return(out)
}




