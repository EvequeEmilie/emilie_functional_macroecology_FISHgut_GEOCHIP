################################################################################
# This function associates each gene/OTU from a community to an abundance quantile
#
# input: 
#  - quantile_matrix : a matrix with genes/OTUs in row and communities 
#         in columns as generated using the f_getAbQuantId.R function
#  - originalTable : the original gene table, used toweight composition by gene abundance
#  - CLM : the common label matrix describing the link between genes and functions
#
# output: 
# a list of the length of the number of communities in the original table
#  each element of the list contains the functions in row and the
#  abundance quantiles in column
#
# - funcLevel = "funct_category"
#  define level of functional resoluation used to aggregate gene abundance
#
# - abWeight = TRUE or FALSE
#  determine whether abundance of functions are weighted by the abundance 
#   of the genes that carry them
#
# - scaled = TRUE or FALSE
#  determine whether the observed function abundance is scaled y the expected
#   abundance considering the number of probes from this function on the microarray
################################################################################



getCompoFromQuantilTable <- function(quantile_matrix, originalTable, CLM,
                                     funcLevel = "funct_category", 
                                     abWeight = TRUE, scaled = TRUE) {

  numQ <- max(quantile_matrix)
  compoList <- list() 
  
  for (i in colnames(quantile_matrix)) {

    vec <- quantile_matrix[, i]
    names(vec) <- CLM[, funcLevel]
    
    # Get the composition
    compoMat <- do.call(cbind, lapply(1:numQ, function(q) {
      table(CLM[which(vec == q), funcLevel])
    }))
    colnames(compoMat) <- paste0("Q", 1:numQ)
    
    # Get the composition weighted by genes abundance
    if (abWeight) {
      compoMat <- do.call(cbind, lapply(1:numQ, function(q) {
        fab <- do.call(c, lapply(levels(CLM[, funcLevel]), function(f) {
          sum(originalTable[which(vec == q & names(vec) == f), i])
        }))
        names(fab) <- levels(CLM[, funcLevel])
        fab
      }))
      colnames(compoMat) <- paste0("Q", 1:numQ)
    }

    # Scale abundance of functions using expected abundance
    if (scaled) {
      scal <- base::table(CLM[, funcLevel]) / nrow(CLM) 
      compoMat <- apply(compoMat, 2, function(x) (x / sum(x)) / scal)
      compoMat[is.na(compoMat)] <- 0
    }
    compoList[[i]] <- compoMat
  }
return(compoList)
}


