#
# Function to determine the type of model and it's shape


determine_model_shape <- function(x, x_pred, pval_thd = 0.01, intercept_thd = 1) { 
  
    # Is it a linear model ?
    if (x$best_model == "mod1") {
      
      # Is there a significant rank effect ?
      if (x$mod1_aov_pval_Rank > pval_thd) { # no
        shape <- c("No relation", "No relation")
      } else { # yes
        if (x$mod1_Rank > 0) { # it is positive
          shape <- c("Linear", "Positive")
        } 
        if (x$mod1_Rank < 0) { # it is negative
          shape <- c("Linear", "Negative") }
      }
    }
    
    # Is it an exponential model ?
    if (x$best_model == "mod3") {
      
      # Is there a significant rank effect ?
      if (x$mod3_aov_pval_Rank > pval_thd) { # no
        shape <- c("No relation", "No relation")
      } else { # yes
        if (x$mod3_Rank > 0) { # it is positive
          shape <- c("Exponential", "Positive")
        } 
        if (x$mod3_Rank < 0) { # it is negative
          shape <- c("Exponential", "Negative") }
      }
    }   
    
    # Is it a quadratic model ?
    if (x$best_model == "mod2") {
      
      # Are both Rank effects non-significant ?
      if (x$mod2_aov_pval_Rank > pval_thd && x$mod2_aov_pval_Ranksq > pval_thd) {
        shape <- c("No relation", "No relation")
      }
      
      # is there at least one significant rank effect ?
      if (sum(c(x$mod2_aov_pval_Rank < pval_thd, 
                x$mod2_aov_pval_Ranksq < pval_thd)) != 0) {
        
        # Is the relationhip convex ?
        
        if (x$mod2_Rank < 0 && x$mod2_Ranksq > 0) {
            
          # is the intercept positive ?
          if (x$mod2_Intercept > intercept_thd) {
            
            if (x_pred$mod2[1] > intercept_thd && x_pred$mod2[10] > intercept_thd) {
              shape <- c("Quadratic", "U shaped")
            }
            # Is the relation negative ?
            if (x_pred$mod2[1] >= intercept_thd && x_pred$mod2[10] < intercept_thd) {
              shape <- c("Quadratic", "Negative")
            }
            # Is the relation positive ?
            if (x_pred$mod2[1] <= intercept_thd && x_pred$mod2[10] > intercept_thd) {
              shape <- c("Quadratic", "Positive")
            }
            
          }
          
          # is the intercept negative ?
          if (x$mod2_Intercept < intercept_thd) {
            
            if (x_pred$mod2[1] > intercept_thd && x_pred$mod2[10] > intercept_thd) {
              shape <- c("Quadratic", "U shaped")
            }
            # Is the relation negative ?
            if (x_pred$mod2[1] >= intercept_thd && x_pred$mod2[10] < intercept_thd) {
              shape <- c("Quadratic", "Negative")
            }
            # Is the relation positive ?
            if (x_pred$mod2[1] <= intercept_thd && x_pred$mod2[10] > intercept_thd) {
              shape <- c("Quadratic", "Positive")
            }
            
          }
          
          
        }
        
        # Is the relationhip concav ?
        if (x$mod2_Rank > 0 && x$mod2_Ranksq < 0) {
          
          # is the intercept negative ?
          if (x$mod2_Intercept < intercept_thd) {
            
            # Is the relation bell shapêd ?
            if (x_pred$mod2[1] < intercept_thd && x_pred$mod2[10] < intercept_thd) {
              shape <- c("Quadratic", "Bell shaped")
            }
            # Is the relation positive ?
            if (x_pred$mod2[1] <= intercept_thd && x_pred$mod2[10] > intercept_thd) {
              shape <- c("Quadratic", "Positive")
            }
            # Is the relation negative ?
            if (x_pred$mod2[1] >= intercept_thd && x_pred$mod2[10] < intercept_thd) {
              shape <- c("Quadratic", "Negative")
            }
            
          }
          
          # is the intercept positive ?
          if (x$mod2_Intercept > intercept_thd) {

            # Is the relation bell shapêd ?
            if (x_pred$mod2[1] < intercept_thd && x_pred$mod2[10] < intercept_thd) {
              shape <- c("Quadratic", "Bell shaped")
            }
            # Is the relation positive ?
            if (x_pred$mod2[1] <= intercept_thd && x_pred$mod2[10] > intercept_thd) {
              shape <- c("Quadratic", "Positive")
            }
            # Is the relation negative ?
            if (x_pred$mod2[1] >= intercept_thd && x_pred$mod2[10] < intercept_thd) {
              shape <- c("Quadratic", "Negative")
            }

          }
        }
          
        if (x$mod2_Rank > 0 && x$mod2_Ranksq > 0) {
          
          # Is the relation positive or negative
          if (x$mod2_Intercept < intercept_thd && 
              x_pred$mod2[1] <= intercept_thd && 
              x_pred$mod2[10] > intercept_thd) {
            shape <- c("Quadratic", "Positive")
          }
          if (x$mod2_Intercept > intercept_thd && 
              x_pred$mod2[1] >= intercept_thd && 
              x_pred$mod2[10] < intercept_thd) {
            shape <- c("Quadratic", "Negative")
          }
          if (x$mod2_Intercept > intercept_thd && 
              x_pred$mod2[1] <= intercept_thd && 
              x_pred$mod2[10] > intercept_thd) {
            shape <- c("Quadratic", "Positive")
          }
          
        }
        
        if (x$mod2_Rank < 0 && x$mod2_Ranksq < 0) {
          
          # Is the relation positive or negative
          if (x$mod2_Intercept < intercept_thd && 
              x_pred$mod2[1] <= intercept_thd && 
              x_pred$mod2[10] > intercept_thd) {
            shape <- c("Quadratic", "Positive")
          }
          if (x$mod2_Intercept > intercept_thd && 
              x_pred$mod2[1] >= intercept_thd && 
              x_pred$mod2[10] < intercept_thd) {
            shape <- c("Quadratic", "Negative")
          }
          if (x$mod2_Intercept > intercept_thd && 
              x_pred$mod2[1] <= intercept_thd && 
              x_pred$mod2[10] > intercept_thd) {
            shape <- c("Quadratic", "Positive")
          }
          
        }
        
     }
    }

  return(shape)
}


