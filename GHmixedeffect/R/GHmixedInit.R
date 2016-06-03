#' measurement error input
#' TODO: generlize the class so that one put in only theta
#'       and then theta arguments should be passed on depending on noise class
#'       
#' @param Y     - list with the observations
#' @param noise - the aviable noise classes
#' @param sigma - measurement variance
#' @param nu    - shape parameter for noise (NIG only)
MeasurementErrorInit <- function(Y,
                                 noise = c("Normal","NIG"),
                                 sigma = 1.,
                                 nu    = 1.){
  
  noise <- match.arg(noise)
  output <- list(noise = noise)
  output$sigma = sigma
  
  if(noise == "NIG"){
    output$nu = nu
    output$Vin <- list()
    for(i in 1:length(Y))
    {
      output$Vin[[i]] <- rep(1, length(Y[[i]]))
    }
  }
  return(output)
}


#' Mixed effect input
#' TODO: generlize the class so that one put in only theta
#'       and then theta arguments should be passed on depending on noise class
#'       
#' @param noise - the aviable noise classes
#' @param B_fixed     - fixed effect
#' @param B_random    - random effect
MixedInit <- function(noise = c("Normal","NIG"),
                      B_fixed  = NULL,
                      B_random = NULL){
  
  noise <- match.arg(noise)
  output <- list(noise = noise)
  if(is.null(B_random) == 0)
    output$B_random = B_random 
  
  if(is.null(B_fixed) == 0)
    output$B_fixed = B_fixed
  
  return(output)
}