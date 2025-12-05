#' @name print
#' 
#' @title Print hanam objects.
#' 
#' @param x an object used to select a method.
#' 

#' @rdname print
#' @export
print.HANE = function(x){
  alpha = 1.0 - x$CI_level
  
  summ = x$coefs
  summ$LB = x$coefs$estimates + qnorm(0.5 * alpha) * x$coefs$se
  summ$UB = x$coefs$estimates + qnorm(1.0 - 0.5 * alpha) * x$coefs$se
  
  colnames(summ) = 
    c("Estimates", "SE", 
      paste0(round(50 * alpha,5),"%"), 
      paste0(round(100 - 50 * alpha,5),"%")
    )
  
  summ = 
    cbind(Variable = row.names(summ),
          summ)
  row.names(summ) = NULL
  
  
  summ$Variable =
    summ$Variable |> 
    gsub("beta.","beta : ",x = _) |> 
    gsub("gamma","gamma : ",x = _)
  
  print(summ)
}

#' @rdname print
#' @export
print.HAND = function(x){
  alpha = 1.0 - x$CI_level
  
  summ = x$coefs
  summ$LB = x$coefs$estimates + qnorm(0.5 * alpha) * x$coefs$se
  summ$UB = x$coefs$estimates + qnorm(1.0 - 0.5 * alpha) * x$coefs$se
  
  colnames(summ) = 
    c("Estimates", "SE", 
      paste0(round(50 * alpha,5),"%"), 
      paste0(round(100 - 50 * alpha,5),"%")
    )
  
  summ = 
    cbind(Variable = row.names(summ),
          summ)
  row.names(summ) = NULL
  
  
  summ$Variable =
    summ$Variable |> 
    gsub("beta.","beta : ",x = _) |> 
    gsub("gamma","gamma : ",x = _)
  
  print(summ)
}