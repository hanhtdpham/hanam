#' @name summary
#' 
#' @title Summary functions for hanam objects
#' 
#' @param object HANE or HAND object
#' @param CI_level Posterior probability covered by credible interval
#' 
#' 
#' @export

#' @rdname summary
#' @export
summary.HANE = function(object,
                        CI_level = 0.95){
  alpha = 1 - CI_level
  summ = object$coefs
  summ$LB = object$coefs$estimates + qnorm(0.5 * alpha) * object$coefs$se
  summ$UB = object$coefs$estimates + qnorm(1.0 - 0.5 * alpha) * object$coefs$se
  
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
      gsub("gamma","gamma : U_",x = _)
  
  
  
  summ_for_printing = summ
  summ_for_printing[,-1] = 
    format(signif(summ_for_printing[,-1], 3), 
           scientific = FALSE)
  cat("\n----------\n\nHomophily-adjusted Network Effects (HANE) Model\n\n")
  print(summ_for_printing)
  cat("\n----------\n\n")
  
  invisible(summ)
}

#' @rdname summary
#' @export

summary.HAND = function(object,
                        CI_level = 0.95){
  alpha = 1 - CI_level
  summ = object$coefs
  summ$LB = object$coefs$estimates + qnorm(0.5 * alpha) * object$coefs$se
  summ$UB = object$coefs$estimates + qnorm(1.0 - 0.5 * alpha) * object$coefs$se
  
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
    gsub("gamma","gamma : U_",x = _)
  
  
  
  summ_for_printing = summ
  summ_for_printing[,-1] = 
    format(signif(summ_for_printing[,-1], 3), 
           scientific = FALSE)
  cat("\n----------\n\nHomophily-adjusted Network Disturbances (HAND) Model\n\n")
  print(summ_for_printing)
  cat("\n----------\n\n")
  
  invisible(summ)
}