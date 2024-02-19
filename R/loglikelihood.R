ll_hanam <- function(theta, y, X, A, Lambda, Omega,
                   modeltype = c("HANE", "HAND")){
  n <- length(y)
  p <- ncol(X)
  D <- ncol(Lambda)

  beta    <- theta[1:p]
  gamma   <- theta[(p+1):(p+D)]
  sigmasq <- theta[p+D+1]
  rho     <- theta[p+D+2]

  M <- qr.solve(diag(1, n)-rho*A)
  ret <- ifelse(modeltype == "HANE",
                dmvnorm(y,
                        (M%*%(X%*%beta + Lambda%*%gamma))[,1],
                        as.matrix(M%*%tcrossprod(crossprod(gamma)[1]*Omega + diag(sigmasq,n), M)),
                        log = T),
                dmvnorm(y,
                        (X%*%beta + Lambda%*%gamma)[,1],
                        as.matrix(crossprod(gamma)[1]*Omega + sigmasq*tcrossprod(M)),
                        log = T))
  return(ret)
}


ll_nam <- function(theta, y, X, A, modeltype = c("NEM", "NDM")){
  n <- length(y)
  p <- ncol(X)

  beta    <- theta[1:p]
  sigmasq <- theta[p+1]
  rho     <- theta[p+2]

  M <- qr.solve(diag(1, n)-rho*A)
  ret <- ifelse(modeltype == "NEM",
                dmvnorm(y,
                        (M%*%X%*%beta)[,1],
                        as.matrix(sigmasq*tcrossprod(M)),
                        log = T),
                dmvnorm(y,
                        (X%*%beta)[,1],
                        as.matrix(sigmasq*tcrossprod(M)),
                        log = T))
  return(ret)
}

