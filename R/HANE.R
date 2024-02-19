#' Fitting a homophily adjusted network effects model
#'
#' \code{HANE} implements a normal approximation to the posterior of the homophily
#' adjusted network effects (HANE) model.
#'
#' @param y a numeric vector of responses.
#' @param X a matrix of covariates.
#' @param A a row-normalized adjacency matrix.
#' @param Usample a sample of latent locations to use as priors; an array of size (size, n, D), where size is the number of sample, n is the number of nodes in the network, and D is the number of latent dimension.
#' @param optim.control optional control parameters for \code{\link{optim}}.
#' @param mu_b a mean vector of the normal prior on \eqn{\beta}.
#' @param s2_b a numeric value for the variance of the normal prior on \eqn{\beta}.
#' @param mu_g a mean vector of the normal prior on \eqn{\gamma}.
#' @param s2_g a numeric value for the variance of the normal prior on \eqn{\gamma}.
#' @param a_s2 two times the shape of the inverse-gamma prior on \ifelse{html}{\out{&sigma;<sup>2</sup>}}{\eqn{\sigma^2}}.
#' @param b_s2 two times the scale of the inverse-gamma prior on \ifelse{html}{\out{&sigma;<sup>2</sup>}}{\eqn{\sigma^2}}.
#' @param mu_rho a mean vector of the truncated normal prior on \eqn{\rho}.
#' @param sd_rho a numeric value for the variance of the truncated normal prior on \eqn{\rho}.
#' @param seed random seed used in \code{\link{optim}}.
#' @param init_vec an initialization vector for \eqn{(\beta, \gamma,} \ifelse{html}{\out{&sigma;<sup>2</sup>}}{\eqn{\sigma^2}}, \eqn{\rho)} in \code{\link{optim}}. If not provided, the initial values are estimated using the two-stage least squares (TSLS) estimator proposed by Kelejian & Prucha (1998).
#' @param aa a numeric value corresponds to the significance level used to report the credible intervals of the parameter estimates.
#' @param verbose a logical value indicating whether to print steps of the algorithm.
#' @param Usample.maxiter a integer value indicating the maximum iteration used in approximating Usample with a matrix normal distribution.
#' @param Usample.eps a numeric value indicating the precision used in approximating Usample with a matrix normal distribution.
#'
#'
#' @details
#' \code{HANE} fits a homophily adjsuted network effects model given by
#' \ifelse{html}{\out{<p>y|A, X &#126; N ( (I - &rho; A)<sup>-1</sup>(&Lambda;&gamma; + X&beta;), (I - &rho; A)<sup>-1</sup> (&gamma;'&gamma; &Omega; + &sigma;<sup>2</sup> I) (I - &rho; A')<sup>-1</sup>)</p>}}{\deqn{
#'   y|A, X ~ N ( (I - \rho A)^{-1}(\Lambda\gamma + X\beta),  (I - \rho A)^{-1} (\gamma'\gamma \Omega + \sigma^2 I)(I - \rho A')^{-1})
#' }}
#' where \eqn{\Lambda} and \eqn{\Omega} are MLE of the sample of U assuming a matrix normal distribution with mean \eqn{\Lambda} and row variance covariance \eqn{\Omega}.
#'
#'
#' @import MASS
#' @import truncnorm
#' @import mvtnorm
#' @import Matrix
#' @export
#'
#'
#' @references
#' Kelejian, H. H., & Prucha, I. R. (1998). A generalized spatial two-stage least squares procedure for estimating a spatial autoregressive model with autoregressive disturbances. The journal of real estate finance and economics, 17, 99-121.
#'
#'
#' @seealso
#' \code{\link{optim}}
#'
#'
HANE <- function(y, X, A, Usample,
                 optim.control = list(maxit = 200),
                 mu_b = rep(0, ncol(X)),
                 s2_b = 2.25^2,
                 mu_g = rep(0, dim(Usample)[3]),
                 s2_g = 2.25^2,
                 a_s2 = 2,
                 b_s2 = 2,
                 mu_rho = 0.36,
                 sd_rho = 0.7,
                 seed = 1,
                 init_vec=NULL,
                 aa = 0.05,
                 verbose = F,
                 Usample.maxiter = 100,
                 Usample.eps = 1e-8){
  # Set up
  n <- length(y)
  p <- ncol(X)
  nUsamp <- dim(Usample)[1]
  D <- dim(Usample)[3]

  # Check if A is a row normalized matrix
  if(! all( rowSums(A) == 1 | rowSums(A) == 0 ) ) {
    stop("Error: A is not a row normalized adjacency matrix.")
  }

  # Check if U is an array
  if(length(dim(Usample)) != 3 | dim(Usample)[2] != n){
    stop("Error: Input of Usample not in the right array format (size, n, D)")
  }

  # Calculate matrix normal approximation to prior of U
  if(verbose) cat('Deriving matrix normal approximation to Usample\n')
  Uapprox <- Uprior_appox(Usample, Usample.maxiter, Usample.eps)
  if(!Uapprox$conv) stop('Approximation to Usample failed! Consider increasing Usample.maxiter or decreasing Usample.eps')
  Lambda  <- Uapprox$Lambda
  Omega   <- Uapprox$Omega
  Psi     <- Uapprox$Psi

  # Fit model
  if(verbose) cat('Fitting HANE model\n')
  set.seed(seed)
  # Initialize
  if(is.null(init_vec)){
    fit2sls  <- lnam2sls_effect(y, cbind(X, Lambda), A)
    init_vec <- c(fit2sls$coefs[1:(p+D)], fit2sls$s2, fit2sls$coefs[p+D+1])
  }
  HANE_result <- optim(init_vec,
                       fn = lPostHANE, gr = gradlPostHANE,
                       method = "L-BFGS-B",
                       control= c(list(fnscale=-1), optim.control),
                       hessian = T,
                       lower = c(rep(-Inf, p+D), 1e-5,  -0.999),
                       upper = c(rep(Inf, p+D+1), 0.999),
                       y=y, X=X, A=A, Lambda=Lambda, Omega=Omega, Psi=Psi,
                       mu_b = mu_b, s2_b = s2_b,
                       mu_g = mu_g, s2_g = s2_g,
                       a_s2 = a_s2, b_s2 = b_s2,
                       mu_rho = mu_rho, sd_rho =sd_rho)

  # Format return of estimates
  covar <- chol2inv(chol(nearPD(-HANE_result$hessian)$mat))
  se    <- sqrt(diag(covar))

  coefs <- data.frame(estimates = c(beta  = HANE_result$par[1:p],
                                    gamma = HANE_result$par[(p+1):(p+D)],
                                    rho   = HANE_result$par[p+D+2]),
                      se        = c(beta  = se[1:p],
                                    gamma = se[(p+1):(p+D)],
                                    rho   = se[p+D+2]))
  coefs <- data.frame(coefs,
                      LB=coefs$estimates+qnorm(aa/2)*coefs$se,
                      UB=coefs$estimates+qnorm(1-aa/2)*coefs$se)

  if(coefs["rho", "LB"] < -1 | coefs["rho", "UB"] > 1){
    coefs["rho", ] <- NA
    warning("Normal approximation to the posterior is not reliable.")
  }

  # Calculate log likelihood
  logLik <- ll_hanam(HANE_result$par, y, X, A,
                     Lambda, Omega, modeltype = "HANE")

  ret <- list(init        = c(beta  = init_vec[1:p],
                              gamma = init_vec[(p+1):(p+D)],
                              s2    = init_vec[p+D+1],
                              rho   = init_vec[p+D+2]),
              Lambda      = Lambda,
              Omega       = Omega,
              Psi         = Psi,
              coefs       = coefs,
              s2          = HANE_result$par[p+D+1],
              s2_se       = se[p+D+1],
              covar       = covar,
              convergence = HANE_result$convergence,
              logLik      = logLik)
  class(ret) <- "HANE"
  return(ret)
}


gradlPostHANE <- function(theta, y, X, A,
                          Lambda, Omega, Psi,
                          mu_b = rep(0, ncol(X)),
                          s2_b = 2.25^2,
                          mu_g = rep(0, ncol(Lambda)),
                          s2_g = 2.25^2,
                          a_s2 = 2,
                          b_s2 = 2,
                          mu_rho = 0.36,
                          sd_rho = 0.7){
  # Set up
  n <- length(y)
  p <- ncol(X)
  D <- ncol(Lambda)

  beta    <- theta[1:p]
  gamma   <- theta[(p+1):(p+D)]
  sigmasq <- theta[p+D+1]
  rho     <- theta[p+D+2]

  M    <- qr.solve(diag(1, n)-rho*A)
  MMt  <- tcrossprod(M)
  Sinv <- chol2inv(chol( crossprod(gamma, Psi%*%gamma)[1]*Omega + diag(sigmasq, n) ))

  # Precompute
  mu     <- ((diag(1, n)-rho*A)%*%y - Lambda%*%gamma - X%*%beta)
  muSinv <- crossprod(mu, Sinv)

  # Derivative of sigmasq
  Sgamma <- lapply(1:D, function(x) 2*(Psi%*%gamma)[x,1]*Omega)

  # Derivative of log posterior
  GradBeta  <- crossprod(X, Sinv)%*%mu - (beta-mu_b)/s2_b
  GradGamma <- sapply(1:D, function(g){
    ( crossprod(Lambda[,g], Sinv)%*%mu - (gamma[g]-mu_g[g])/s2_g -
        0.5*sum(diag( Sinv%*%Sgamma[[g]] )) + 0.5*muSinv%*%tcrossprod(Sgamma[[g]], muSinv) )[1]
  })
  Gradsigmasq <- -0.5*sum(diag( Sinv )) + 0.5*tcrossprod(muSinv) -
    (0.5*a_s2 + 1)/sigmasq + 0.5*b_s2/sigmasq^2
  GradRho   <- -sum(diag(A%*%M)) +
    crossprod(A%*%y, Sinv%*%mu) - (rho - mu_rho)/(sd_rho^2)

  out <- c(GradBeta, GradGamma, Gradsigmasq, GradRho)
  return(unlist(out))
}

lPostHANE <- function(theta, y, X, A,
                      Lambda, Omega, Psi,
                      mu_b = rep(0, ncol(X)),
                      s2_b = 2.25^2,
                      mu_g = rep(0, ncol(Lambda)),
                      s2_g = 2.25^2,
                      a_s2 = 2,
                      b_s2 = 2,
                      mu_rho = 0.36,
                      sd_rho = 0.7){
  n <- length(y)
  p <- ncol(X)
  D <- ncol(Lambda)

  beta    <- theta[1:p]
  gamma   <- theta[(p+1):(p+D)]
  sigmasq <- theta[p+D+1]
  rho     <- theta[p+D+2]

  M    <- qr.solve(diag(1, n)-rho*A)

  dmvnorm(y,
          (M%*%(X%*%beta + Lambda%*%gamma))[,1],
          as.matrix(M%*%tcrossprod(crossprod(gamma, Psi%*%gamma)[1]*Omega + diag(sigmasq,n), M)),
          log = T)+
    dmvnorm(beta , mu_b, diag(s2_b, p), log = T) +
    dmvnorm(gamma, mu_g, diag(s2_g, D), log = T) -
    (0.5*a_s2+1)*log(sigmasq) - 0.5*b_s2/sigmasq +
    log(dtruncnorm(rho, -1, 1, mu_rho, sd_rho))
}


lnam2sls_effect = function(y,X,A){
  ZZ = cbind(X,A%*%y)
  if(mean(X[,1] == 1) == 1){
    X1 = X[,-1]
  }else{
    X1 = X
  }
  AX = A%*%X1
  HH = model.matrix(~X1+AX+A%*%AX)
  PP = tcrossprod(HH%*%chol2inv(chol(crossprod(HH))),HH)
  ZHat = PP%*%ZZ
  ZtZInv= qr.solve(crossprod(ZHat))
  Coefs = drop(ZtZInv%*%crossprod(ZHat,y))

  s2Hat = drop(crossprod(y-ZZ%*%Coefs))/nrow(X)

  return(list(coefs=Coefs,s2=s2Hat))
}
