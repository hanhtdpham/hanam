#' Fitting a network disturbances model
#'
#' \code{NDM} implements a normal approximation to the posterior of the network disturbances model.
#'
#' @param y a numeric vector of responses.
#' @param X a matrix of covariates.
#' @param A a row-normalized adjacency matrix.
#' @param optim.control optional control parameters for \code{\link{optim}}.
#' @param mu_b a mean vector of the normal prior on \eqn{\beta}.
#' @param s2_b a numeric value for the variance of the normal prior on \eqn{\beta}.
#' @param a_s2 two times the shape of the inverse-gamma prior on \ifelse{html}{\out{&sigma;<sup>2</sup>}}{\eqn{\sigma^2}}.
#' @param b_s2 two times the scale of the inverse-gamma prior on \ifelse{html}{\out{&sigma;<sup>2</sup>}}{\eqn{\sigma^2}}.
#' @param mu_rho a mean vector of the truncated normal prior on \eqn{\rho}.
#' @param sd_rho a numeric value for the variance of the truncated normal prior on \eqn{\rho}.
#' @param seed random seed used in \code{\link{optim}}.
#' @param init_vec an initialization vector for \eqn{(\beta, \gamma,} \ifelse{html}{\out{&sigma;<sup>2</sup>}}{\eqn{\sigma^2}}, \eqn{\rho)} in \code{\link{optim}}. If not provided, the initial values are estimated using the two-stage least squares (TSLS) estimator proposed by Kelejian & Prucha (1998).
#' @param aa a numeric value corresponds to the significance level used to report the credible intervals of the parameter estimates.
#'
#'
#' #' @details
#' \code{NDM} fits a network disturbances model given by
#' \ifelse{html}{\out{<p>y|A, X &#126; N ( X&beta;, &sigma;<sup>2</sup>(I - &rho; A)<sup>-1</sup>(I - &rho; A')<sup>-1</sup>)</p>}}{\deqn{
#'   y|A, X ~ N ( X\beta,  \sigma^2 (I - \rho A)^{-1}(I - \rho A')^{-1})
#' }}
#'
#' @import MASS
#' @import truncnorm
#' @import mvtnorm
#' @import Matrix
#'
#'
#' @export
#'
#'
#' @references
#' Doreian, P. (1980). Linear models with spatially distributed data: Spatial disturbances or spatial effects?. Sociological Methods & Research, 9(1), 29-60.
#' Kelejian, H. H., & Prucha, I. R. (1998). A generalized spatial two-stage least squares procedure for estimating a spatial autoregressive model with autoregressive disturbances. The journal of real estate finance and economics, 17, 99-121.
#'
#'
#' @seealso
#' \code{\link{optim}}
#'
#'
NDM <- function(y, X, A,
                optim.control = list(maxit = 200),
                mu_b = rep(0, ncol(X)),
                s2_b = 2.25^2,
                a_s2 = 2,
                b_s2 = 2,
                mu_rho = 0.36,
                sd_rho = 0.7,
                seed = 1,
                init_vec=NULL,
                aa = 0.05){
  # Set up
  n <- length(y)
  p <- ncol(X)

  # Check if A is a row normalized matrix
  if(! all( rowSums(A) == 1 | rowSums(A) == 0 ) ) {
    stop("Error: A is not a row normalized adjacency matrix.")
  }

  # Fit model
  set.seed(seed)
  if(is.null(init_vec)){
    fit2sls  <- lnam2sls_disturbance(y, X, A)
    init_vec <- c(fit2sls$coefs, fit2sls$s2, fit2sls$rho)
  }
  NDM_result <- optim(init_vec,
                      fn = lPostNDM, gr = gradlPostNDM,
                      method = "L-BFGS-B",
                      control= c(list(fnscale=-1), optim.control),
                      hessian = T,
                      lower = c(rep(-Inf, p), 1e-5,  -0.999),
                      upper = c(rep(Inf, p+1), 0.999),
                      y=y, X=X, A=A,
                      mu_b = mu_b, s2_b = s2_b,
                      a_s2 = a_s2, b_s2 = b_s2,
                      mu_rho = mu_rho, sd_rho =sd_rho)

  # Format return of estimates
  covar <- chol2inv(chol(nearPD(-NDM_result$hessian)$mat))
  se    <- sqrt(diag(covar))

  coefs <- data.frame(estimates = c(beta  = NDM_result$par[1:p],
                                    rho   = NDM_result$par[p+2]),
                      se        = c(beta  = se[1:p],
                                    rho   = se[p+2]))
  coefs <- data.frame(coefs,
                      LB=coefs$estimates+qnorm(aa/2)*coefs$se,
                      UB=coefs$estimates+qnorm(1-aa/2)*coefs$se)

  if(coefs["rho", "LB"] < -1 | coefs["rho", "UB"] > 1){
    coefs["rho", c("LB", "UB")] <- NA
    warning("Normal approximation to the posterior is not reliable.")
  }

  # Calculate log likelihood
  logLik <- ll_nam(NDM_result$par, y, X, A, modeltype = "NDM")

  ret <- list(init        = c(beta  = init_vec[1:p],
                              gamma = init_vec[(p+1):(p)],
                              s2    = init_vec[p+1],
                              rho   = init_vec[p+2]),
              coefs       = coefs,
              s2          = NDM_result$par[p+1],
              s2_se       = se[p+1],
              covar       = covar,
              convergence = NDM_result$convergence,
              logLik      = logLik)
  class(ret) <- "NAM"
  return(ret)
}


gradlPostNDM <- function(theta, y, X, A,
                         mu_b = rep(0, ncol(X)),
                         s2_b = 2.25^2,
                         a_s2 = 2,
                         b_s2 = 2,
                         mu_rho = 0.36,
                         sd_rho = 0.7){
  # Set up
  n <- length(y)
  p <- ncol(X)

  beta   <- theta[1:p]
  sigma  <- theta[p+1]
  rho    <- theta[p+2]

  M    <- qr.solve(diag(1,n) - rho*A)
  MMt_inv <- crossprod(diag(1,n) - rho*A)

  # Precompute
  mu       <- y - X%*%beta

  # Derivative of log posterior
  GradBeta  <- crossprod(X, MMt_inv)%*%mu/sigma - (beta-mu_b)/s2_b
  GradSigma <- -(0.5*n + 0.5*a_s2 + 1)/sigma +
    (0.5*crossprod(mu, MMt_inv)%*%mu + 0.5*b_s2)/sigma^2
  GradRho   <- -sum(diag( A%*%M )) +
    crossprod(mu, crossprod((diag(1,n) - rho*A), A))%*%mu/sigma -
    (rho - mu_rho)/(sd_rho^2)

  out <- c(GradBeta, GradSigma, GradRho)
  return(unlist(out))
}

lPostNDM <- function(theta, y, X, A,
                     mu_b = rep(0, ncol(X)),
                     s2_b = 2.25^2,
                     a_s2 = 2,
                     b_s2 = 2,
                     mu_rho = 0.36,
                     sd_rho = 0.7){
  n <- length(y)
  p <- ncol(X)

  beta     <- theta[1:p]
  sigmasq  <- theta[p+1]
  rho      <- theta[p+2]

  M        <- qr.solve(diag(1, n)-rho*A)

  dmvnorm(y,
          (X%*%beta)[,1],
          as.matrix(sigmasq*tcrossprod(M)),
          log = T)+
    dmvnorm(beta , mu_b, diag(s2_b, p), log = T) -
    (0.5*a_s2+1)*log(sigmasq) - 0.5*b_s2/sigmasq +
    log(dtruncnorm(rho, -1, 1, mu_rho, sd_rho))
}

