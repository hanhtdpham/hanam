#' Fitting a Homophily Adjusted Network Disturbances Model
#'
#' \code{HAND} implements a normal approximation to the posterior of the homophily
#' adjusted network disturbances (HAND) model.
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
#' \code{HAND} fits a homophily adjsuted network disturbances model given by
#' \ifelse{html}{\out{<p>y|A, X &#126; N (&Lambda;&gamma; + X&beta;,  &gamma;'&gamma; &Omega; + &sigma;<sup>2</sup> (I - &rho; A)<sup>-1</sup>(I - &rho; A')<sup>-1</sup>)</p>}}{\deqn{
#'   y|A, X ~ N (\Lambda\gamma + X\beta,  \gamma'\gamma \Omega + \sigma^2 (I - \rho A)^{-1}(I - \rho A')^{-1})
#' }}
#' where \eqn{\Lambda} and \eqn{\Omega} are MLE of the sample of U assuming a matrix normal distribution with mean \eqn{\Lambda} and row variance covariance \eqn{\Omega}.
#'
#'
#' @import MASS
#' @import truncnorm
#' @import mvtnorm
#' @import Matrix
#' @import numDeriv
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
HAND <- function(y, X, A, Usample,
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
    stop('Error: A is not a row normalized adjacency matrix.')
  }

  # Check if U is an array
  if(length(dim(Usample)) != 3 | dim(Usample)[2] != n){
    stop('Error: Input of Usample not in the right array format (size, n, D)')
  }

  # Calculate matrix normal approximation to prior of U
  if(verbose) cat('Deriving matrix normal approximation to Usample\n')
  Uapprox <- Uprior_appox(Usample, Usample.maxiter, Usample.eps)
  if(!Uapprox$conv) stop('Approximation to Usample failed! Consider increasing Usample.maxiter or decreasing Usample.eps')
  Lambda  <- Uapprox$Lambda
  Omega   <- Uapprox$Omega
  Psi     <- Uapprox$Psi

  # Fit model
  if(verbose) cat('Fitting HAND model\n')
  set.seed(seed)
  if(is.null(init_vec)){
    fit2sls  <- lnam2sls_disturbance(y, cbind(X, Lambda), A)
    init_vec <- c(fit2sls$coefs, fit2sls$s2, fit2sls$rho)
  }
  HAND_result <- optim(init_vec,
                       fn = lPostHAND, gr = gradlPostHAND,
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
  covar <- chol2inv(chol(nearPD(-HAND_result$hessian)$mat))
  se    <- sqrt(diag(covar))

  coefs <- data.frame(estimates = c(beta  = HAND_result$par[1:p],
                                    gamma = HAND_result$par[(p+1):(p+D)],
                                    rho   = HAND_result$par[p+D+2]),
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
  logLik <- ll_hanam(HAND_result$par, y, X, A,
                     Lambda, Omega, modeltype = "HAND")

  ret <- list(init        = c(beta  = init_vec[1:p],
                              gamma = init_vec[(p+1):(p+D)],
                              s2    = init_vec[p+D+1],
                              rho   = init_vec[p+D+2]),
              Lambda      = Lambda,
              Omega       = Omega,
              Psi         = Psi,
              coefs       = coefs,
              s2          = HAND_result$par[p+D+1],
              s2_se       = se[p+D+1],
              covar       = covar,
              convergence = HAND_result$convergence,
              logLik      = logLik)
  class(ret) <- "HAND"
  return(ret)
}


gradlPostHAND <- function(theta, y, X, A,
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

  beta   <- theta[1:p]
  gamma  <- theta[(p+1):(p+D)]
  sigma  <- theta[p+D+1]
  rho    <- theta[p+D+2]

  M    <- qr.solve(diag(1,n) - rho*A)
  Mt   <- t(M)
  MMt  <- tcrossprod(M)
  Sinv <- chol2inv(chol( crossprod(gamma, Psi%*%gamma)[1]*Omega + sigma*MMt ))

  # Precompute
  mu       <- y - Lambda%*%gamma - X%*%beta
  muSinv   <- crossprod(mu, Sinv)

  # Derivate of Sigma
  Srho   <- sigma*(M%*%A%*%MMt + MMt%*%crossprod(A, Mt))
  Sgamma <- lapply(1:D, function(x) 2*(Psi%*%gamma)[x,1]*Omega)

  # Derivative of log posterior
  GradBeta  <- crossprod(X, Sinv)%*%mu - (beta-mu_b)/s2_b
  GradGamma <- sapply(1:D, function(g){
    ( crossprod(Lambda[,g], Sinv)%*%mu - (gamma[g]-mu_g[g])/s2_g -
        0.5*sum(diag( Sinv%*%Sgamma[[g]] )) + 0.5*muSinv%*%tcrossprod(Sgamma[[g]], muSinv) )[1]
  })
  GradSigma <- -0.5*sum(diag( Sinv%*%MMt )) + 0.5*muSinv%*%tcrossprod(MMt, muSinv) -
    (0.5*a_s2 + 1)/sigma + 0.5*b_s2/sigma^2
  GradRho   <- -0.5*sum(diag( Sinv%*%Srho )) + 0.5*muSinv%*%tcrossprod(Srho, muSinv) -
    (rho - mu_rho)/(sd_rho^2)

  out <- c(GradBeta, GradGamma, GradSigma, GradRho)
  return(unlist(out))
}

lPostHAND <- function(theta, y, X, A,
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

  beta   <- theta[1:p]
  gamma  <- theta[(p+1):(p+D)]
  sigmasq  <- theta[p+D+1]
  rho    <- theta[p+D+2]

  M    <- qr.solve(diag(1, n)-rho*A)

  dmvnorm(y,
          (X%*%beta + Lambda%*%gamma)[,1],
          as.matrix(crossprod(gamma, Psi%*%gamma)[1]*Omega +
                      M%*%tcrossprod(diag(sigmasq,n), M)),
          log = T)+
    dmvnorm(beta , mu_b, diag(s2_b, p), log = T) +
    dmvnorm(gamma, mu_g, diag(s2_g, D), log = T) -
    (0.5*a_s2+1)*log(sigmasq) - 0.5*b_s2/sigmasq +
    log(dtruncnorm(rho, -1, 1, mu_rho, sd_rho))
}


lnam2sls_disturbance = function(y,X,A,optim.init=NULL, seed=NULL){
  ZZ = X
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
  DeltaTilde = ZtZInv%*%crossprod(ZHat,y)

  u       = drop(y - ZZ%*%DeltaTilde)
  uBar    = drop(A%*%u)
  uBarBar = drop(A%*%A%*%u)
  G = matrix(c(
    2*mean(u*uBar)          ,       -mean(uBar^2),                 1/nrow(X),
    2*mean(uBar*uBarBar)    ,    -mean(uBarBar^2),  mean(diag(crossprod(A))),
    mean(u*uBarBar + uBar^2), -mean(uBar*uBarBar),                         0
  ), nrow = 3, byrow = T)
  g = c(mean(u^2), mean(uBar^2), mean(u*uBar))
  objFun <- function(x){
    alpha <- c(x[1], x[1], x[2])
    drop(crossprod(g - G%*%alpha))
  }
  if(is.null(optim.init)){
    if(!is.null(seed)){
      set.seed(seed)
    }
    optim.init <- c(runif(1, -1 , 1), rgamma(1, 2))
  }
  objFunMin <- optim(optim.init, objFun, method = "L-BFGS-B",
                     lower = c(-0.999, 1e-3),
                     upper = c(0.999, Inf))
  rhoHat     = objFunMin$par[1]
  s2Hat      = objFunMin$par[2]

  ZHatStar = PP %*% (ZZ - rhoHat*A%*%ZZ)
  DeltaHat = drop(qr.solve(crossprod(ZHatStar))%*%crossprod(ZHatStar,y - rhoHat*A%*%y))

  return(list(coefs=DeltaHat,s2=s2Hat,rho=rhoHat))
}


# Derive matrix normal approximation to prior of U
Uprior_appox <- function(Usample, maxiter = 100, eps = 1e-8){
  K <- dim(Usample)[1]
  n <- dim(Usample)[2]
  D <- dim(Usample)[3]

  Lambda <- apply(Usample, MARGIN=c(2,3), FUN=mean)

  Psi    <- diag(1, D)
  Omega  <- diag(1, n)

  ll    <- numeric(maxiter)
  ll[1] <- llMatNorm(Usample, Lambda, Omega, Psi)

  for(it in 2:maxiter){
    Psi_inv <- chol2inv(chol(Psi))
    S   <- matrix(0, nrow=n, ncol=n)
    for(j in 1:K){
      S <- S + tcrossprod((Usample[j,,] - Lambda)%*%Psi_inv,
                          Usample[j,,] - Lambda)
    }
    Omega  <- S/S[1,1]
    eta    <- S[1,1]/(K*D)
    Omega[-1,-1] <- eta*Omega[-1,-1] + (1-eta)*tcrossprod(Omega[-1, 1])
    Omega_inv <- chol2inv(chol(Omega))

    Psi <- matrix(0, nrow=D, ncol=D)
    for(j in 1:K){
      Psi <- Psi + crossprod((Usample[j,,] - Lambda),
                             Omega_inv%*%Usample[j,,] - Lambda)/(K*n)
    }

    ll[it] <- llMatNorm(Usample, Lambda, Omega, Psi)
    if( abs((ll[it] - ll[it-1])/ll[it-1]) < eps ) break
  }


  return(list(Lambda = Lambda,
              Omega  = Omega,
              Psi    = Psi,
              conv   = (it <= maxiter)))
}

llMatNorm <- function(X, M, U, V){
  K <- dim(X)[1]
  n <- dim(X)[2]
  D <- dim(X)[3]

  U_inv <- chol2inv(chol(U))
  V_inv <- chol2inv(chol(V))

  0.5*n*K*determinant(V, logarithm = T)$modulus +
    0.5*D*K*determinant(U, logarithm = T)$modulus -
    0.5*sum( apply(X, 1,
                   function(x){sum(diag(tcrossprod(V_inv, x - M)%*%
                                          U_inv%*%(x - M)))}) )
}
