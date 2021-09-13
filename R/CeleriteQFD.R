#' CeleriteQFD detrending and flare detection with stan
#' @export
#' @param tt the time vector
#' @param flux the flux data, make sure there is no NA's
#' @param kernel the celerite kernel to be used, currently can only be one of "Rot" for rotation, "SHO" for SHO kernel and "QP" for quasi-periodic kernel. Default "Rot".
#' @param priors a list of priors, see details for details on different kernels
#' @param ... extra argument to stan::sampling
#' @return a list with stan's output in `$stan_out` and majority decoding of states at `$viterbi` and trend in `$trend`
CeleriteQFD <- function(tt, flux, kernel = "Rot", priors=NULL, ...){
  if(!kernel %in% c("Rot","SHO","GP")){
    stop("kernel can only be Rot, SHO or GP")
  }
  tt <- na.omit(tt)
  flux <- na.omit(flux)
  if(length(tt)!=length(flux)){
    stop("length mismatch of tt and flux, check for dimensions and NA's")
  }
  N <- length(tt)
  priors_given <- names(priors)
  if(kernel == "Rot"){
    QFD_data <- QFD_data <- list(N=N, t = tt,
                                 y = flux,
                                 sigma_prior = c(-8,8),
                                 Q0_prior = c(-8,8),
                                 dQ_prior = c(-8,8),
                                 period_prior = c(-8,8),
                                 f_prior = c(1e-6,1-1e-6),
                                 alpha_quiet = c(1,.1),
                                 alpha_firing = c(1,1),
                                 alpha_decay = c(1,.1,1),
                                 mu0_quiet = 0,
                                 lambda_quiet = .01,
                                 gamma_noise = c(0.01,0.01),
                                 mu0_rate_firing = 0,
                                 sigma_rate_firing = 1e3,
                                 mu0_rate_decay = 0,
                                 sigma_rate_decay = 1e3,
                                 diag = rep(1e-6,N)
    )
    for(pp in priors_given){
      `[`(QFD_data,pp) <- `[`(priors,pp)
    }

    stan_out <- rstan::sampling(stanmodels$CeleriteQFDexN, data = QFD_data, ...)
    offsets <- 22
  }
  if(kernel == "SHO"){
    QFD_data <- list(N=N, t = tt,
         y = flux,
         S0_prior = c(-5,5),
         w0_prior = c(-5,5),
         Q_prior = c(-5,5),
         alpha_quiet = c(1,.1),
         alpha_firing = c(1,1),
         alpha_decay = c(1,.1,1),
         mu0_quiet = 0,
         lambda_quiet = 1e-3,
         gamma_noise = c(0.01,0.01),
         mu0_rate_firing = 0,
         sigma_rate_firing = 1e3,
         mu0_rate_decay = 0,
         sigma_rate_decay = 1e3,
         diag = rep(1e-6,N)
    )
    for(pp in priors_given){
      `[`(QFD_data,pp) <- `[`(priors,pp)
    }
    stan_out <- rstan::sampling(stanmodels$CeleriteSHOQFDexN, data = QFD_data, ...)
    offsets <- 19
  }
  if(kernel == "QP"){
    QFD_data <- list(N=N, t = tt,
                     y = flux,
                     B_prior = c(-10,0),
                     L_prior = c(1.5,5),
                     P_prior = c(-3,5),
                     C_prior = c(-10,10),
                     alpha_quiet = c(1,.1),
                     alpha_firing = c(1,1),
                     alpha_decay = c(1,.1,1),
                     mu0_quiet = 0,
                     lambda_quiet = 1e-3,
                     gamma_noise = c(0.01,0.01),
                     mu0_rate_firing = 0,
                     sigma_rate_firing = 1e3,
                     mu0_rate_decay = 0,
                     sigma_rate_decay = 1e3,
                     diag = rep(1e-6,N)
    )
    for(pp in priors_given){
      `[`(QFD_data,pp) <- `[`(priors,pp)
    }
    stan_out <- rstan::sampling(stanmodels$CeleriteQPQFDexN, data = QFD_data, ...)
    QFD_samples <- as.data.frame(stan_out)
    offsets <- 21
  }
  summQFD <- summary(fitQFD)
  QFD_samples <- as.data.frame(stan_out)
  Viterbi_raw <- QFD_samples[,1:(N-1) + (N + offsets)]
  Viterbi_max <- apply(Viterbi_raw,2,majority)
  trend <- summQFD[[1]][1:N+2*N+offsets-1, 1]
  rm(QFD_samples)
  rm(Viterbi_raw)
  rm(summQFD)
  return(list(stan_out = stan_out,viterbi = Viterbi_max,kernel = kernel, trend = trend))
}

#' Celerite detrending with stan
#' @export
#' @param tt the time vector
#' @param flux the flux data, make sure there is no NA's
#' @param kernel the celerite kernel to be used, currently can only be one of "Rot" for rotation, "SHO" for SHO kernel and "QP" for quasi-periodic kernel. Default "Rot".
#' @param priors a list of priors, see details for details on different kernels
#' @param ... extra argument to stan::sampling
#' @return a list with stan's output in `$stan_out` and trend in `$trend`
Celerite <- function(tt, flux, kernel = "Rot", priors=NULL, ...){
  if(!kernel %in% c("Rot","SHO","GP")){
    stop("kernel can only be Rot, SHO or GP")
  }
  tt <- na.omit(tt)
  flux <- na.omit(flux)
  if(length(tt)!=length(flux)){
    stop("length mismatch of tt and flux, check for dimensions and NA's")
  }
  N <- length(tt)
  priors_given <- names(priors)
  if(kernel == "Rot"){
    QFD_data <- QFD_data <- list(N=N, t = tt,
                                 y = flux,
                                 sigma_prior = c(-8,8),
                                 Q0_prior = c(-8,8),
                                 dQ_prior = c(-8,8),
                                 period_prior = c(-8,8),
                                 f_prior = c(1e-6,1-1e-6),
                                 err_prior = c(.01,.01),
                                 diag = rep(1e-6,N)
    )
    for(pp in priors_given){
      `[`(QFD_data,pp) <- `[`(priors,pp)
    }

    stan_out <- rstan::sampling(stanmodels$celerite, data = QFD_data, ...)
    offsets <- 23
  }
  if(kernel == "SHO"){
    QFD_data <- list(N=N, t = tt,
                     y = flux,
                     S0_prior = c(-5,5),
                     w0_prior = c(-5,5),
                     Q_prior = c(-5,5),
                     err_prior = c(.01,.01),
                     diag = rep(1e-6,N)
    )
    for(pp in priors_given){
      `[`(QFD_data,pp) <- `[`(priors,pp)
    }
    stan_out <- rstan::sampling(stanmodels$celeriteSHO, data = QFD_data, ...)
    offsets <- 7
  }
  if(kernel == "QP"){
    QFD_data <- list(N=N, t = tt,
                     y = flux,
                     B_prior = c(-10,0),
                     L_prior = c(1.5,5),
                     P_prior = c(-3,5),
                     C_prior = c(-10,10),
                     error_prior = c(.01,.01),
                     diag = rep(1e-6,N)
    )
    for(pp in priors_given){
      `[`(QFD_data,pp) <- `[`(priors,pp)
    }
    stan_out <- rstan::sampling(stanmodels$celeriteQP, data = QFD_data, ...)
    QFD_samples <- as.data.frame(stan_out)
    offsets <- 9
  }
  summQFD <- summary(fitQFD)
  trend <- summQFD[[1]][[1]][1:N + (N+offsets),1]
  rm(summQFD)
  return(list(stan_out = stan_out,kernel = kernel, trend = trend))
}

