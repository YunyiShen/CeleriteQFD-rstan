#' Sample from Pareto distribution
#' Take n samples from the Pareto distribution, kind of power law distribution
#' @param n number of samples
#' @param xm lower bound of the distribution, default 50, should be >0
#' @param alpha the power it decays, default 1
#' @return a vector of n samples
#' @export
rPareto <- function(n,xm=50,alpha=1){
    U <- runif(n)
    xm/((U)^(1/alpha))
}

kepler_raising <- function(x){
    1+1.941 * x  - 0.175 * x^2-2.246 * x^3 -1.125 * x^4
}

kepler_decay <- function(x){
    0.6890 * exp(-1.6 * x) + 0.3030 * exp(-0.2783 * x)
}

#' Simulate n Kepler flares
#' Routine simulating n Kepler flares in the given time vector tt
#' @param tt the given time vector where flares were put
#' @param t_half the characteristic time scale of the flares
#' @param n number of flares to be simulated
#' @param flux_dist a function that peak value of the flare was taken, default is Pareto
#' @param ... extra arguments to flux_dist
#' @return a list with vector of (absolute) flare flux in `$flux` and states (i.e. normal=1, increasing=2 and decay=3) in `$states`
#' @details The characteristic time scale is used for all flares, for each individual flare, its own local scale was the product of its peak value and the overall characteristic timescale. Peak flux was sampled from `flux_dist`.
#' @export
kepler_flare <- function(tt,t_half=.00005, n=5 ,flux_dist = rPareto, ...){
    flare <- 0 * tt
    states <- 1 + flare # store states
    flux_all <- flux_dist(n = n, ...)
    peak_time_all <- sample(tt, n)
    peak_time_all <- sort(peak_time_all)
    #browser()
    for(i in 1:n){
        flux_loc <- flux_all[i]
        t_half_loc <- t_half * flux_loc
        peak_time_loc <- peak_time_all[i]

        raising_phase <- which(peak_time_loc-tt >=0 & peak_time_loc-tt <= t_half_loc )
        decaying_phase <- which(tt - peak_time_loc >0 & tt - peak_time_loc <= 10 * t_half_loc )

        states[raising_phase] <- 2
        states[decaying_phase] <- 3

        raising_time_loc <- (tt[raising_phase]-peak_time_loc)/t_half_loc
        decaying_time_loc <- (tt[decaying_phase]-peak_time_loc)/t_half_loc

        raising_flux <- flux_loc * kepler_raising(raising_time_loc)
        decaying_flux <- flux_loc * kepler_decay(decaying_time_loc)
        #browser()
        flare[raising_phase] <- raising_flux
        flare[decaying_phase] <- decaying_flux
    }

    return(list(flare = flare, states = states))

}
