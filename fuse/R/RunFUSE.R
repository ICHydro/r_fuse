#' Combine FUSE soil moisture accounting module and routing module in one function
#'
#' @param DATA data.frame containing observations. It consists of 3 columns:
#' Rainfall (P), Potential Evapo-Transpiration (E) and Streamflow (Q)
#' @param ParameterSet list of parameters
#' @param deltim observation time step (days)
#' @param mid model id number in Model List 2011(see below for details)
#'
#' @details The list of parameters can be generated as follows: ParameterSet <- GeneratePsetsFUSE(1).
#'
#' @return Simulated streamflow discharge
#'

RunFUSE <- function(DATA, ParameterSet, deltim, mid){

  U <- fusesma.sim(DATA,
                   mid,
                   deltim,
                   ParameterSet$frchzne,
                   ParameterSet$fracten,
                   ParameterSet$maxwatr_1,
                   ParameterSet$percfrac,
                   ParameterSet$fprimqb,
                   ParameterSet$qbrate_2a,
                   ParameterSet$qbrate_2b,
                   ParameterSet$qb_prms,
                   ParameterSet$maxwatr_2,
                   ParameterSet$baserte,
                   ParameterSet$rtfrac1,
                   ParameterSet$percrte,
                   ParameterSet$percexp,
                   ParameterSet$sacpmlt,
                   ParameterSet$sacpexp,
                   ParameterSet$iflwrte,
                   ParameterSet$axv_bexp,
                   ParameterSet$sareamax,
                   ParameterSet$loglamb,
                   ParameterSet$tishape,
                   ParameterSet$qb_powr)

  Q <- fuserouting.sim(U,
                       mid,
                       deltim,
                       timedelay=ParameterSet$timedelay)

  return(Q)

}
