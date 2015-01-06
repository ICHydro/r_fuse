## FUSE Soil Moisture Accounting model
# Author: Claudia Vitolo
#
# Args:
#   DATA:                          matrix containing 3 columns called: P (precipitation), E (potential evapotranspiration) and Q (observed streamflow discharge, optional)
#   mid:                           model id
#   modlist:                       list of model structures ordered by model id
#   states:                        boolean. If states=TRUE, the output contains the list of state variables
#   fluxes:                        boolean. If fluxes=TRUE, the output contains the list of fluxes (last element of the list is U)
#   fracstate...qb_powr:           model parameters (collected as mparam)
#
# Returns:
#   U:                             Instantaneous runoff
#   s:                             (optional) list of state variables
#   f:                             (optional) list of fluxes (containing also U)

fusesma.sim <- function(DATA,mid,deltim,
                        frchzne,fracten,maxwatr_1,percfrac,
                        fprimqb,qbrate_2a,qbrate_2b,
                        qb_prms,maxwatr_2,baserte,rtfrac1,percrte,
                        percexp,sacpmlt,
                        sacpexp,iflwrte,axv_bexp,sareamax,loglamb,tishape,
                        qb_powr,
                        fracstate0 = 0.25,
                        rferr_add = 0,
                        rferr_mlt = 1,
                        absError = 10^(-4),         # absolute solver error (default)
                        relError = 10^(-4),         # relative solver error (default)
                        StatesFluxes = FALSE,
                        correctNegVols = FALSE){    # TODO: add module for correction of negative values

    ## Keep attributes, but work with raw matrix
    inAttr <- attributes(DATA[, 1])
    DATA<-coredata(DATA)

    stopifnot(c("P","E") %in% colnames(DATA))
    P <- DATA[,"P"]
    E <- DATA[,"E"]

    ## skip over missing values
    #bad <- is.na(P) | is.na(E)
    #P[bad] <- 0
    #E[bad] <- 0

    # load list of availabe models
    load(system.file("data/modlist.rda", package = "fuse"))
    modlist <- modlist # to remove NOTE in R CMD check

    # Make sure mid is an integer (needed for compatibility with hydromad)
    mid <- round(mid,0)

    # Read model structure [LIST]
    smodl <- list("rferr"=modlist[mid,2],
                "arch1"=modlist[mid,3],
                "arch2"=modlist[mid,4],
                "qsurf"=modlist[mid,5],
                "qperc"=modlist[mid,6],
                "esoil"=modlist[mid,7],
                "qintf"=modlist[mid,8],
                "q_tdh"=modlist[mid,9])

    res <- .C(  "fusesmaSimR",
                as.double(P),
                as.integer(length(P)),
                as.double(E),
                as.integer(length(E)),
                as.integer(smodl), # A (1d) array of ints (the model params)
                as.integer(length(smodl)),
                as.double(deltim),
                as.double(frchzne),
                as.double(fracten),
                as.double(maxwatr_1),
                as.double(percfrac),
                as.double(fprimqb),
                as.double(qbrate_2a),
                as.double(qbrate_2b),
                as.double(qb_prms),
                as.double(maxwatr_2),
                as.double(baserte),
                as.double(rtfrac1),
                as.double(percrte),
                as.double(percexp),
                as.double(sacpmlt),
                as.double(sacpexp),
                as.double(iflwrte),
                as.double(axv_bexp),
                as.double(sareamax),
                as.double(loglamb),
                as.double(tishape),
                as.double(qb_powr),
                stateResults = double( 10 * length(P) ),
                fluxResults = double( 20 * length(P) ),
                as.double(absError),
                as.double(relError),
                as.logical(correctNegVols),
                as.double(fracstate0),
                as.double(rferr_add),
                as.double(rferr_mlt),
                status = integer(1) , PACKAGE="fuse" )

    if (res$status != 0){
      stop("fusesmaSimR Failed!")
    }

    s <- res$stateResults
    dim(s) <- c(length(P),10)
    s <- data.frame( matrix(s,nrow=length(P),ncol=10,byrow=T) )

    f <- res$fluxResults
    dim(f) <- c(length(P),20)
    f <- data.frame( matrix(f,nrow=length(P),ncol=20,byrow=T) )

    results <- cbind(s,f)
    names(results) <- c("tens_1a","tens_1b","tens_1","free_1",
                        "watr_1","tens_2","free_2a","free_2b","watr_2","free_2",
                        "eff_ppt","satarea","qsurf","evap_1a","evap_1b","evap_1",
                        "evap_2","rchr2excs","tens2free_1","tens2free_2","qintf_1",
                        "qperc_12","qbase_2","qbase_2a","qbase_2b","oflow_1","oflow_2",
                        "oflow_2a","oflow_2b","U")

    if (StatesFluxes == FALSE) {
      results <- results$U
      attributes(results) <- inAttr
    }

    return(results)
}

fusesma.ranges <- function() {
    list("rferr_add" = 0,              # additive rainfall error (mm)
         "rferr_mlt" = 1,              # multiplicative rainfall error (-)
         "maxwatr_1" = c(25, 500),     # depth of the upper soil layer (mm)
         "maxwatr_2" = c(50, 5000),    # depth of the lower soil layer (mm)
         "fracten"   = c(0.05, 0.95),  # fraction total storage in tension storage (-)
         "frchzne"   = c(0.05, 0.95),  # fraction tension storage in recharge zone (-)
         "fprimqb"   = c(0.05, 0.95),  # fraction storage in 1st baseflow reservoir (-)
         "rtfrac1"   = c(0.05, 0.95),  # fraction of roots in the upper layer (-)
         "percrte"   = c(0.01, 1000),  # percolation rate (mm day-1)
         "percexp"   = c(1, 20),       # percolation exponent (-)
         "sacpmlt"   = c(1, 250),      # SAC model percltn mult for dry soil layer (-)
         "sacpexp"   = c(1, 5),        # SAC model percltn exp for dry soil layer (-)
         "percfrac"  = c(0.05, 0.95),  # fraction of percltn to tension storage (-)
         "iflwrte"   = c(0.01, 1000),  # interflow rate (mm day-1)
         "baserte"   = c(0.001, 1000), # baseflow rate (mm day-1)
         "qb_powr"   = c(1, 10),       # baseflow exponent (-)
         "qb_prms"   = c(0.001, 0.25), # baseflow depletion rate (day-1)
         "qbrate_2a" = c(0.001, 0.25), # baseflow depletion rate 1st reservoir (day-1)
         "qbrate_2b" = c(0.001, 0.25), # baseflow depletion rate 2nd reservoir (day-1)
         "sareamax"  = c(0.05, 0.95),  # maximum saturated area (-)
         "axv_bexp"  = c(0.001, 3),    # ARNO/VIC "b" exponent (-)
         "loglamb"   = c(5, 10),       # mean value of the topographic index (m)
         "tishape"   = c(2, 5))        # shape param for the topo index Gamma dist (-)
}

