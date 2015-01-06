## FUSE Routing model
# Author: Claudia Vitolo
#
# Args:
#   U:  effective rainfall (sum of surface runoff, overflow, interflow, and baseflow)
#   mid:  model id number
#   modlist:  list of models
#   deltim: data time step: deltim = 1 (daily time step), 1/24 (hourly time step), 1/24/4 (15 min time step)
#   timedelay:  time delay in runoff (days) mparam$timedelay (use a gamma distribution with shape parameter <- 2.5)
#
# Returns:
#   X:                             routed runoff

fuserouting.sim <- function(U, mid, deltim, timedelay) {    
  
   inAttr <- attributes(U)
   U <- coredata(U)
   
   # load list of availabe models
   load(system.file("data/modlist.rda", package = "fuse"))
   modlist <- modlist # to remove NOTE in R CMD check
   
   # Read model structure
   smodl <- c(modlist[mid,2],    # rferr
              modlist[mid,3],    # arch1
              modlist[mid,4],    # arch2
              modlist[mid,5],    # qsurf
              modlist[mid,6],    # qperc
              modlist[mid,7],    # esoil
              modlist[mid,8],    # qintf
              modlist[mid,9])    # q_tdh
   
   res <- .C(  "fuseRoutingSimR",
               as.double(U),
               as.integer(length(U)),
               as.integer(smodl), # A (1d) array of ints (the model params)
               as.integer(length(smodl)),
               as.double(timedelay),
               as.double(deltim),
               X = double( length(U) ),
               as.integer(length(U)),
               status = integer(1), PACKAGE="fuse" )
   
   if (res$status != 0){
     stop("fuseRoutingSimR Failed!")
   }
   
   X <- res$X
   attributes(X) <- inAttr
 
   return(X)
}

fuserouting.ranges <- function() {
   list("timedelay" = c(0.01, 5))        # time delay in runoff (days)
}
    
