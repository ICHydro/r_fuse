\name{fuserouting.sim}
\alias{fuserouting.sim}
\title{Implementation of the framework for hydrological modelling FUSE - Routing model}
\description{Routing module derived from the Fortran version of FUSE by Martyn Clark (2011)}
\usage{fuserouting.sim(U, mid, deltim, timedelay)}
\arguments{
  \item{U}{Instantaneous runoff computed using fusesma.sim}
  \item{mid}{model id number in Model List 2011(see below for details)}
  \item{timedelay}{time delay in runoff (days)}
  \item{deltim}{observation time step (days)}
}

\details{fuserouting.sim() is a routing module based on a two parameter Gamma distribution. It takes in input the instantaneous discharge and returns the routed discharge. It is compatible with the HYDROMAD framework (see hydromad package: http://hydromad.catchment.org/).}

\value{The function returns an array of simulated "routed" discharges. It can be used after calculating instantaneous discharges with \code{fusesma.sim}.}

\references{Clark M. P., SlaterA. G., Rupp D. E., Woods R. A., Vrugt J. A., Gupta H. V., Wagener T. and Hay L. E. (2008), Framework for Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models, Water Resour. Res. 44 p. 91-94 

Clark M. P., McMillan H. K., Collins D. B. G., Kavetski D. and Woods R. A. (2011), Hydrological field data from a modeller's perspective: Part 2: process-based evaluation of model hypotheses. Hydrological Processes, 25: 523-543. doi: 10.1002/hyp.7902}

\author{Claudia Vitolo, Imperial College London}
\seealso{\code{\link{fusesma.sim}}}

\examples{
data(DATA)

# returns the instantaneous runoff
Qinst <- fusesma.sim (DATA, mid = 650, deltim = 1, 
                      rferr_add = 0, rferr_mlt = 1, 
                      maxwatr_1 = 458.1, maxwatr_2 = 3147.06, 
                      fracten = 0.267488, frchzne = 0.888808, 
                      fprimqb = 0.642572, rtfrac1 = 0.477619, percrte = 358.852,
                      percexp = 10.8906, sacpmlt = 25.6919,sacpexp = 3.18839,
                      percfrac = 0.119675, iflwrte = 346.623,
                      baserte = 493.714, qb_powr = 7.29671, qb_prms = 0.210137,
                      qbrate_2a = 0.152566, qbrate_2b = 0.0501708, 
                      sareamax = 0.349326, axv_bexp = 2.08099, 
                      loglamb = 5.45848, tishape = 4.16578)
                        
# returns the routed runoff
Qrout <- fuserouting.sim(Qinst, mid=5, deltim=1, timedelay = 0.462784)

}

\keyword{models}
\keyword{fuse}
\keyword{routing}

