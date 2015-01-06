\name{modlist}
\docType{data}
\alias{modlist}
\title{List of model structures as per Clark et al., 2011}
\description{
'modlist' is the list of model structures implemented by Clark et al. (2011).
}
\usage{data}
\format{data frame}
\details{
In modlist each row identifies a model structure, and each column identifies a different element of the model structure:

\tabular{ll}{
first column   \tab mid   = model id number\cr
second column  \tab rferr = type of rainfall error (optional)\cr
third column   \tab arch1 = architecture of the upper soil layer\cr
fourth column  \tab arch2 = architecture of the lower soil layer\cr
fifth column 	 \tab qsurf = surface runoff\cr
sixth column 	 \tab qperc = vertical drainage\cr
seventh column \tab esoil = evapotranspiration\cr
eighth column	 \tab qintf = interflow\cr
ninth column 	 \tab q_tdh = routing\cr
}

For each element of the model structure, several model decisions can be made (see Clark et al. 2011 for details):
\itemize{
	\item Rainfall error (rferr) \itemize{
	           \item additive = 11
	           \item multiplicative = 12
	           }
	\item Architecture of the upper soil layer (arch1) \itemize{
	            \item Single state = 21
	            \item Separate tension storage = 22
	            \item Cascading buckets = 23
	            }
	\item Architecture of the lower soil layer (arch2) \itemize{
	            \item Baseflow reservoir of fixed size = 31
	            \item Tension reservoir plus two parallel tanks = 32
	            \item Baseflow reservoir of unlimited size (frac rate) = 33
	            \item Baseflow reservoir of unlimited size (power recession) = 34
	            %% \item Topmodel exponential reservoir (no more in the list) = 35
	            }
	\item Runoff (qsurf) \itemize{
	            \item Unsaturated zone Pareto = 41
	            \item Unsaturated zone linear = 42
	            \item Saturated zone topographic = 43
	            }
	\item Percolation (qperc) \itemize{
	            \item Drainage above field capacity = 51
	            \item Gravity drainage = 52
	            \item Saturated zone control = 53
	            }
	\item Evaporation (esoil) \itemize{
	            \item Root weighting = 61
	            \item Sequential = 62
	            }
	\item Interflows (qintf) \itemize{
	            \item Interflow denied = 71
	            \item Interflow allowed = 72
	            }
	\item Routing (q_tdh) \itemize{
	            \item Routing denied = 81
	            \item Routing allowed using Gamma distribution = 82
	            }
}

For instance, model 5 is identified by the following string (from the table above reading first the row number then column number for each model decision):
5 11 22 33 41 51 62 71 82
}
%\source{}
\references{
Clark M. P., SlaterA. G., Rupp D. E., Woods R. A., Vrugt J. A., Gupta H. V., Wagener T. and Hay L. E. (2008), Framework for 	Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models, Water Resour. Res. 44 p. 91-94

Clark M. P., McMillan H. K., Collins D. B. G., Kavetski D. and Woods R. A. (2011), Hydrological field data from a modeller's perspective: Part 2: process-based evaluation of model hypotheses. Hydrological Processes, 25: 523-543. doi: 10.1002/hyp.7902
}
\keyword{datasets}

