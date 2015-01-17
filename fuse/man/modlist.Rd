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
	           \item additive_e (additive rainfall error) = 11
	           \item multiplc_e (multiplicative rainfall error) = 12
	           }
	\item Architecture of the upper soil layer (arch1) \itemize{
	            \item onestate_1 (Single state: upper layer defined by a single state variable) = 21
	            \item tension1_1 (Separate tension storage: upper layer broken up into tension and free storage) = 22
	            \item tension2_1 (Cascading buckets: tension storage sub-divided into recharge and excess) = 23
	            }
	\item Architecture of the lower soil layer (arch2) \itemize{
	            \item fixedsiz_2 (baseflow reservoir of fixed size) = 31
	            \item tens2pll_2 (tension reservoir plus two parallel tanks) = 32
	            \item unlimfrc_2 (baseflow reservoir of unlimited size (0-HUGE), frac rate) = 33
	            \item unlimpow_2 (baseflow reservoir of unlimited size (0-HUGE), power recession) = 34
	            %% \item topmdexp_2 (topmodel exponential reservoir (no more in the list) = 35
	            }
	\item Runoff (qsurf) \itemize{
	            \item arno_x_vic (Unsaturated zone Pareto: ARNO/Xzang/VIC parameterization (upper zone control)) = 41
	            \item prms_varnt (Unsaturated zone linear: PRMS variant (fraction of upper tension storage)) = 42
	            \item tmdl_param (Saturated zone topographic: TOPMODEL parameterization (only valid for TOPMODEL qb)) = 43
	            }
	\item Percolation (qperc) \itemize{
	            \item perc_f2sat (Drainage above field capacity: water from (field cap to sat) avail for percolation) = 51
	            \item perc_lower (Gravity drainage: percolation defined by moisture content in lower layer (SAC)) = 52
	            \item perc_w2sat (Saturated zone control: water from wilt point to saturation is available for percolation) = 53
	            }
	\item Evaporation (esoil) \itemize{
	            \item rootweight (root weighting) = 61
	            \item sequential (sequential evaporation model) = 62
	            }
	\item Interflows (qintf) \itemize{
	            \item intflwnone (interflow denied) = 71
	            \item intflwsome (interflow allowed) = 72
	            }
	\item Routing (q_tdh) \itemize{
	            \item no_routing (routing denied) = 81
	            \item rout_gamma (routing allowed using Gamma distribution with shape parameter = 2.5) = 82
	            }
}

For instance, model 5 is identified by the following string (from the table above reading first the row number then column number for each model decision):
5 11 22 33 41 51 62 71 82

The function \code{readmd()} returns the model components corresponding to a given id. For example:

readmd(5)

"[11] rferr = additive_e (additive rainfall error)"
"[22] arch1 = tension1_1 (upper layer broken up into tension and free storage)"
"[33] arch2 = unlimfrc_2 (baseflow resvr of unlimited size (0-HUGE), frac rate)"
"[41] qsurf = arno_x_vic (ARNO/Xzang/VIC parameterization (upper zone control))"
"[51] qperc = perc_f2sat (water from (field cap to sat) avail for percolation)"
"[62] esoil = sequential (sequential evaporation model)"
"[71] qintf = intflwnone (no interflow)"
"[82] q_tdh = rout_gamma (use a Gamma distribution with shape parameter = 2.5)"

A similar function (\code{readmd2var()}) reads the components as either character or numeric vector. For example:

readmd2var(5)
       rferr        arch1        arch2        qsurf        qperc        esoil 
"additive_e" "tension1_1" "unlimfrc_2" "arno_x_vic" "perc_f2sat" "sequential" 
       qintf        q_tdh 
"intflwnone" "rout_gamma" 

readmd2var(5, number = TRUE)
rferr arch1 arch2 qsurf qperc esoil qintf q_tdh 
   11    22    33    41    51    62    71    82
   
On Clark's model list, the parent models (Topmodel, VIC, PRMS and Sacramento) have id = 60, 230, 342 and 426 respectively. 

}
%\source{}
\references{
Clark M. P., SlaterA. G., Rupp D. E., Woods R. A., Vrugt J. A., Gupta H. V., Wagener T. and Hay L. E. (2008), Framework for 	Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models, Water Resour. Res. 44 p. 91-94

Clark M. P., McMillan H. K., Collins D. B. G., Kavetski D. and Woods R. A. (2011), Hydrological field data from a modeller's perspective: Part 2: process-based evaluation of model hypotheses. Hydrological Processes, 25: 523-543. doi: 10.1002/hyp.7902
}
\keyword{datasets}

