Framework for Understanding Structural Errors (FUSE, R package)
======

Implementation of the framework for hydrological modelling FUSE, based on the Fortran version described in Clark et al. (2008). The package consists of two modules: Soil Moisture Accounting module (fusesma.sim) and Gamma routing module (fuserouting.sim). It also contains default parameter ranges (fusesma.ranges and fuserouting.ranges) and three data objects: DATA, parameters and modlist.

### Basics
Install dependencies
```R
install.packages(c("zoo", "tgp", "BH", "Rcpp", "devtools"))
```

Install and load the FUSE package
```R
library(devtools)
install_github("cvitolo/r_fuse", subdir = "fuse")
library(fuse)
```

Load sample data (daily time step)
```R
data(DATA)
myDELTIM <- 1
```

Define parameter ranges
```R
DefaultRanges <- data.frame(t(data.frame(c(fusesma.ranges(),
                                           fuserouting.ranges()))))
names(DefaultRanges) <- c("Min","Max")
```

Sample parameter set using Latin Hypercube method
```R
numberOfRuns <- 100
parameters <- lhs( numberOfRuns, as.matrix(DefaultRanges) )
parameters <- data.frame(parameters)
names(parameters) <- row.names(DefaultRanges)
```

Alternatively, sample parameter set using built-in function
```R
parameters <- GenerateFUSEParameters(100)
```

### Example usage with 1 model structure
Define the model to use, e.g. TOPMODEL (MID = 60)
```R
myMID <- 60
```

Use the built-in function to run FUSE for the 1st sampled parameter set
```R
x <- RunFUSE(DATA, parameters[1,], myDELTIM, myMID)

plot(x,xlab="",ylab="Streamflow [mm/day]")
```

Run FUSE for all the sampled parameter sets 
```R
plot(DATA$Q,type="l",xlab="",ylab="Streamflow [mm/day]")
allQ <- data.frame(matrix(NA,ncol=numberOfRuns,nrow=dim(DATA)[1]))
for (i in 1:numberOfRuns){
  allQ[,i] <- RunFUSE(DATA, parameters[i,], myDELTIM, myMID)
  lines(zoo(allQ[,i],order.by=index(DATA)),col="gray",lwd=0.1)
}
lines(DATA$Q,col="black")
```

### Ensemble example usage
Define a group of model structures to use
```R 
mids <- c(60, 230, 342, 426)
```
 
Run a multi-model calibration using the Nash-Sutcliffe efficiency as objective function
```R
library(qualV)
indices <- rep(NA,4*numberOfRuns)
discharges <- matrix(NA,ncol=4*numberOfRuns,nrow=dim(DATA)[1])
kCounter <- 0

for (m in 1:4){

  myMID <- mids[m]

  for (pid in 1:numberOfRuns){

    kCounter <- kCounter + 1
    ParameterSet <- as.list(parameters[pid,])
    
    Qrout <- RunFUSE(DATA, parameters[pid,], myDELTIM, myMID)
 
    indices[kCounter] <- EF(DATA$Q,Qrout)  
    discharges[,kCounter] <- Qrout
    
    }
}
```

Compare results
```R 
bestRun <- which(indices == max(indices))
 
bestModel <- function(runNumber){
 if (runNumber<(numberOfRuns+1)) myBestModel <- "TOPMODEL"
 if (runNumber>(numberOfRuns+1) & runNumber<(2*numberOfRuns+1)) myBestModel <- "ARNOXVIC"
 if (runNumber>(2*numberOfRuns+1) & runNumber<(3*numberOfRuns+1)) myBestModel <- "PRMS"
 if (runNumber>(3*numberOfRuns+1) & runNumber<(4*numberOfRuns+1)) myBestModel <- "SACRAMENTO"
 return(myBestModel)
}
bestModel(bestRun)
 
plot(coredata(DATA$Q),type="l",xlab="",ylab="Streamflow [mm/day]", lwd=0.5)
 
for(pid in 1:(4*numberOfRuns)){
 lines(discharges[,pid], col="gray", lwd=3)
}
 
lines(coredata(DATA$Q),col="black", lwd=1)
lines(discharges[,bestRun],col="red", lwd=1)
```

How the best simulation of each model structure compare to each other?
```R 
bestRun0060 <- which(indices[1:numberOfRuns] == max(indices[1:numberOfRuns]))
bestRun0230 <- numberOfRuns + which(indices[(numberOfRuns+1):(2*numberOfRuns)] == max(indices[(numberOfRuns+1):(2*numberOfRuns)]))
bestRun0342 <- 2*numberOfRuns + which(indices[(2*numberOfRuns+1):(3*numberOfRuns)] == max(indices[(2*numberOfRuns+1):(3*numberOfRuns)]))
bestRun0426 <- 3*numberOfRuns + which(indices[(3*numberOfRuns+1):(4*numberOfRuns)] == max(indices[(3*numberOfRuns+1):(4*numberOfRuns)]))
 
plot(coredata(DATA$Q),type="l",xlab="",ylab="Streamflow [mm/day]", lwd=1)
lines(discharges[,bestRun0060], col="green", lwd=1)
lines(discharges[,bestRun0230], col="blue", lwd=1)
lines(discharges[,bestRun0342], col="pink", lwd=1)
lines(discharges[,bestRun0426], col="orange", lwd=1)
 
legend("top", 
        c("TOPMODEL", "ARNOXVIC", "PRMS","SACRAMENTO"), 
        col = c("green", "blue", "pink", "orange"),
        lty = c(1, 1, 1, 1))
```

# Leave your feedback
I would greatly appreciate if you could leave your feedbacks via email (cvitolodev@gmail.com).
