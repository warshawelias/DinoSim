library(evoTS)

## Create a multivariate dataset
Rmatrix = as.matrix(read.csv("2mod.csv", header = F))

##Custom function
#annotations at modifications from original
sim.multi.URW.update <- function (ns = 30, anc = c(0, 0), R = matrix(c(0.5, 0, 0, 0.5), 
                                                                     nrow = 2, byrow = TRUE), vp = 0.1, nn = rep(30, ns), 
                                  tt = 0:(ns - 1)) 
{
  m <- ncol(R) #number of traits
  MM <- matrix(nrow = ncol(R), ncol = ns)
  mm <- matrix(nrow = ncol(R), ncol = ns)
  vv <- matrix(nrow = ncol(R), ncol = ns)
  time <- tt/max(tt)
  dt <- diff(time)
  Chol <- chol(R)
  temp <- MASS::mvrnorm(n = (ns), mu = rep(0, m), Sigma = ((t(Chol) %*% 
                                                              Chol) * (dt[1])))
  for (i in 1:m) {
    MM[i, c(1:ns)] <- cumsum(c(temp[, i]))
    mm[i, c(1:ns)] <- MM[i, c(1:ns)] + rnorm(ns, 0, sqrt(vp/nn))
    vv[i, c(1:ns)] <- rep(vp, (ns)) + rnorm(ns, 0, 0.5) ##add stochastic component to variance
  }
  MM <- MM + anc
  mm <- mm + anc
  List <- list()
  for (i in 1:m) {
    List[[i]] <- paleoTS::as.paleoTS(mm = mm[i, ], vv = vv[i, 
    ], nn = nn, tt = time, MM = MM[i, ], label = "Created by sim.multi.BM", 
    reset.time = FALSE)
  }
  if (m == 2) 
    yy <- make.multivar.evoTS(List[[1]], List[[2]])
  if (m == 3) 
    yy <- make.multivar.evoTS(List[[1]], List[[2]], List[[3]])
  if (m == 4) 
    yy <- make.multivar.evoTS(List[[1]], List[[2]], List[[3]], 
                              List[[4]])
  if (m == 5) 
    yy <- make.multivar.evoTS(List[[1]], List[[2]], List[[3]], 
                              List[[4]], List[[5]])
  if (m == 10) # allow 10 traits
    yy <- make.multivar.evoTS(List[[1]], List[[2]], List[[3]], 
                              List[[4]], List[[5]], List[[6]], 
                              List[[7]], List[[8]], List[[9]], 
                              List[[10]])
  return(yy)
}
data_set<-sim.multi.URW.update(ns=50, anc=c(1:10),R = Rmatrix, vp=2, nn=rep(30,50),tt=(0:49))

## plot the data
plotevoTS.multivariate(data_set)
