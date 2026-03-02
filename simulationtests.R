####FROM SCRATCH
sim.pop <- function(n, Bs, Bm, Bx) { #ontogeny with a single segment/individual
 y <- matrix(nrow =  n, ncol = 4)
 colnames(y) <- c("Juvenile Shape", "Slope", "Duration", "Adult Shape")
  for (i in 1:n) {
y[i,1] <- Bs + rnorm(1, mean = 0, sd = 1) #mean juvenile shape plus variation
y[i,2] <- Bm + rnorm(1, mean = 5, sd = 1) #mean slope plus variation
y[i,3] <- Bx + rnorm(1, mean = 5, sd = 1) #mean duration plus variation
y[i,4] <- y[i,1] + y[i,2] * y[i,3] #set of adult morphologies
  }
  return(y)
}

sim.pop.multiT <- function(n, Bs, Bm, Bx, t) {  #multiple distinct segments of trajectory
  Y <-vector("list", length = t)
  Y[[1]] <- matrix(nrow = n, ncol = 4)
  colnames(Y[[1]]) <- c("Juvenile Shape", "Slope", "Duration", "End Shape")
  for (i in 1:n) {
    Y[[1]][i,1] <- Bs + rnorm(1, mean = 0, sd = 1)
    Y[[1]][i,2] <- Bm[[1]] + rnorm(1, mean = 0, sd = 1)
    Y[[1]][i,3] <- Bx[[1]] + rnorm(1, mean = 0, sd = 1)
    Y[[1]][i,4] <- Y[[1]][i,1] + Y[[1]][i,2] * Y[[1]][i,3]
  }
  for (e in 2:(t)) {
    Y[[e]] <- matrix(nrow = n, ncol = 4)
    colnames(Y[[e]]) <- c("Start Shape", "Slope", "Duration", "End Shape")
    for (i in 1:n) {
      Y[[e]][i,1] <- Y[[e-1]][i,4]
      Y[[e]][i,2] <- Bm[[e]] + rnorm(1, mean = 0, sd = 1)
      Y[[e]][i,3] <- Bx[[e]] + rnorm(1, mean = 0, sd = 1)
      Y[[e]][i,4] <- Y[[e]][i,1] + Y[[e]][i,2] * Y[[e]][i,3]
    }
  }
  colnames(Y[[t]])[4] <- c("Adult Shape")
  return(Y)
}

make.multiTplotset <- function(n,t,x) { #isolate simulated data necessary for plotting
  X  <- matrix(ncol = t+1,  nrow = n)
  for (i in 1:t){
    X[1:n,i] <- x[[i]][1:n,1]
  }
  X[1:n,t+1] <- x[[t]][1:n,4]
  Y <- matrix(ncol = t+1, nrow = n)
  Y[1:n,1]<- rep(0,n)
  for(i in 2:(t+1)){
    Y[1:n,i] <- x[[i-1]][1:n,3] + Y[1:n,(i-1)]
  }
  set <- cbind(X,Y)
  return(set)
}
testset <- make.multiTplotset(30,3,testmultiT)

plotbouquet <-function(data, pcolor, bgcolor){ #plot data with differently colored points for each segment and lines connecting total trajectories  
  x <- t(data[1:nrow(data),1:(ncol(data)/2)])
  y <- t(data[1:nrow(data),((ncol(data)/2)+1):ncol(data)])
  plot(x,y, xlab = "Shape", ylab = "Age")
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray85") # Color
  grid(nx=NULL, ny=NULL, col="white", lty="dotted", lwd=2)
  par(new=TRUE)
  matlines(x,y,type="l",lwd=1, col = "black",lty=2)
  for (i in 1:length(bgcolor)) {
    points(x[i,1:(nrow(data))], y[i,1:nrow(data)], pch = 21, col = pcolor, bg = bgcolor[i])
  }
}

plotbouquet(testset,"darkblue", c("darkslategray2","darkslategray3","darkslategray4","darkslategrey"))

baseline <- sim.pop.multiT(n = 100, Bs = 1, Bm = c(5,5,5), Bx = c(5,5,5), t = 3)
baseplotset <- make.multiTplotset(100,3,baseline)
plotbouquet(baseplotset,"darkblue",c("darkslategray2","darkslategray3","darkslategray4","darkslategrey"))

steep <- sim.pop.multiT(n = 100, Bs = 1, Bm = c(10,10,10), Bx = c(5,5,5), t = 3) #steep ontogenetic trajectories
steepplotset <- make.multiTplotset(100,3,steep)
plotbouquet(steepplotset,"darkblue",c("darkslategray2","darkslategray3","darkslategray4","darkslategrey"))

long <- sim.pop.multiT(n = 100, Bs = 1, Bm = c(5,5,5), Bx = c(10,10,10), t = 3) #long ontogenetic trajectories
longplotset <- make.multiTplotset(100,3,long)
plotbouquet(longplotset,"darkblue",c("darkslategray2","darkslategray3","darkslategray4","darkslategrey"))
