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

make.multiTplotset <- function(x) { #isolate simulated data necessary for plotting
  n <- nrow(x[[1]])
  t <- length(x)
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

plotbouquet <-function(data, pcolor, bgcolor, xlim=NULL, ylim=NULL, cex=NULL){ #plot data with differently colored points for each segment and lines connecting total trajectories  
  x <- t(data[1:nrow(data),1:(ncol(data)/2)])
  y <- t(data[1:nrow(data),((ncol(data)/2)+1):ncol(data)])
  plot(x,y, xlab = "Shape", ylab = "Age", xlim = xlim, ylim = ylim)
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray85") # Color
  grid(nx=NULL, ny=NULL, col="white", lty="dotted", lwd=2)
  par(new=TRUE)
  matlines(x,y,type="l",lwd=1, col = "black",lty=2)
  for (i in 1:length(bgcolor)) {
    points(x[i,1:(nrow(data))], y[i,1:nrow(data)], pch = 21, col = pcolor, bg = bgcolor[i], cex=cex)
  }
}

baseline <- sim.pop.multiT(n = 100, Bs = 1, Bm = c(5,5,5), Bx = c(5,5,5), t = 3)
baseplotset <- make.multiTplotset(baseline)
plotbouquet(baseplotset,"darkblue",c("darkslategray2","darkslategray3","darkslategray4","darkslategrey"), xlim =c(0,200), ylim = c(0,40),cex=2)

steep <- sim.pop.multiT(n = 100, Bs = 1, Bm = c(10,10,10), Bx = c(5,5,5), t = 3) #steep ontogenetic trajectories
steepplotset <- make.multiTplotset(steep)
plotbouquet(steepplotset,"darkblue",c("darkslategray2","darkslategray3","darkslategray4","darkslategrey"), xlim =c(0,200), ylim = c(0,40),cex=2)

long <- sim.pop.multiT(n = 100, Bs = 1, Bm = c(5,5,5), Bx = c(10,10,10), t = 3) #long ontogenetic trajectories
longplotset <- make.multiTplotset(long)
plotbouquet(longplotset,"darkblue",c("darkslategray2","darkslategray3","darkslategray4","darkslategrey"), xlim =c(0,200), ylim = c(0,40),cex=2)

extract.bouquet.var <- function(x, #adult trait data set
                                t,  #number of ontogeny segments
                                sd = 1, #standard deviation of stochastic components
                                Bmr, #range of slopes to test
                                ){
  
  Bm <- sample(Bmr, length(Bmr), size=1, replace=TRUE)
  sim.pop.multiT(n, 0, Bm, Bx, t)
}

bouquet.statmodel.SLOPE <- function(x, #adult trait data set
                                    t,  #number of ontogeny segments
                                    sd = 1, #standard deviation of stochastic components
                                    Bm0, #slope estimate
                                    iter #number of simulated data sets to create
                                    ){
  n = nrow(x)
  Bx = mean(x[1:n,2])
  simdata <- matrix(nrow = iter, ncol = 2)
  for(i in 1:iter){
    Bm <- Bm0 + rnorm(1, mean = 0, sd = 1)
    bouquet <- sim.pop.multiT(n,0,rep(Bm,t),rep(Bx,t),t)
    simdata[i,1] <- var(bouquet[[t]][1:n,4])
    simdata[i,2] <- Bm
  }
  regression <- lm(simdata[,1]~simdata[,2]) #correlation between ontogeny slope and adult trait variance
  estimate <- (var(x[,1]) - regression[["coefficients"]][1])/regression[["coefficients"]][2] #estimate of ontogeny slope from empirical data
  print(estimate)
  return(simdata)
}

statmodeltest <- matrix(nrow =  100, ncol = 2)
statmodeltest[1:100,1] <- baseline[[3]][1:100,4]
statmodeltest[1:100,2] <- baseline[[3]][1:100,3]

testmodel <- bouquet.statmodel.SLOPE(statmodeltest,3,Bm0=5,iter=1000)
