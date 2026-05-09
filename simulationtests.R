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

###Parameter Estimation
library(tidyverse)
library(tidybayes)
library(rethinking)
library(cmdstanr)
#rebuild_cmdstan(cores=2)
#set_cmdstan_path(path="C:/Users/warsh/.cmdstan/cmdstan-2.38.0")

ulambaseline <- data.frame(baseline[[3]][,4],(baseline[[1]][,3] + baseline[[2]][,3] + baseline[[3]][,3]))
colnames(ulambaseline) <- c("Shape","Age")
ulamsteep <- data.frame(steep[[3]][,4],(steep[[1]][,3] + steep[[2]][,3] + steep[[3]][,3]))
colnames(ulamsteep) <- c("Shape","Age")
ulamlong <- data.frame(long[[3]][,4],(long[[1]][,3] + long[[2]][,3] + long[[3]][,3]))
colnames(ulamlong) <- c("Shape","Age")

flist1T <- alist(Shape  ~ dnorm(a2,sigma),
                 a2 <- a + (Age * Slope),
                 sigma ~ dunif(5,15),
                 a ~ dnorm(0,1),
                 Slope ~ dnorm(10,5)    
)

test.ulambaseline <- ulam(flist1T,data=ulambaseline,cores = 4, chains = 4, iter = 2000)
precis(test.ulambaseline)
plot(test.ulambaseline)

test.ulamsteep <- ulam(flist1T, data=ulamsteep, cores = 4, chains =4, iter = 2000)
precis(test.ulamsteep)
plot(test.ulamsteep)

test.ulamlong <- ulam(flist1T, data = ulamlong, cores = 4, chains = 4, iter = 2000)
precis(test.ulamlong)
plot(test.ulamlong)

###Empirical Tests
library(geomorph)


croc.landmarks <- readland.tps(file = "Extant_Croc_Dorsal.tps", specID = "imageID")
croc.landmarks
crocgpa <- gpagen(croc.landmarks)
crocgpa
coords <- crocgpa$coords
coords
PCA <- gm.prcomp(coords)
PCAscores <- as.data.frame(PCA$x)
PCAscores

crocmetadata <-  read.csv(file = "Extant_Croc_Covariates.csv")
crocmetadata

croctotal <- cbind(crocmetadata,PCAscores,log(crocgpa$Csize))
colnames(croctotal)[colnames(croctotal) == "log(crocgpa$Csize)"] <- "logCsize"
niloticus <- filter(croctotal,Species == "C. niloticus")

niloticus.ulam <- niloticus[,14:38] #PC scores for all nile croc specimens
niloticus.adults <- filter(niloticus, Age == "Adult" | Age =="Subadult")
niloticus.adultsulam <- niloticus.adults[,14:38]

nilelistUV <- alist(Comp1 ~ dnorm(a2,sigma),
                    a2 <- a + (logCsize * Slope),
                    sigma ~ dunif(0,0.1),
                    a ~ dnorm(-0.02,0.5),
                    Slope ~ dnorm(10,5)
                    )
test.ulamnileUV <- ulam(nilelistUV, data=niloticus.adultsulam, cores = 4, chains = 4, iter = 2000)
precis(test.ulamnileUV)

nilelistMV <- alist(Comp1:Comp2 ~ multi_normal(adult, diag(2),sigma),
                    adult <- a + (logCsize * Slope),
                    sigma ~ dunif(0,1),
                    a ~ multi_normal(rep(-0.02,2),diag(2),rep(0.5,2)),
                    Slope ~ multi_normal(rep(5,2),diag(2),rep(2,2))
                    )
test.ulamnileMV <- ulam(nilelistMV, data=niloticus.ulam, cores = 4, chains = 4, iter = 2000)
precis(test.ulamnileMV)