####FROM SCRATCH
sim.pop <- function(n, Bs, Bm, Bx) {
 y <- matrix(nrow =  n, ncol = 4)
 colnames(y) <- c("Juvenile Shape", "Slope", "Duration", "Adult Shape")
  for (i in 1:n) {
y[i,1] <- Bs + rnorm(1, mean = 0, sd = 1) #mean juvenile shape plus variation
y[i,2] <- Bm + rnorm(1, mean = 0, sd = 1) #mean slope plus variation
y[i,3] <- Bx + rnorm(1, mean = 0, sd = 1) #mean duration plus variation
y[i,4] <- y[i,1] + y[i,2] * y[i,3] #set of adult morphologies
  }
  return(y)
}

test <- sim.pop(30, 1, 2, 6)
test
test2 <-  sim.pop(30,1,2,3)
test2

testlist <- vector("list", length=2)
testlist[[1]] <- test
testlist[[2]] <- test2

testlist[[1]][8,2]

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

testmultiT  <- sim.pop.multiT(n = 30, Bs = 1, Bm = c(1,2,3), Bx = c(1,2,3), t = 3) #slope and duration arguments must be lists of length t (number of segments)
testmultiT
