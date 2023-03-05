brady <- read.csv("brady.csv")

### QUESTION 1 ###

# a)
hist(brady$LNG, main = "Histogram of Brady's longest passes",
     xlab = "Longest Passes in yards")

# b)
# calculate E[Y] and E[Y^2] for y=LNG
ey <- mean(brady$LNG)
ey2 <- (sum(brady$LNG))^2 / 346

# g)
y <- brady$LNG
createGammaPsiFn <- function(y) { 
  function(theta) {
    alpha <- theta[1]
    beta <- theta[2]
    c(-346*(digamma(alpha)+log(beta)) + sum(log(y)), 
      -346*(alpha/beta) + (1/beta^2)*(sum(y)))
  }
}

createGammePsiPrimeFn <- function(y) {
  function(theta) {
    alpha <- theta[1]
    beta <- theta[2]
    val <- matrix(0, nrow = 2, ncol = 2)
    val[1,1] = -346*(trigamma(alpha))
    val[1,2] = -346*(1/beta)
    val[2,1] = -346*(1/beta) 
    val[2,2] = -346*(alpha/beta^2) - (2/beta^3)*(sum(y))
    return(val)
  }
}

# h)

testConvergence <- function(thetaNew, thetaOld, tolerance = 1e-10, relative = FALSE) { 
  sum(abs(thetaNew - thetaOld)) < if (relative)
  tolerance * sum(abs(thetaOld)) else tolerance 
}

NewtonRaphson <- function(theta, psiFn, psiPrimeFn, dim, testConvergenceFn = testConvergence, 
                          maxIterations = 100, tolerance = 1e-06, relative = FALSE) {
  if (missing(theta)) {
    ## need to figure out the dimensionality
    if (missing(dim)) {
      dim <- length(psiFn()) 
    }
    theta <- rep(0, dim) 
  }
  converged <- FALSE 
  i <- 0
  
  while (!converged & i <= maxIterations) {
    thetaNew <- theta - solve(psiPrimeFn(theta), psiFn(theta))
    converged <- testConvergenceFn(thetaNew, theta, tolerance = tolerance,
                                   relative = relative)
    theta <- thetaNew
    i <- i + 1 
  }
  ## Return last value and whether converged or not
  list(theta = theta, converged = converged, iteration = i, fnValue = psiFn(theta)) 
}

nr <- NewtonRaphson(theta = c(0.00289855, 13586.62), 
                    psiFn = createGammaPsiFn(y), psiPrimeFn = createGammePsiPrimeFn(y))
print(nr)


# h)
xlen <- seq(0,100)
hist(brady$LNG, prob = TRUE, main = "Density Histogram of Brady's LNG",
     xlab = "LNG", ylim=c(0, 0.035))
lines(dgamma(x=xlen, shape=0.00289855, scale = 13586.62), col = "red")
lines(dgamma(x=xlen, shape=7.185768, scale = 5.499906), col = "blue")
legend("topright", legend=c("alpha = 0.00289855, beta = 13586.62",
                            "alpha = 7.185768, beta = 5.499906"), 
       col=c("red","blue"), lty=c(1,1), bty="n",
       cex=0.6)


### QUESTION 2 ###

irls <- function(y, x, theta, rhoPrimeFn, dim = 2, delta = 1E-10,
                 testConvergenceFn = testConvergence,
                 maxIterations = 100, tolerance = 1E-6, relative = FALSE
                 # maximum number of iterations
                 # parameters for the test
                 # for convergence function
) { 
  if (missing(theta)) {theta <- rep(0, dim)} ## Initialize
  converged <- FALSE
  i <- 0
  N <- length(y)
  wt <- rep(1,N)
  ## LOOP
  while (!converged & i <= maxIterations) {
    ## get residuals
    resids <- getResids(y, x, wt, theta)
    ## update weights (should check for zero resids) wt <- getWeights(resids, rhoPrimeFn, delta)
    ## solve the least squares problem
    thetaNew <- getTheta(y, x, wt)
    ##
    ## Check convergence
    converged <- testConvergenceFn(thetaNew, theta,
                                   tolerance = tolerance, relative = relative)
                                   ## Update iteration
    theta <- thetaNew
    i <- i + 1 
  }
  ## Return last value and whether converged or not
  list(theta = theta, converged = converged, iteration = i) 
}

getWeights <- function(resids, rhoPrimeFn, delta = 1e-12) {
  ## for calculating weights, minimum |residual| will be delta 
  smallResids <- abs(resids) <= delta
  ## take care to preserve sign (in case rhoPrime not symmetric) 
  resids[smallResids] <- delta * sign(resids[smallResids])
  ## calculate and return weights
  rhoPrimeFn(resids)/resids
}

getTheta <- function(y, x, wt) {
  theta <- numeric(length = 2)
  ybarw <- sum(wt * y)/sum(wt)
  xbarw <- sum(wt * x)/sum(wt)
  theta[1] <- ybarw - (sum(wt * (x - xbarw) * y)/sum(wt * (x - xbarw)^2))*(xbarw - mean(x))
  theta[2] <- sum(wt * (x - xbarw) * y)/sum(wt * (x - xbarw)^2)
  ## return theta
  theta
}
getResids <- function(y, x, wt, theta) {
  xbar <- mean(x)  
  alpha <- theta[1]
  beta <- theta[2]
  ## resids are
  y - alpha - beta * (x - xbar)
}

LSrhoPrime <- function(resid) {
  2 * resid 
}

LADrhoPrime <- function(resid) {
  sign(resid)
}

huber.fn.prime <- function(resid, k = 1.345) { 
  val = resid
  subr = abs(resid) > k
  val[subr] = k * sign(resid[subr]) 
  return(val)
}

tukey.fn.prime <- function(resid, k=4.685) {
  val = resid - (2 * resid^3)/(k^2) + (resid^5)/(k^4) 
  subr = abs(resid) > k
  val[subr] = 0
  return(val)
}

# a)
# irls function from class
IRLSresult <- irls(brady$YDS, brady$AVG, theta = c(1, 1), rhoPrimeFn = LSrhoPrime,
                   tolerance = 1e-10, maxIterations = 1000)
print(IRLSresult)

# b)
# irls function from class
IRLSresult2 <- irls(brady$YDS, brady$AVG, theta = c(264.89306, 28.43653), rhoPrimeFn = LADrhoPrime,
                    tolerance = 1e-10, maxIterations = 1000)
print(IRLSresult2)

# c)
# irls function from class
IRLSresult3 <- irls(brady$YDS, brady$AVG, theta = c(264.89306, 28.43653), rhoPrimeFn = huber.fn.prime,
                    tolerance = 1e-10, maxIterations = 1000)
print(IRLSresult3)

# d)
# irls function from class
IRLSresult4 <- irls(brady$YDS, brady$AVG, theta = c(264.89306, 28.43653), rhoPrimeFn = tukey.fn.prime,
                    tolerance = 1e-10, maxIterations = 1000)
print(IRLSresult4)

# e) 
plot(x = brady$AVG, y= brady$YDS, main = "Scatter plot of Brady's YDS vs AVG",
     xlab = "AVG", ylab = "YDS", pch = 16)
abline(a = IRLSresult$theta[1] - IRLSresult$theta[2] * mean(brady$AVG),
       b = IRLSresult$theta[2], col = "firebrick", lwd = 2, lty = 1)
abline(a = IRLSresult2$theta[1] - IRLSresult2$theta[2] * mean(brady$AVG),
       b = IRLSresult2$theta[2], col = "navyblue", lwd = 2, lty = 2)
abline(a = IRLSresult3$theta[1] - IRLSresult3$theta[2] * mean(brady$AVG),
       b = IRLSresult3$theta[2], col = "green", lwd = 2, lty = 3)
abline(a = IRLSresult4$theta[1] - IRLSresult4$theta[2] * mean(brady$AVG),
       b = IRLSresult4$theta[2], col = "purple", lwd = 2, lty = 4)
legend("bottomright", c("LS line", "LAD line", "Huber line", "Tukey line"),
       col = c("firebrick", "navyblue", "green", "purple"),lty = 1:4, lwd = 2,
       cex = 0.75, bty = "n")


### QUESTION 3 ###

# a)
brady2 <- read.csv("brady.csv", header = TRUE)
bradypw <- brady2[brady2$PS == 1 & brady2$RESULT == "W", ]

pop <- bradypw$YDS
samples5 <- combn(x = pop, m = 5)
print(paste0("Mean of 1st sample: ", mean(samples5[, 1]))) # mean of 1st sample
print(paste0("Mean of 3rd sample: ",mean(samples5[, 3]))) # mean of 3rd sample
print(paste0("Mean of 10th sample: ",mean(samples5[, 10]))) # mean of 10th sample

# b)

samples30 <- combn(x = pop, m = 30)
print(paste0("Mean of 1st sample: ", mean(samples30[, 1]))) # mean of 1st sample
print(paste0("Mean of 3rd sample: ",mean(samples30[, 3]))) # mean of 3rd sample
print(paste0("Mean of 10th sample: ",mean(samples30[, 10]))) # mean of 10th sample


# c)

# n = 5
avesSamp5 <- apply(samples5, MARGIN = 2, FUN = function(s) {
  mean(samples5[s]) 
})

sampleErrorsa5 <- avesSamp - avePop

# n = 30
avesSamp30 <- apply(samples30, MARGIN = 2, FUN = function(s) {
  mean(samples30[s]) 
})

# histogram plot
avePop <- mean(pop)
M5 <- ncol(samples5)
M30 <- ncol(samples30)

sampleErrorsa5 <- avesSamp5 - avePop
sampleErrorsa30 <- avesSamp30 - avePop

# average and std dev of sample errors (n=5)
# avg = -57.32353, std dev = 40.04048
meanSamErra5 <- mean(sampleErrorsa5)
sdSamErra5 <- sd(sampleErrorsa5)

# average and std dev of sample errors (n=30)
# avg = 9.352941, std dev = 5.574544
meanSamErra30 <- mean(sampleErrorsa30)
sdSamErra30 <- sd(sampleErrorsa30)

par(mfrow = c(1,2))
hist(sampleErrorsa5, prob = TRUE, main = "All possible sample errors (n=5) \n attribute of interest is the average", 
     cex.main=0.5,xlim=c(-175, 100), xlab = "Sample Error", col = adjustcolor("purple", 0.3))
abline(v = 0, col = "red")
legend("topleft", c("Avg: -57.32353", "Std Dev: 40.04048"), bty = "n", cex=0.33)
hist(sampleErrorsa30, prob = TRUE, main = "All possible sample errors (n=30) \n attribute of interest is the average", 
     cex.main=0.5, xlim=c(-175, 100), xlab = "Sample Error")
abline(v = 0, col = "red")
legend("topleft", c("Avg: 9.352941", "Std Dev: 5.574544"), bty = "n", cex=0.4)


# d)

# n = 5
medSamp5 <- apply(samples5, MARGIN = 2, FUN = function(s) {
  median(samples5[s]) 
})

# n = 30
medSamp30 <- apply(samples30, MARGIN = 2, FUN = function(s) {
  median(samples30[s]) 
})

# histogram plot
medPop <- median(pop)

sampleErrorsb5 <- medSamp5 - medPop
sampleErrorsb30 <- medSamp30 - medPop

# average and std dev of sample errors (n=5)
# avg = -68.99081, std dev = 71.8316
meanSamErrb5 <- mean(sampleErrorsb5)
sdSamErrb5 <- sd(sampleErrorsb5)

# average and std dev of sample errors (n=30)
# avg = 18.66129, std dev = 7.281544
meanSamErrb30 <- mean(sampleErrorsb30)
sdSamErrb30 <- sd(sampleErrorsb30)

par(mfrow = c(1,2))
hist(sampleErrorsb5, prob = TRUE, main = "All possible sample errors (n=5) \n attribute of interest is the median", 
     cex.main=0.5, xlim=c(-200, 100), xlab = "Sample Error", col = adjustcolor("purple", 0.3))
abline(v = 0, col = "red")
legend("topleft", c("Avg: -68.99081", "Std Dev: 71.8316"), bty = "n", cex=0.5)
hist(sampleErrorsb30, prob = TRUE, main = "All possible sample errors (n=30) \n attribute of interest is the median", 
     cex.main=0.5, xlim=c(-200, 100), xlab = "Sample Error")
abline(v = 0, col = "red")
legend("topleft", c("Avg: 18.66129", "Std Dev: 7.281544"), bty = "n", cex=0.4)


# e)

# n = 5
sdSamp5 <- apply(samples5, MARGIN = 2, FUN = function(s) {
  sd(samples5[s]) 
})

# n = 30
sdSamp30 <- apply(samples30, MARGIN = 2, FUN = function(s) {
  sd(samples30[s]) 
})

# histogram plot
sdPop <- sd(pop)

sampleErrorsc5 <- sdSamp5 - sdPop
sampleErrorsc30 <- sdSamp30 - sdPop

# average and std dev of sample errors (n=5)
# avg = 13.62639, std dev = 23.51835
meanSamErrc5 <- mean(sampleErrorsc5)
sdSamErrc5 <- sd(sampleErrorsc5)

# average and std dev of sample errors (n=30)
# avg = 8.515346, std dev = 3.76611
meanSamErrc30 <- mean(sampleErrorsc30)
sdSamErrc30 <- sd(sampleErrorsc30)

par(mfrow = c(1,2))
hist(sampleErrorsc5, prob = TRUE, main = "All possible sample errors (n=5) \n attribute of interest is the std dev", 
     cex.main=0.5, xlim=c(-175, 100), xlab = "Sample Error", col = adjustcolor("purple", 0.3))
abline(v = 0, col = "red")
legend("topleft", c("Avg: 13.62639", "Std Dev: 23.51835"), bty = "n", cex=0.4)
hist(sampleErrorsc30, prob = TRUE, main = "All possible sample errors (n=30) \n attribute of interest is the std dev", 
     cex.main=0.5, xlim=c(-175, 100), xlab = "Sample Error")
abline(v = 0, col = "red")
legend("topleft", c("Avg: 8.515346", "Std Dev: 3.76611"), bty = "n", cex=0.4)


# f)

# n = 5
rangeSamp5 <- apply(samples5, MARGIN = 2, FUN = function(s) {
  max(samples5[s]) - min(samples5[s]) 
})

# n = 30
rangeSamp30 <- apply(samples30, MARGIN = 2, FUN = function(s) {
  max(samples30[s]) - min(samples30[s])
})

# histogram plot
rangePop <- max(pop) - min(pop)

sampleErrorsd5 <- rangeSamp5 - rangePop
sampleErrorsd30 <- rangeSamp30 - rangePop

# average and std dev of sample errors (n=5)
# avg = -139.189, std dev = 50.58669
meanSamErrd5 <- mean(sampleErrorsd5)
sdSamErrd5 <- sd(sampleErrorsd5)

# average and std dev of sample errors (n=30)
# avg = -1.18833, std dev = 8.929258
meanSamErrd30 <- mean(sampleErrorsd30)
sdSamErrd30 <- sd(sampleErrorsd30)

par(mfrow = c(1,2))
hist(sampleErrorsd5, prob = TRUE, main = "All possible sample errors (n=5) \n attribute of interest is the range", 
     cex.main=0.5, xlim=c(-400, 50), xlab = "Sample Error", col = adjustcolor("purple", 0.3))
abline(v = 0, col = "red")
legend("topleft", c("Avg: -139.189", "Std Dev: 50.58669"), bty = "n", cex=0.4)
hist(sampleErrorsd30, prob = TRUE, main = "All possible sample errors (n=30) \n attribute of interest is the range", 
     cex.main=0.5, xlim=c(-400, 50), xlab = "Sample Error")
abline(v = 0, col = "red")
legend("topleft", c("Avg: -1.18833", "Std Dev: 8.929258"), bty = "n", cex=0.4)

