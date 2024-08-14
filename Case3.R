# Case 3: Detecting rare species through monitoring

###########################
#Fixed-sample-size approach
###########################

# Type II error (Eq. 10b in the text and used for figure 3)

beta_fun <- function(m, n) {
  exp(-n*m)
}

#####################################################
# Sequential test of Bayesian posterior probabilities 
#####################################################

# Sequential test of Bayesian posterior probabilities (Eq. 11 in the text)

CP.v3C <- function(dat, prd, prior) {
  lld1 <- function(x) {
    prod(dpois(dat, lambda = x))
  }
  
  Pr3A <- (lld1(prd) * prior) /
    (integrate(Vectorize(lld1), lower = prd, upper = Inf)$value * (1 - prior) + (lld1(prd) * prior))
  Pr3A
}

# Procedure to simulate Bayesian posterior probabilities

simu_SCPT_C <- function(s, ns, prior = 0.5) {
  
  pob <- rpois(100000, lambda = s)
  
  pord_obs <- function(ns, s) {
    samD <- matrix(NA, ns, 20)
    for(i in 1: 20) {
      samD[, i] <- sample(pob, ns, replace = FALSE)
    }
    return(list(regular = samD))
  }
  
  test1 <- pord_obs(ns = ns, s = s)$regular
  
  pr2W <- rep(NA, 20)
  
  pr2W[1] <- CP.v3C(dat = test1[, 1], prd = 0, prior = prior)
  for(i in 1: 19) {
    if(pr2W[i] < 0.0001) pr2W[i] <- 0.0001
    if(pr2W[i] > 0.9999) pr2W[i] <- 0.9999
    pr2W[i + 1] <- CP.v3C(dat = test1[, i + 1], prd = 0, prior = pr2W[i])
  }
  
  abo <- which(pr2W > 0.999)
  belo <- which(pr2W < 0.001)
  
  if(is.na(abo[1]) & !is.na(belo[1])) {
    resp <- 1
    len <- belo[1]
  }
  
  if(!is.na(abo[1]) & is.na(belo[1])) {
    resp <- 0
    len <- abo[1]
  }
  
  if(is.na(abo[1]) & is.na(belo[1])) {
    resp <- 0
    len <- 20
  }
  
  if(!is.na(abo[1]) & !is.na(belo[1])) {
    len <- c(abo[1], belo[1])[which.min(c(abo[1], belo[1]))]
    if(length(abo) < length(belo)) resp <- 1 else resp <- 0
  }
  
  return(list(result  = resp, bouts = len, pr = pr2W, data = test1, track = abo))
}


#############
# Simulations
#############

means_det <- c(0.01, 0.05, 0.1, 0.15, 0.2) # tested means

# Bayesian posterior probabilities

# Type II error

beta1 <- rep(NA, 5)

for(i in 1: 5) {
  beta1[i] <- length(which(replicate(1000, simu_SCPT_C(s = means_det[i], ns = 1)$result) > 0) )/1000
}


beta3 <- rep(NA, 5)

for(i in 1: 5) {
  beta3[i] <- length(which(replicate(1000, simu_SCPT_C(s = means_det[i], ns = 3)$result) > 0) )/1000
}


beta5 <- rep(NA, 5)

for(i in 1: 5) {
  beta5[i] <- length(which(replicate(1000, simu_SCPT_C(s = means_det[i], ns = 5)$result) > 0) )/1000
}

beta10 <- rep(NA, 5)

for(i in 1: 5) {
  beta10[i] <- length(which(replicate(1000, simu_SCPT_C(s = means_det[i], ns = 10)$result) > 0) )/1000
}

# Sample size

size1 <- rep(NA, 5)

for(i in 1: 5) {
  size1[i] <- mean(replicate(1000, simu_SCPT_C(s = means_det[i], ns = 1)$bouts))
}

size3 <- rep(NA, 5)

for(i in 1: 5) {
  size3[i] <- mean(replicate(1000, simu_SCPT_C(s = means_det[i], ns = 3)$bouts))
}


size5 <- rep(NA, 5)

for(i in 1: 5) {
  size5[i] <- mean(replicate(1000, simu_SCPT_C(s = means_det[i], ns = 5)$bouts))
}

size10 <- rep(NA, 5)

for(i in 1: 5) {
  size10[i] <- mean(replicate(1000, simu_SCPT_C(s = means_det[i], ns = 10)$bouts))
}


# Fixed-sample-size approach

betaGreen30 <- rep(NA, 5)

for(i in 1: 5) {
  betaGreen30[i] <- length(which(replicate(1000, sum(sample(rpois(lambda = means_det[i], n = 100000), size = 30, replace = F))) > 0))/1000
}

betaGreen20 <- rep(NA, 5)

for(i in 1: 5) {
  betaGreen20[i] <- length(which(replicate(1000, sum(sample(rpois(lambda = means_det[i], n = 100000), size = 20, replace = F))) > 0))/1000
}

betaGreen10 <- rep(NA, 5)

for(i in 1: 5) {
  betaGreen10[i] <- length(which(replicate(1000, sum(sample(rpois(lambda = means_det[i], n = 100000), size = 10, replace = F))) > 0))/1000
}


size1_0 <- mean(replicate(1000, simu_SCPT_C(s = 0, ns = 1)$bouts))
size3_0 <- mean(replicate(1000, simu_SCPT_C(s = 0, ns = 3)$bouts))
size5_0 <- mean(replicate(1000, simu_SCPT_C(s = 0, ns = 5)$bouts))
size10_0 <- mean(replicate(1000, simu_SCPT_C(s = 0, ns = 10)$bouts))


s1alt <- simu_SCPT_C(s = 0, ns = 1, prior = 0.1)$bouts
s3alt <- simu_SCPT_C(s = 0, ns = 3, prior = 0.1)$bouts
s5alt <- simu_SCPT_C(s = 0, ns = 5, prior = 0.1)$bouts
s10alt <- simu_SCPT_C(s = 0, ns = 10, prior = 0.1)$bouts

s1a <- simu_SCPT_C(s = 0, ns = 1, prior = 0.9)$bouts
s3a <- simu_SCPT_C(s = 0, ns = 3, prior = 0.9)$bouts
s5a <- simu_SCPT_C(s = 0, ns = 5, prior = 0.9)$bouts
s10a <- simu_SCPT_C(s = 0, ns = 10, prior = 0.9)$bouts
