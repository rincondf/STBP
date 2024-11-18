# Case 2: Testing dynamic population sizes through group sequential sampling

#######
#T-SPRT
#######

# Endemic and Outbreak trajectories (Table S1) (from Pedigo and Schaik 1984)

m0 <- c(2, 3, 4, 7, 8, 6, 3, 2, 1)
m1 <- c(4, 5, 16, 18, 23, 38, 34, 26, 25)

# upper stop threshold (Eq. 6a in the text)

upper <- function(DDs, m1, m0, ns) {
  log((1 - 0.01)/(0.01)) + (1.16 *
                              cumsum(ns * log((1.16 + m1) / 
                                                ((1.16 + m0)))))
}

# lower stop threshold (Eq. 6b in the text)

lower <- function(DDs, m1, m0, ns) {
  log((0.01)/(1 - 0.01)) + (1.16 *
                              cumsum(ns * log((1.16 + m1) / 
                                                ((1.16 + m0)))))
}

# count weights (Eq. 7 in the text)

calc_weight <- function(m1, m0) {
  log(m1 / m0) - log((1.16 + m1) / (1.16 + m0))
}

# testing trajectories (Eq. 8 in the text)

test_traj <- function(s) {
  m0^(1-s) * m1^(s)
}

# Generation of testing trajectories from a NB distribution

require(matrixStats)

# produce_obs <- function(ns, s) {
#   samD <- matrix(NA, ns, 9)
#   mu <- test_traj(s)
#   for(i in 1: 9) {
#     samD[, i] <- sample(rnbinom(10000, size = 1.16, 
#                                 mu = mu[i]), 10)
#   }
#   return(list(regular = samD, cumulative = rowCumsums(samD)))
# }


# procedure to simulate the T-SPRT

simu_SPRT <- function(s, ns) {
  produce_obs <- function(ns, s) {
    samD <- matrix(NA, ns, 9)
    mu <- test_traj(s)
    for(i in 1: 9) {
      samD[, i] <- sample(rnbinom(10000, size = 1.16, 
                                  mu = mu[i]), ns)
    }
    return(list(regular = samD, cumulative = rowCumsums(samD)))
  }

  test1 <- produce_obs(ns = ns, s = s)$regular
  upper_criterion <- upper(seq(1, 9), m1 = m1, m0 = m0, ns = ns)
  lower_criterion <- lower(seq(1, 9), m1 = m1, m0 = m0, ns = ns)
  weight_count <- cumsum(calc_weight(m1, m0) * colSums(test1))

  # find index where weight count is greater than the upper criterion
  above <- which(weight_count > upper_criterion)
  # find index where weight count is less than the upper criterion
  below <- which(weight_count < lower_criterion)

  # If there is no value above, but there is one below, response is 1
  # set len to the index where weight count became less than the lower criterion
  if(is.na(above[1]) && !is.na(below[1])) {
    resp <- 1
    len <- below[1]
  }

  # If there is a value above, but there is no value below, response is 0
  # set len to the index where weight count became greater than the upper criterion
  if(!is.na(above[1]) && is.na(below[1])) {
    resp <- 0
    len <- above[1]
  }

  # If there are no values either above or below, response is 0
  # set len to the total number of bouts, 9
  if(is.na(above[1]) && is.na(below[1])) {
    resp <- 0
    len <- 9
  }

  # If there are values both above and below...
  if(!is.na(above[1]) && !is.na(below[1])) {
    # Select whichever index the weight count crossed a criteria first
    len <- c(above[1], below[1]) |> min()
    # The response is 1 if there are less values above than below, otherwise the response is 0
    if(length(above) < length(below)) resp <- 1 else resp <- 0
  }

  return(list(result  = resp, bouts = len))
}


#####################################################
# Sequential test of Bayesian posterior probabilities 
#####################################################

# Sequential test of Bayesian posterior probabilities (Eq. 9 in the text)

CP.clo <- function(dat, prd, prior) {
  lld1 <- function(x) {
    prod(dnbinom(dat, mu = x, size = 1.16))
  }
  
  Pr3A <- (integrate(Vectorize(lld1), lower = prd, upper = Inf)$value * prior) /
    (integrate(Vectorize(lld1), lower = 0, upper = prd)$value * (1 - prior) + (integrate(Vectorize(lld1), lower = prd, upper = Inf)$value * prior))
  Pr3A
}

# Procedure to simulate Sequential test of Bayesian posterior probabilities 

simu_SCPTA <- function(s, ns, prior1 = 0.5) {
  m0 <- c(2, 3, 4, 7, 8, 6, 3, 2, 1)
  m1 <- c(4, 5, 16, 18, 23, 38, 34, 26, 25)
  
  test_traj <- function(s) {
    m0^(1-s) * m1^(s)
  }
  
  mid <- test_traj(0.5)
  
  pord_obs <- function(ns, s) {
    
    samD <- matrix(NA, ns, 9)
    mi <- test_traj(s)
    for(i in 1: 9) {
      samD[, i] <- sample(rnbinom(10000, size = 1.16, 
                                  mu = mi[i]), ns)
    }
    return(list(regular = samD, cumulative = rowCumsums(samD)))
  }
  
  test1 <- pord_obs(ns = ns, s = s)$regular
  
  pr2W <- rep(NA, 9)
  
  pr2W[1] <- CP.clo(dat = test1[, 1], prd = mid[1], prior = prior1)
  for(i in 1: 8) {
    if(pr2W[i] < 0.0001) pr2W[i] <- 0.0001
    if(pr2W[i] > 0.9999) pr2W[i] <- 0.9999
    pr2W[i + 1] <- CP.clo(dat = test1[, i + 1], prd = mid[i + 1], prior = pr2W[i])
  }
  
  abo <- which(pr2W > 0.9)
  belo <- which(pr2W < 0.1)
  
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
    len <- 9
  }
  
  if(!is.na(abo[1]) & !is.na(belo[1])) {
    len <- c(abo[1], belo[1])[which.min(c(abo[1], belo[1]))]
    if(length(abo) < length(belo)) resp <- 1 else resp <- 0
  }
  
  return(list(result  = resp, bouts = len, prs = pr2W))
}


#############
# Simulations
#############

# Decisions

repl_SPRT <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- rep(NA, 10)
  
  
  for(i in 1: 10) {
    result[i] <- (sum(replicate(1000, simu_SPRT(levels[i], ns)$result)) / 1000)
  }
  result
}

result5 <- repl_SPRT(5)
result10 <- repl_SPRT(10)
result20 <- repl_SPRT(20)
result30 <- repl_SPRT(30)


repl_SCPTA1 <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- rep(NA, 10)
  priors <- levels
  
  for(i in 1: 10) {
    result[i] <- (sum(replicate(1000, simu_SCPTA(levels[i], ns, prior1 = priors[i])$result)) / 1000)
  }
  result
}


result5CPA1 <- repl_SCPTA1(5)
result10CPA1 <- repl_SCPTA1(10)
result20CPA1 <- repl_SCPTA1(20)
result30CPA1 <- repl_SCPTA1(30)


repl_SCPTA <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- rep(NA, 10)
  
  
  for(i in 1: 10) {
    result[i] <- (sum(replicate(1000, simu_SCPTA(levels[i], ns)$result)) / 1000)
  }
  result
}


result5CPA <- repl_SCPTA(5)
result10CPA <- repl_SCPTA(10)
result20CPA <- repl_SCPTA(20)
result30CPA <- repl_SCPTA(30)



repl_SCPTA2 <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- rep(NA, 10)
  priors <- 1-levels
  
  for(i in 1: 10) {
    result[i] <- (sum(replicate(1000, simu_SCPTA(levels[i], ns, prior1 = priors[i])$result)) / 1000)
  }
  result
}


result5CPA2 <- repl_SCPTA2(5)
result10CPA2 <- repl_SCPTA2(10)
result20CPA2 <- repl_SCPTA2(20)
result30CPA2 <- repl_SCPTA2(30)


# Sample size

repl_SPRTs <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- rep(NA, 10)
  
  
  for(i in 1: 10) {
    result[i] <- (sum(replicate(1000, simu_SPRT(levels[i], ns)$bouts)) / 1000)
  }
  result
}

result5s <- repl_SPRTs(5)
result10s <- repl_SPRTs(10)
result20s <- repl_SPRTs(20)
result30s <- repl_SPRTs(30)


repl_SCPTAs1 <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- rep(NA, 10)
  priors <- levels
  
  for(i in 1: 10) {
    result[i] <- (sum(replicate(1000, simu_SCPTA(levels[i], ns, prior1 = priors[i])$bouts)) / 1000)
  }
  result
}

result5CPAs1 <- repl_SCPTAs1(5)
result10CPAs1 <- repl_SCPTAs1(10)
result20CPAs1 <- repl_SCPTAs1(20)
result30CPAs1 <- repl_SCPTAs1(30)



repl_SCPTAs <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- rep(NA, 10)
  
  
  for(i in 1: 10) {
    result[i] <- (sum(replicate(1000, simu_SCPTA(levels[i], ns)$bouts)) / 1000)
  }
  result
}

result5CPAs <- repl_SCPTAs(5)
result10CPAs <- repl_SCPTAs(10)
result20CPAs <- repl_SCPTAs(20)
result30CPAs <- repl_SCPTAs(30)



repl_SCPTAs2 <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- rep(NA, 10)
  priors <- 1-levels
  
  for(i in 1: 10) {
    result[i] <- (sum(replicate(1000, simu_SCPTA(levels[i], ns, prior1 = priors[i])$bouts)) / 1000)
  }
  result
}

result5CPAs2 <- repl_SCPTAs2(5)
result10CPAs2 <- repl_SCPTAs2(10)
result20CPAs2 <- repl_SCPTAs2(20)
result30CPAs2 <- repl_SCPTAs2(30)

correct2 <- c(rep(1, 5), rep(0, 5))
