# Case 2: Testing dynamic population sizes through group sequential sampling

#######
#T-SPRT
#######

# Endemic and Outbreak trajectories (Table S1) (from Pedigo and Schaik 1984)
m0 <- c(2, 3, 4, 7, 8, 6, 3, 2, 1)
m1 <- c(4, 5, 16, 18, 23, 38, 34, 26, 25)

# upper stop threshold (Eq. 6a in the text)
upper <- function(DDs, m1, m0, ns) {
  log((1 - 0.01)/(0.01)) +
  (1.16 * cumsum(ns * log((1.16 + m1) / ((1.16 + m0)))))
}

# lower stop threshold (Eq. 6b in the text)
lower <- function(DDs, m1, m0, ns) {
  log((0.01)/(1 - 0.01)) +
  (1.16 * cumsum(ns * log((1.16 + m1) / ((1.16 + m0)))))
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
produce_obs <- function(ns, s) {
  samD <- matrix(NA, ns, 9)
  mu <- test_traj(s)
  for(i in 1: 9) {
    samD[, i] <- sample(rnbinom(10000, size = 1.16,
                                mu = mu[i]), ns, replace = FALSE)
  }
  return(list(regular = samD, cumulative = rowCumsums(samD)))
}

require(matrixStats)
source("STBP.R")

# Simulation of the T-SPRT

simu_SPRT <- function(s, ns) {

  samples <- produce_obs(ns = ns, s = s)$regular
  upper_criterion <- upper(seq(1, 9), m1 = m1, m0 = m0, ns = ns)
  lower_criterion <- lower(seq(1, 9), m1 = m1, m0 = m0, ns = ns)
  weight_count <- cumsum(calc_weight(m1, m0) * colSums(samples))

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
  # set len to the index where weight count
  # became greater than the upper criterion
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
    # The response is 1 if the first value is above,
    # otherwise the response is 0
    if(above[1] < below[1]) resp <- 1 else resp <- 0
  }

  return(list(result  = resp, bouts = len))
}


#####################################################
# Sequential test of Bayesian posterior probabilities
#####################################################

# Procedure to simulate Sequential test of Bayesian posterior probabilities

STBP_case2 <- function(s, ns, prior1 = 0.5) {
  trajectory <- test_traj(0.5)
  samples <- produce_obs(ns = ns, s = s)$regular

  likelihood_func <- function(data, mu) {
    dnbinom(data, mu = mu, size = 1.16)
  }

  test <- stbp(samples,
               trajectory,
               likelihood_func,
               prior = prior1,
               lower_criterion = 0.05,
               upper_criterion = 0.95)
  posteriors = test$probabilities
  len = test$num_iterations
  response = test$recommendation

  return(list(result  = response, bouts = len, prs = posteriors))
}


#############
# Simulations
#############

# For the sake of efficiency, this code runs simulations with future’s parallel 
# processing capabilities using the package furrr.

# Simulations can also be run with conventional, sequential processing. 
# Code provided in Seq_simulations.R

require(furrr)
set.seed(123)

ncores <- 13 # Set the number of available cores

plan(multisession, workers = ncores)

# Decisions

repl_SPRT <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- levels |> future_map_dbl(
              ~simu_SPRT(., ns)$result |>
                replicate(n = 1000) |>
                mean(),
              .options = furrr_options(seed = 123)
            )
  result
}

result5 <- repl_SPRT(5)
result10 <- repl_SPRT(10)
result20 <- repl_SPRT(20)
result30 <- repl_SPRT(30)


repl_SCPTA1 <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- levels |> future_map_dbl(
              ~STBP_case2(., ns, prior1 = .)$result |>
                replicate(n = 1000) |>
                mean(),
              .options = furrr_options(seed = 123)
            )
  result
}


result5CPA1 <- repl_SCPTA1(5)
result10CPA1 <- repl_SCPTA1(10)
result20CPA1 <- repl_SCPTA1(20)
result30CPA1 <- repl_SCPTA1(30)


repl_SCPTA <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- levels |> future_map_dbl(
              ~STBP_case2(., ns)$result |>
                replicate(n = 1000) |>
                mean(),
              .options = furrr_options(seed = 123)
            )
  result
}


result5CPA <- repl_SCPTA(5)
result10CPA <- repl_SCPTA(10)
result20CPA <- repl_SCPTA(20)
result30CPA <- repl_SCPTA(30)



repl_SCPTA2 <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- levels |> future_map_dbl(
              ~STBP_case2(., ns, prior1 = 1 - .)$result |>
                replicate(n = 1000) |>
                mean(),
              .options = furrr_options(seed = 123)
            )
  result
}


result5CPA2 <- repl_SCPTA2(5)
result10CPA2 <- repl_SCPTA2(10)
result20CPA2 <- repl_SCPTA2(20)
result30CPA2 <- repl_SCPTA2(30)


# Sample size

repl_SPRTs <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- levels |> future_map_dbl(
              ~simu_SPRT(., ns)$bouts |>
                replicate(n = 1000) |>
                mean(),
              .options = furrr_options(seed = 123)
            )
  result
}

result5s <- repl_SPRTs(5)
result10s <- repl_SPRTs(10)
result20s <- repl_SPRTs(20)
result30s <- repl_SPRTs(30)


repl_SCPTAs1 <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- levels |> future_map_dbl(
              ~STBP_case2(., ns, prior1 = .)$bouts |>
                replicate(n = 1000) |>
                mean(),
              .options = furrr_options(seed = 123)
            )
  result
  result
}

result5CPAs1 <- repl_SCPTAs1(5)
result10CPAs1 <- repl_SCPTAs1(10)
result20CPAs1 <- repl_SCPTAs1(20)
result30CPAs1 <- repl_SCPTAs1(30)



repl_SCPTAs <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- levels |> future_map_dbl(
              ~STBP_case2(., ns)$bouts |>
                replicate(n = 1000) |>
                mean(),
              .options = furrr_options(seed = 123)
            )
  result
}

result5CPAs <- repl_SCPTAs(5)
result10CPAs <- repl_SCPTAs(10)
result20CPAs <- repl_SCPTAs(20)
result30CPAs <- repl_SCPTAs(30)



repl_SCPTAs2 <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- levels |> future_map_dbl(
              ~STBP_case2(., ns, prior1 = 1 - .)$bouts |>
                replicate(n = 1000) |>
                mean(),
              .options = furrr_options(seed = 123)
            )
  result
  result
}

result5CPAs2 <- repl_SCPTAs2(5)
result10CPAs2 <- repl_SCPTAs2(10)
result20CPAs2 <- repl_SCPTAs2(20)
result30CPAs2 <- repl_SCPTAs2(30)

correct2 <- c(rep(1, 5), rep(0, 5))

plan(sequential) # back to sequential computing (housekeeping)

#########
# Metrics
#########

# Overall error rate

# T-SPRT
mean(c(mean(1 - (1 - abs(correct2 - result5))),
       mean(1 - (1 - abs(correct2 - result10))),
       mean(1 - (1 - abs(correct2 - result20))),
       mean(1 - (1 - abs(correct2 - result30)))))

# STBP with naive init priors
mean(c(mean(1 - (1 - abs(correct2 - result5CPA))),
       mean(1 - (1 - abs(correct2 - result10CPA))),
       mean(1 - (1 - abs(correct2 - result20CPA))),
       mean(1 - (1 - abs(correct2 - result30CPA)))))

# STBP with correct init priors
mean(c(mean(1 - (1 - abs(correct2 - result5CPA1))),
       mean(1 - (1 - abs(correct2 - result10CPA1))),
       mean(1 - (1 - abs(correct2 - result20CPA1))),
       mean(1 - (1 - abs(correct2 - result30CPA1)))))

# STBP with incorrect init priors
mean(c(mean(1 - (1 - abs(correct2 - result5CPA2))),
       mean(1 - (1 - abs(correct2 - result10CPA2))),
       mean(1 - (1 - abs(correct2 - result20CPA2))),
       mean(1 - (1 - abs(correct2 - result30CPA2)))))


# Overall error rate for STPB excluding n = 5 

# STBP with naive init priors
mean(c(
  mean(1 - (1 - abs(correct2 - result10CPA))),
  mean(1 - (1 - abs(correct2 - result20CPA))),
  mean(1 - (1 - abs(correct2 - result30CPA)))))

# STBP with correct init priors
mean(c(
  mean(1 - (1 - abs(correct2 - result10CPA1))),
  mean(1 - (1 - abs(correct2 - result20CPA1))),
  mean(1 - (1 - abs(correct2 - result30CPA1)))))

# STBP with incorrect init priors
mean(c(
  mean(1 - (1 - abs(correct2 - result10CPA2))),
  mean(1 - (1 - abs(correct2 - result20CPA2))),
  mean(1 - (1 - abs(correct2 - result30CPA2)))))



# Type I error

# T-SPRT
mean(c(mean(1 - (1 - abs(correct2[1:4] - result5[1:4]))),
       mean(1 - (1 - abs(correct2[1:4] - result10[1:4]))),
       mean(1 - (1 - abs(correct2[1:4] - result20[1:4]))),
       mean(1 - (1 - abs(correct2[1:4] - result30[1:4])))))

# STBP with naive init priors
mean(c(mean(1 - (1 - abs(correct2[1:4] - result5CPA[1:4]))),
       mean(1 - (1 - abs(correct2[1:4] - result10CPA[1:4]))),
       mean(1 - (1 - abs(correct2[1:4] - result20CPA[1:4]))),
       mean(1 - (1 - abs(correct2[1:4] - result30CPA[1:4])))))

# STBP with correct init priors
mean(c(mean(1 - (1 - abs(correct2[1:4] - result5CPA1[1:4]))),
       mean(1 - (1 - abs(correct2[1:4] - result10CPA1[1:4]))),
       mean(1 - (1 - abs(correct2[1:4] - result20CPA1[1:4]))),
       mean(1 - (1 - abs(correct2[1:4] - result30CPA1[1:4])))))

# STBP with incorrect init priors
mean(c(mean(1 - (1 - abs(correct2[1:4] - result5CPA2[1:4]))),
       mean(1 - (1 - abs(correct2[1:4] - result10CPA2[1:4]))),
       mean(1 - (1 - abs(correct2[1:4] - result20CPA2[1:4]))),
       mean(1 - (1 - abs(correct2[1:4] - result30CPA2[1:4])))))



# Type I error for STPB excluding n = 5

# STBP with naive init priors
mean(c(
  mean(1 - (1 - abs(correct2[1:4] - result10CPA[1:4]))),
  mean(1 - (1 - abs(correct2[1:4] - result20CPA[1:4]))),
  mean(1 - (1 - abs(correct2[1:4] - result30CPA[1:4])))))

# STBP with correct init priors
mean(c(
  mean(1 - (1 - abs(correct2[1:4] - result10CPA1[1:4]))),
  mean(1 - (1 - abs(correct2[1:4] - result20CPA1[1:4]))),
  mean(1 - (1 - abs(correct2[1:4] - result30CPA1[1:4])))))

# STBP with incorrect init priors
mean(c(
  mean(1 - (1 - abs(correct2[1:4] - result10CPA2[1:4]))),
  mean(1 - (1 - abs(correct2[1:4] - result20CPA2[1:4]))),
  mean(1 - (1 - abs(correct2[1:4] - result30CPA2[1:4])))))



# Type II error

# T-SPRT
mean(c(mean(1 - (1 - abs(correct2[5:10] - result5[5:10]))),
       mean(1 - (1 - abs(correct2[5:10] - result10[5:10]))),
       mean(1 - (1 - abs(correct2[5:10] - result20[5:10]))),
       mean(1 - (1 - abs(correct2[5:10] - result30[5:10])))))

# STBP with naive init priors
mean(c(mean(1 - (1 - abs(correct2[5:10] - result5CPA[5:10]))),
       mean(1 - (1 - abs(correct2[5:10] - result10CPA[5:10]))),
       mean(1 - (1 - abs(correct2[5:10] - result20CPA[5:10]))),
       mean(1 - (1 - abs(correct2[5:10] - result30CPA[5:10])))))

# STBP with correct init priors
mean(c(mean(1 - (1 - abs(correct2[5:10] - result5CPA1[5:10]))),
       mean(1 - (1 - abs(correct2[5:10] - result10CPA1[5:10]))),
       mean(1 - (1 - abs(correct2[5:10] - result20CPA1[5:10]))),
       mean(1 - (1 - abs(correct2[5:10] - result30CPA1[5:10])))))

# STBP with incorrect init priors
mean(c(mean(1 - (1 - abs(correct2[5:10] - result5CPA2[5:10]))),
       mean(1 - (1 - abs(correct2[5:10] - result10CPA2[5:10]))),
       mean(1 - (1 - abs(correct2[5:10] - result20CPA2[5:10]))),
       mean(1 - (1 - abs(correct2[5:10] - result30CPA2[5:10])))))

# Overall average number of sampling bouts

# T-SPRT
mean(c(mean(result5s),
       mean(result10s),
       mean(result20s),
       mean(result30s)))

# STBP with naive init priors
mean(c(mean(result5CPAs),
       mean(result10CPAs),
       mean(result20CPAs),
       mean(result30CPAs)))

# STBP with correct init priors
mean(c(mean(result5CPAs1),
       mean(result10CPAs1),
       mean(result20CPAs1),
       mean(result30CPAs1)))

# STBP with incorrect init priors
mean(c(mean(result5CPAs2),
       mean(result10CPAs2),
       mean(result20CPAs2),
       mean(result30CPAs2)))

# Average number of sampling bouts excluding n = 5

# STBP with naive init priors
mean(c(mean(result10CPAs),
       mean(result20CPAs),
       mean(result30CPAs)))

# STBP with correct init priors
mean(c(mean(result10CPAs1),
       mean(result20CPAs1),
       mean(result30CPAs1)))

# STBP with incorrect init priors
mean(c(mean(result10CPAs2),
       mean(result20CPAs2),
       mean(result30CPAs2)))
