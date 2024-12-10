# Case 1: Testing static population sizes through purely sequential sampling

#####
#SPRT
#####


require(truncdist)
source("STBP.R")

# Function to estimate parameter k for NB distributions 
# (from Rincon et al. 2021)

estimate_k <- function(mean) {
  a = 1.830012
  b = 1.218041
  (mean^2) / ((a * mean^(b)) - mean)
}

# same but for simulation adding an stochastic component the standard
# deviation of the random component is the mean squared error of the mean-
# variance model fit by Rincon et al. (2021)

estimate_k_stoch <- function(mean) {
  a <- 1.830012
  b <- 1.218041
  times <- length(mean)
  a1 <- rep(NA, times)
  for (i in 1:times) {
    a1[i] <- (mean[i]^2) /
      ((a * mean[i]^(b) *
          exp(rtrunc(1, "norm", a = log(1 / (a * mean[i]^(b - 1))),
                     b = Inf, mean = 0, sd = 0.3222354)))
       - mean[i])
  }
  a1
}

# value for k at the threshold

k_9 <- estimate_k(9)

# low intercept for stop line for negative binomial distribution
# (from Binns, Nyrop and Werf, 2000)

low_int_nb <- function(alpha, beta, mu0, mu1, k_est){
  (log(beta / (1 - alpha))) /
  (log((mu1 * (mu0 + k_est)) / (mu0 * (mu1 + k_est))))
}

lower_criterion_intercept <- low_int_nb(alpha = 0.1,
                                        beta = 0.1,
                                        mu0 = 8,
                                        mu1 = 10,
                                        k_est = k_9)

# hi intercept for stop line

hi_int_nb <- function(alpha, beta, mu0, mu1, k_est){
  (log((1 - beta) / (alpha))) / 
  (log((mu1 * (mu0 + k_est)) / (mu0 * (mu1 + k_est))))
}


higher_criterion_intercept <- hi_int_nb(alpha = 0.1,
                                        beta = 0.1,
                                        mu0 = 8,
                                        mu1 = 10,
                                        k_est = k_9)

# slope for both lines

criterion_slope_nb <- function(alpha, beta, mu0, mu1, k_est){
  (k_est * log((mu1 + k_est) / (mu0 + k_est))) /
    (log((mu1 * (mu0 + k_est)) / (mu0 * (mu1 + k_est))))
}

criteria_slope <- criterion_slope_nb(alpha = 0.1,
                                     beta = 0.1,
                                     mu0 = 8,
                                     mu1 = 10,
                                     k_est = k_9)

# Functions for stop lines

low_criterion_line <- function(x){
  criteria_slope * x + lower_criterion_intercept
}

hi_criterion_line <- function(x){
  criteria_slope * x + higher_criterion_intercept
}

# procedure to simulate SPRT

SPRT_case1 <- function(d){
  samples <- rep(NA, 100)
  pool <- rnbinom(mu = d, size = estimate_k_stoch(d), n = 6000)
  mean_Re <- mean(pool)
  for(i in 1:100){
    samples[i] <- sample(pool, size = 1, replace = FALSE)
    pool <- pool[-match(samples[i], pool)]
    dat <- cumsum(samples)[i]
    if ((dat < low_criterion_line(i)) || (dat > hi_criterion_line(i))) break
  }
  if (dat < low_criterion_line(i)) {
    resp <- 1
  } else {
    resp <- 0
  }
  return(list(estimate = dat / i,
              samples = i,
              recommendation = resp,
              col_data = samples,
              mean = mean_Re)
             )
}


#####################################################
# Sequential test of Bayesian posterior probabilities
#####################################################

# Procedure to simulate Sequential test of Bayesian posterior probabilities

STBP_case1 <- function(pop_mean, prior){
  samples <- rep(NA, 100)
  pool <- rnbinom(mu = pop_mean, size = estimate_k_stoch(pop_mean), n = 6000)
  mean_Re <- mean(pool)
  for(i in 1:100){
    samples[i] <- sample(pool, size = 1, replace = FALSE)
    pool <- pool[-match(samples[i], pool)]
  }
  likelihood_func <- function(data, mu) {
    dnbinom(data,
            mu = mu,
            size = if(estimate_k(mu) < 0 | is.nan(estimate_k(mu))) 0
                   else estimate_k(mu)
     )
  }
  test <- stbp(data= samples,
               hypotheses = 9,
               likelihood_func = likelihood_func,
               prior = prior,
               lower_criterion = 0.01,
               upper_criterion = 0.99)
  return(list(
    Probabilities = test$probabilities,
    samples = test$num_iterations,
    recommendation = test$recommendation,
    col_data = samples,
    mean = mean_Re
  ))
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

correct_prior <- function(hypothesis) {
  if (hypothesis == 9) return(0.5)
  ifelse(hypothesis < 9, 0.1, 0.9)
}
incorrect_prior <- function(hypothesis) {
  ifelse(hypothesis >= 9, 0.1, 0.9)
}

# Decisions

SPRTA <- (1:13) |>
          future_map_dbl(
            ~SPRT_case1(.)$recommendation |>
               replicate(n = 1000) |>
               mean(),
            .options = furrr_options(seed = 123)
          )

STCHA <- (1:13) |>
          future_map_dbl(
            ~STBP_case1(., correct_prior(.))$recommendation |>
               replicate(n = 1000) |>
               mean(),
            .options = furrr_options(seed = 123)
          )

STCHAa <- (1:13) |>
          future_map_dbl(
            ~STBP_case1(., 0.5)$recommendation |>
               replicate(n = 1000) |>
               mean(),
            .options = furrr_options(seed = 123)
          )

STCHAb <- (1:13) |>
          future_map_dbl(
            ~STBP_case1(., incorrect_prior(.))$recommendation |>
               replicate(n = 1000) |>
               mean(),
            .options = furrr_options(seed = 123)
          )

correct1 <- c(rep(1, 8), rep(0, 5))


# Sample size


SPRTAs <- (1:13) |>
          future_map_dbl(
            ~SPRT_case1(.)$samples |>
               replicate(n = 1000) |>
               mean(),
            .options = furrr_options(seed = 123)
          )

STCHAs <- (1:13) |>
          future_map_dbl(
            ~STBP_case1(., correct_prior(.))$samples |>
               replicate(n = 1000) |>
               mean(),
            .options = furrr_options(seed = 123)
          )

STCHAas <- (1:13) |>
          future_map_dbl(
            ~STBP_case1(., 0.5)$samples |>
               replicate(n = 1000) |>
               mean(),
            .options = furrr_options(seed = 123)
          )

STCHAbs <- (1:13) |>
          future_map_dbl(
            ~STBP_case1(., incorrect_prior(.))$samples |>
               replicate(n = 1000) |>
               mean(),
            .options = furrr_options(seed = 123)
          )

plan(sequential) # back to sequential computing (housekeeping)

#########
# Metrics
#########

# Overall error rate
sum(1 - (1 - abs(correct1 - SPRTA))) / 13 # for SPRT
sum(1 - (1 - abs(correct1 - STCHA))) / 13 # for STBP with correct init priors
sum(1 - (1 - abs(correct1 - STCHAa))) / 13 # for STBP with naive init priors
sum(1 - (1 - abs(correct1 - STCHAb))) / 13 # for STBP with incorrect init priors

# Type II error
sum(1 - (1 - abs(correct1[9:13] - SPRTA[9:13]))) / 5 # for SPRT
sum(1 - (1 - abs(correct1[9:13] - STCHA[9:13]))) / 5 # for STBP with correct init priors
sum(1 - (1 - abs(correct1[9:13] - STCHAa[9:13]))) / 5 # for STBP with naive init priors
sum(1 - (1 - abs(correct1[9:13] - STCHAb[9:13]))) / 5 # for STBP with incorrect init priors

# Type I error
sum(1 - (1 - abs(correct1[1:8] - SPRTA[1:8]))) / 8 # for SPRT
sum(1 - (1 - abs(correct1[1:8] - STCHA[1:8]))) / 8 # for STBP with correct init priors
sum(1 - (1 - abs(correct1[1:8] - STCHAa[1:8]))) / 8 # for STBP with naive init priors
sum(1 - (1 - abs(correct1[1:8] - STCHAb[1:8]))) / 8 # for STBP with incorrect init priors

# Mean sample sizes required
mean(SPRTAs)
mean(STCHAs)
mean(STCHAas)
mean(STCHAbs)
