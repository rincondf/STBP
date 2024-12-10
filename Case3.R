# Case 3: Detecting rare species through monitoring

source("STBP.R")

###########################
#Fixed-sample-size approach
###########################

# Type II error (Eq. 11b in the text and used for figure 3a)

beta_fun <- function(m, n) {
  exp(-n * m)
}

#####################################################
# Sequential test of Bayesian posterior probabilities
#####################################################

# Procedure to simulate Bayesian posterior probabilities

STBP_case3 <- function(s, ns, prior = 0.5) {
  # generate population pool to sample from
  pool <- rpois(100000, lambda = s)
  # function to sample from the pool
  produce_obs <- function(ns) {
    samD <- matrix(NA, ns, 20)
    for(i in 1: 20) {
      samD[, i] <- sample(pool, ns, replace = FALSE)
    }
    return(list(regular = samD))
  }
  samples <- produce_obs(ns = ns)$regular
  
  likelihood_func <- function(data, lambda) {
    dpois(data, lambda)
  }
  test <- STBP_simpleH(samples,
                       0,
                       likelihood_func,
                       prior = prior,
                       lower_criterion = 0.0001,
                       upper_criterion = 0.9999)
  posteriors = test$probabilities
  len = test$num_iterations
  response = test$recommendation
  above = test$decision_indices$above
  
  return(list(result = response,
              bouts = len,
              pr = posteriors,
              data = samples,
              track = above))
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

means_det <- c(0.01, 0.05, 0.1, 0.15, 0.2) # tested means

# Bayesian posterior probabilities

# Type II error

beta1 <- means_det |> future_map_dbl(
            ~((STBP_case3(s = ., ns = 1)$result) |>
              replicate(n = 1000) > 0) |>
              which() |>
              length() /
              1000,
            .options = furrr_options(seed = 123)
          )


beta3 <- means_det |> future_map_dbl(
            ~((STBP_case3(s = ., ns = 3)$result) |>
              replicate(n = 1000) > 0) |>
              which() |>
              length() /
              1000,
            .options = furrr_options(seed = 123)
          )



beta5 <- means_det |> future_map_dbl(
          ~((STBP_case3(s = ., ns = 5)$result) |>
            replicate(n = 1000) > 0) |>
            which() |>
            length() /
            1000,
            .options = furrr_options(seed = 123)
          )

beta10 <- means_det |> future_map_dbl(
            ~((STBP_case3(s = ., ns = 10)$result) |>
              replicate(n = 1000) > 0) |>
              which() |>
              length() /
              1000,
            .options = furrr_options(seed = 123)
          )

# Sample size

size1 <- means_det |> future_map_dbl(
            ~STBP_case3(s = ., ns = 1)$bouts |>
              replicate(n = 1000) |>
              mean(),
            .options = furrr_options(seed = 123)
          )


size3 <- means_det |> future_map_dbl(
            ~STBP_case3(s = ., ns = 3)$bouts |>
              replicate(n = 1000) |>
              mean(),
            .options = furrr_options(seed = 123)
          )

size5 <- means_det |> future_map_dbl(
          ~STBP_case3(s = ., ns = 5)$bouts |>
            replicate(n = 1000) |>
            mean(),
          .options = furrr_options(seed = 123)
        )

size10 <- means_det |> future_map_dbl(
            ~STBP_case3(s = ., ns = 10)$bouts |>
              replicate(n = 1000) |>
              mean(),
            .options = furrr_options(seed = 123)
          )


# Fixed-sample-size approach

betaGreen30 <- means_det |> future_map_dbl(
                ~(rpois(lambda = ., n = 100000) |>
                    sample(size = 30, replace = FALSE) |>
                    sum() |>
                    replicate(n = 1000) > 0) |>
                  which() |>
                  length() /
                  1000,
                .options = furrr_options(seed = 123)
              )

betaGreen20 <- means_det |> future_map_dbl(
                ~(rpois(lambda = ., n = 100000) |>
                    sample(size = 20, replace = FALSE) |>
                    sum() |>
                    replicate(n = 1000) > 0) |>
                  which() |>
                  length() /
                  1000,
                .options = furrr_options(seed = 123)
              )

betaGreen10 <- means_det |> future_map_dbl(
                ~(rpois(lambda = ., n = 100000) |>
                    sample(size = 10, replace = FALSE) |>
                    sum() |>
                    replicate(n = 1000) > 0) |>
                  which() |>
                  length() /
                  1000,
                .options = furrr_options(seed = 123)
              )

plan(sequential) # back to sequential computing (housekeeping)

# Sample sizes when the species is actually absent

# With naive priors
size1_0 <- STBP_case3(s = 0, ns = 1)$bouts
size3_0 <- STBP_case3(s = 0, ns = 3)$bouts
size5_0 <- STBP_case3(s = 0, ns = 5)$bouts
size10_0 <- STBP_case3(s = 0, ns = 10)$bouts

# With good reasons to expect the species present (low credibility for H1)
s1alt <- STBP_case3(s = 0, ns = 1, prior = 0.1)$bouts
s3alt <- STBP_case3(s = 0, ns = 3, prior = 0.1)$bouts
s5alt <- STBP_case3(s = 0, ns = 5, prior = 0.1)$bouts
s10alt <- STBP_case3(s = 0, ns = 10, prior = 0.1)$bouts

# With good reasons to expect the species absent (high credibility for H1)
s1a <- STBP_case3(s = 0, ns = 1, prior = 0.9)$bouts
s3a <- STBP_case3(s = 0, ns = 3, prior = 0.9)$bouts
s5a <- STBP_case3(s = 0, ns = 5, prior = 0.9)$bouts
s10a <- STBP_case3(s = 0, ns = 10, prior = 0.9)$bouts
