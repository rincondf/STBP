# Case 3: Detecting rare species through monitoring

source("STBP.R")

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
calc_posterior_case3 <- function(data,
                           hypothesis,
                           likelihood_func,
                           prior,
                           upper_bnd = Inf) {
  likelihood <- function(x) {
    prod(likelihood_func(data, x))
  }

  H1 <- prior *
            likelihood(hypothesis)
  H0 <- (1 - prior) *
            integrate(
               Vectorize(likelihood),
               lower = hypothesis,
               upper = Inf
            )$value
  posterior <- H1 / (H0 + H1)
  posterior
}

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
  test <- stbp(samples,
               0,
               likelihood_func,
               prior = prior,
               lower_criterion = 0.0001,
               upper_criterion = 0.9999)
  posteriors = test$probabilities
  len = test$decision_indices$first
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

means_det <- c(0.01, 0.05, 0.1, 0.15, 0.2) # tested means

# Bayesian posterior probabilities

# Type II error

beta1 <- rep(NA, 5)

for(i in 1: 5) {
  beta1[i] <- (1000 |>
                replicate(STBP_case3(s = means_det[i], ns = 1)$result) > 0) |>
                which() |>
                length() / 1000
}


beta3 <- rep(NA, 5)

for(i in 1: 5) {
  beta3[i] <- (1000 |>
                replicate(STBP_case3(s = means_det[i], ns = 3)$result) > 0) |>
                which() |>
                length() / 1000
}


beta5 <- rep(NA, 5)

for(i in 1: 5) {
  beta5[i] <- (1000 |>
                replicate(STBP_case3(s = means_det[i], ns = 5)$result) > 0) |>
                which() |>
                length() / 1000
}

beta10 <- rep(NA, 5)

for(i in 1: 5) {
  beta10[i] <- (1000 |>
                replicate(STBP_case3(s = means_det[i], ns = 10)$result) > 0) |>
                which() |>
                length() / 1000
}

# Sample size

size1 <- rep(NA, 5)

for(i in 1: 5) {
  size1[i] <- replicate(1000, STBP_case3(s = means_det[i], ns = 1)$bouts) |>
              mean()
}

size3 <- rep(NA, 5)

for(i in 1: 5) {
  size3[i] <- replicate(1000, STBP_case3(s = means_det[i], ns = 3)$bouts) |>
              mean()
}


size5 <- rep(NA, 5)

for(i in 1: 5) {
  size5[i] <- replicate(1000, STBP_case3(s = means_det[i], ns = 5)$bouts) |>
              mean()
}

size10 <- rep(NA, 5)

for(i in 1: 5) {
  size10[i] <- replicate(1000, STBP_case3(s = means_det[i], ns = 10)$bouts) |>
               mean()
}


# Fixed-sample-size approach

betaGreen30 <- rep(NA, 5)

for(i in 1: 5) {
  betaGreen30[i] <- (1000 |>
                     replicate(rpois(lambda = means_det[i], n = 100000) |>
                               sample(size = 30, replace = FALSE) |>
                               sum()) > 0
                              ) |>
                     which() |>
                     length() / 1000
}

betaGreen20 <- rep(NA, 5)

for(i in 1: 5) {
  betaGreen20[i] <- (1000 |>
                     replicate(rpois(lambda = means_det[i], n = 100000) |>
                               sample(size = 20, replace = FALSE) |>
                               sum()) > 0
                              ) |>
                     which() |>
                     length() / 1000
}

betaGreen10 <- rep(NA, 5)

for(i in 1: 5) {
  betaGreen10[i] <- (1000 |>
                     replicate(rpois(lambda = means_det[i], n = 100000)|>
                               sample(size = 10, replace = FALSE) |>
                               sum()) > 0
                              ) |>
                     which() |>
                     length() / 1000
}


size1_0 <- replicate(1000, STBP_case3(s = 0, ns = 1)$bouts) |>
           mean()
size3_0 <- replicate(1000, STBP_case3(s = 0, ns = 3)$bouts) |>
           mean()
size5_0 <- replicate(1000, STBP_case3(s = 0, ns = 5)$bouts) |>
           mean()
size10_0 <- replicate(1000, STBP_case3(s = 0, ns = 10)$bouts) |>
            mean()


s1alt <- STBP_case3(s = 0, ns = 1, prior = 0.1)$bouts
s3alt <- STBP_case3(s = 0, ns = 3, prior = 0.1)$bouts
s5alt <- STBP_case3(s = 0, ns = 5, prior = 0.1)$bouts
s10alt <- STBP_case3(s = 0, ns = 10, prior = 0.1)$bouts

s1a <- STBP_case3(s = 0, ns = 1, prior = 0.9)$bouts
s3a <- STBP_case3(s = 0, ns = 3, prior = 0.9)$bouts
s5a <- STBP_case3(s = 0, ns = 5, prior = 0.9)$bouts
s10a <- STBP_case3(s = 0, ns = 10, prior = 0.9)$bouts
