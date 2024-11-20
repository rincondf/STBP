# Case 1: Testing static population sizes through purely sequential sampling

#####
#SPRT
#####


require(truncdist)

# function to estimate k parameter from NB distributions (from Rincon et al. 2021)

estimate_k <- function(mean) {
  a = exp(0.6043225)
  b = 1.218041
  (mean^2) / ((a * mean^(b)) - mean)
}

# same but for simulation adding an stochastic component

estimate_k_stoch <- function(mean) {
  a <- exp(0.6043225)
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

# low intercept for stop line for negative binomial distribution (from Binns, Nyrop and Werf, 2000)

low_int_nb <- function(alpha, beta, mu0, mu1, k_est){
  (log(beta / (1 - alpha))) / (log((mu1 * (mu0 + k_est)) / (mu0 * (mu1 + k_est))))
}

lower_criterion_intercept <- low_int_nb(alpha = 0.1, beta = 0.1, mu0 = 8, mu1 = 10, k_est = k_9)

# hi intercept for stop line

hi_int_nb <- function(alpha, beta, mu0, mu1, k_est){
  (log((1 - beta) / (alpha))) / (log((mu1 * (mu0 + k_est)) / (mu0 * (mu1 + k_est))))
}


higher_criterion_intercept <- hi_int_nb(alpha = 0.1, beta = 0.1, mu0 = 8, mu1 = 10, k_est = k_9)

# slope for both lines

criterion_slope_nb <- function(alpha, beta, mu0, mu1, k_est){
  (k_est * log((mu1 + k_est) / (mu0 + k_est))) /
    (log((mu1 * (mu0 + k_est)) / (mu0 * (mu1 + k_est))))
}

criteria_slope <- criterion_slope_nb(alpha = 0.1, beta = 0.1, mu0 = 8, mu1 = 10, k_est = k_9)

# Functions for stop lines

low_criterion_line <- function(x){
  criteria_slope * x + lower_criterion_intercept
}

hi_criterion_line <- function(x){
  criteria_slope * x + higher_criterion_intercept
}

# procedure to simulte SPRT

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

# Sequential test of Bayesian posterior probabilities (Eq. 5 in the text)

calc_posterior <- function(data, hypothesis, prior) {
  likelihood <- function(x) {
    prod(dnbinom(data,
                 mu = x,
                 size = if(estimate_k(x) < 0 || is.nan(estimate_k(x))) 0
                        else estimate_k(x)
                ))
  }

  null <- prior * integrate(
                    Vectorize(likelihood),
                    lower = hypothesis,
                    upper = Inf
                  )$value
  alt <- (1 - prior) * integrate(
                         Vectorize(likelihood),
                         lower = 0,
                         upper = hypothesis
                       )$value
  posterior <- null / (alt + null)
  posterior
}


# procedure to simulate Sequential test of Bayesian posterior probabilities

STBP_case1 <- function(pop_mean, prior){
  samples <- rep(NA, 100)
  pool <- rnbinom(mu = pop_mean, size = estimate_k_stoch(pop_mean), n = 6000)
  mean_Re <- mean(pool)
  posteriors <- c(prior, rep(NA, 99))
  for(i in 1:100){
    samples[i] <- sample(pool, size = 1, replace = FALSE)
    pool <- pool[-match(samples[i], pool)]
    posteriors[i + 1] <- calc_posterior(data = samples[i],
                                        hypothesis = 9,
                                        prior = posteriors[i]
                                        )
    if ((length(!is.na(samples)) > 3) &&
        ((posteriors[i + 1] < 0.01) ||
        (posteriors[i + 1] > 0.99))
       ) break
  }
  if (posteriors[i + 1] < 0.01) {
    resp <- 1
  } else {
    resp <- 0
  }
  return(list(
    Probabilities = posteriors,
    samples = i,
    recommendation = resp,
    col_data = samples,
    mean = mean_Re
  ))
}

#############
# Simulations
#############


# Decisions

SPRTA <- c(sum(replicate(1000, SPRT_case1(1)$recommendation))/1000,
           sum(replicate(1000, SPRT_case1(2)$recommendation))/1000,
           sum(replicate(1000, SPRT_case1(3)$recommendation))/1000,
           sum(replicate(1000, SPRT_case1(4)$recommendation))/1000,
           sum(replicate(1000, SPRT_case1(5)$recommendation))/1000,
           sum(replicate(1000, SPRT_case1(6)$recommendation))/1000,
           sum(replicate(1000, SPRT_case1(7)$recommendation))/1000,
           sum(replicate(1000, SPRT_case1(8)$recommendation))/1000,
           sum(replicate(1000, SPRT_case1(9)$recommendation))/1000,
           sum(replicate(1000, SPRT_case1(10)$recommendation))/1000,
           sum(replicate(1000, SPRT_case1(11)$recommendation))/1000,
           sum(replicate(1000, SPRT_case1(12)$recommendation))/1000,
           sum(replicate(1000, SPRT_case1(13)$recommendation))/1000)

STCHA <- c(sum(replicate(1000, STBP_case1(1, 0.1)$recommendation))/1000,
           sum(replicate(1000, STBP_case1(2, 0.1)$recommendation))/1000,
           sum(replicate(1000, STBP_case1(3, 0.1)$recommendation))/1000,
           sum(replicate(1000, STBP_case1(4, 0.1)$recommendation))/1000,
           sum(replicate(1000, STBP_case1(5, 0.1)$recommendation))/1000,
           sum(replicate(1000, STBP_case1(6, 0.1)$recommendation))/1000,
           sum(replicate(1000, STBP_case1(7, 0.1)$recommendation))/1000,
           sum(replicate(1000, STBP_case1(8, 0.1)$recommendation))/1000,
           sum(replicate(1000, STBP_case1(9, 0.5)$recommendation))/1000,
           sum(replicate(1000, STBP_case1(10, 0.9)$recommendation))/1000,
           sum(replicate(1000, STBP_case1(11, 0.9)$recommendation))/1000,
           sum(replicate(1000, STBP_case1(12, 0.9)$recommendation))/1000,
           sum(replicate(1000, STBP_case1(13, 0.9)$recommendation))/1000)


STCHAa <- c(sum(replicate(1000, STBP_case1(1, 0.5)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(2, 0.5)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(3, 0.5)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(4, 0.5)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(5, 0.5)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(6, 0.5)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(7, 0.5)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(8, 0.5)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(9, 0.5)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(10, 0.5)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(11, 0.5)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(12, 0.5)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(13, 0.5)$recommendation))/1000)


STCHAb <- c(sum(replicate(1000, STBP_case1(1, 0.9)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(2, 0.9)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(3, 0.9)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(4, 0.9)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(5, 0.9)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(6, 0.9)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(7, 0.9)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(8, 0.9)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(9, 0.1)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(10, 0.1)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(11, 0.1)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(12, 0.1)$recommendation))/1000,
            sum(replicate(1000, STBP_case1(13, 0.1)$recommendation))/1000)

correct1 <- c(rep(1, 8), rep(0, 5))


# Sample size


SPRTAs <- c(sum(replicate(1000, SPRT_case1(1)$samples))/1000,
            sum(replicate(1000, SPRT_case1(2)$samples))/1000,
            sum(replicate(1000, SPRT_case1(3)$samples))/1000,
            sum(replicate(1000, SPRT_case1(4)$samples))/1000,
            sum(replicate(1000, SPRT_case1(5)$samples))/1000,
            sum(replicate(1000, SPRT_case1(6)$samples))/1000,
            sum(replicate(1000, SPRT_case1(7)$samples))/1000,
            sum(replicate(1000, SPRT_case1(8)$samples))/1000,
            sum(replicate(1000, SPRT_case1(9)$samples))/1000,
            sum(replicate(1000, SPRT_case1(10)$samples))/1000,
            sum(replicate(1000, SPRT_case1(11)$samples))/1000,
            sum(replicate(1000, SPRT_case1(12)$samples))/1000,
            sum(replicate(1000, SPRT_case1(13)$samples))/1000)

STCHAs <- c(sum(replicate(1000, STBP_case1(1, 0.1)$samples))/1000,
            sum(replicate(1000, STBP_case1(2, 0.1)$samples))/1000,
            sum(replicate(1000, STBP_case1(3, 0.1)$samples))/1000,
            sum(replicate(1000, STBP_case1(4, 0.1)$samples))/1000,
            sum(replicate(1000, STBP_case1(5, 0.1)$samples))/1000,
            sum(replicate(1000, STBP_case1(6, 0.1)$samples))/1000,
            sum(replicate(1000, STBP_case1(7, 0.1)$samples))/1000,
            sum(replicate(1000, STBP_case1(8, 0.1)$samples))/1000,
            sum(replicate(1000, STBP_case1(9, 0.5)$samples))/1000,
            sum(replicate(1000, STBP_case1(10, 0.9)$samples))/1000,
            sum(replicate(1000, STBP_case1(11, 0.9)$samples))/1000,
            sum(replicate(1000, STBP_case1(12, 0.9)$samples))/1000,
            sum(replicate(1000, STBP_case1(13, 0.9)$samples))/1000)


STCHAas <- c(sum(replicate(1000, STBP_case1(1, 0.5)$samples))/1000,
             sum(replicate(1000, STBP_case1(2, 0.5)$samples))/1000,
             sum(replicate(1000, STBP_case1(3, 0.5)$samples))/1000,
             sum(replicate(1000, STBP_case1(4, 0.5)$samples))/1000,
             sum(replicate(1000, STBP_case1(5, 0.5)$samples))/1000,
             sum(replicate(1000, STBP_case1(6, 0.5)$samples))/1000,
             sum(replicate(1000, STBP_case1(7, 0.5)$samples))/1000,
             sum(replicate(1000, STBP_case1(8, 0.5)$samples))/1000,
             sum(replicate(1000, STBP_case1(9, 0.5)$samples))/1000,
             sum(replicate(1000, STBP_case1(10, 0.5)$samples))/1000,
             sum(replicate(1000, STBP_case1(11, 0.5)$samples))/1000,
             sum(replicate(1000, STBP_case1(12, 0.5)$samples))/1000,
             sum(replicate(1000, STBP_case1(13, 0.5)$samples))/1000)


STCHAbs <- c(sum(replicate(1000, STBP_case1(1, 0.9)$samples))/1000,
             sum(replicate(1000, STBP_case1(2, 0.9)$samples))/1000,
             sum(replicate(1000, STBP_case1(3, 0.9)$samples))/1000,
             sum(replicate(1000, STBP_case1(4, 0.9)$samples))/1000,
             sum(replicate(1000, STBP_case1(5, 0.9)$samples))/1000,
             sum(replicate(1000, STBP_case1(6, 0.9)$samples))/1000,
             sum(replicate(1000, STBP_case1(7, 0.9)$samples))/1000,
             sum(replicate(1000, STBP_case1(8, 0.9)$samples))/1000,
             sum(replicate(1000, STBP_case1(9, 0.1)$samples))/1000,
             sum(replicate(1000, STBP_case1(10, 0.1)$samples))/1000,
             sum(replicate(1000, STBP_case1(11, 0.1)$samples))/1000,
             sum(replicate(1000, STBP_case1(12, 0.1)$samples))/1000,
             sum(replicate(1000, STBP_case1(13, 0.1)$samples))/1000)
