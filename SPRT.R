#SPRT

require(truncdist)

k_est1 <- function(me) {
  a = exp(0.6043225)
  b = 1.218041
  (me^2) / ((a * me^(b)) - me)
}

ka_es_r <- function(m) {
  a <- exp(0.6043225)
  b <- 1.218041
  times <- length(m)
  a1 <- rep(NA, times)
  for(i in 1:times) {
    a1[i] <- (m[i]^2) /
      ((a * m[i]^(b) *
          exp(rtrunc(1, "norm", a = log(1 / (a * m[i]^(b - 1))),
                     b = Inf, mean = 0, sd = 0.3222354)))
       - m[i])
  }
  a1
}

k_um <- k_est1(9)

low_int_nb <- function(al, be, me0, me1, k_est){
  (log(be / (1 - al))) / (log((me1 * (me0 + k_est)) / (me0 * (me1 + k_est))))
}

low_um_p <- low_int_nb(al = 0.1, be = 0.1, me0 = 8, me1 = 10, k_est = k_um)


hi_int_nb <- function(al, be, me0, me1, k_est){
  (log((1 - be) / (al))) / (log((me1 * (me0 + k_est)) / (me0 * (me1 + k_est))))
}


hi_um_p <- hi_int_nb(al = 0.1, be = 0.1, me0 = 8, me1 = 10, k_est = k_um)

ll_sl_nb <- function(al, be, me0, me1, k_est){
  (k_est * log((me1 + k_est) / (me0 + k_est))) /
    (log((me1 * (me0 + k_est)) / (me0 * (me1 + k_est))))
}

pen_p <- ll_sl_nb(al = 0.1, be = 0.1, me0 = 8, me1 = 10, k_est = k_um)

low_seq_c <- function(x){
  pen_p * x + low_um_p
}

hi_seq_c <- function(x){
  pen_p * x + hi_um_p
}


proced_SPRT_cont <- function(d){
  acu <- rep(NA, 100)
  fi <- rnbinom(mu = d, size = ka_es_r(d), n = 6000)
  mean_Re <- mean(fi)
  for(i in 1:100){
    acu[i] <- sample(fi, size = 1, replace = F)
    fi <- fi[-match(acu[i], fi)]
    dat <- cumsum(acu)[i]
    if ((dat < low_seq_c(i)) | (dat > hi_seq_c(i))) break
  }
  if (dat < low_seq_c(i)) {
    resp <- 1
  } else {
    resp <- 0
  }
  return(list(estimate = dat / i, samples = i, recommendation = resp, col_data = acu,
              mean = mean_Re))
}

SPRTA <- c(sum(replicate(1000, proced_SPRT_cont(1)$recommendation))/1000,
           sum(replicate(1000, proced_SPRT_cont(2)$recommendation))/1000,
           sum(replicate(1000, proced_SPRT_cont(3)$recommendation))/1000,
           sum(replicate(1000, proced_SPRT_cont(4)$recommendation))/1000,
           sum(replicate(1000, proced_SPRT_cont(5)$recommendation))/1000,
           sum(replicate(1000, proced_SPRT_cont(6)$recommendation))/1000,
           sum(replicate(1000, proced_SPRT_cont(7)$recommendation))/1000,
           sum(replicate(1000, proced_SPRT_cont(8)$recommendation))/1000,
           sum(replicate(1000, proced_SPRT_cont(9)$recommendation))/1000,
           sum(replicate(1000, proced_SPRT_cont(10)$recommendation))/1000,
           sum(replicate(1000, proced_SPRT_cont(11)$recommendation))/1000,
           sum(replicate(1000, proced_SPRT_cont(12)$recommendation))/1000,
           sum(replicate(1000, proced_SPRT_cont(13)$recommendation))/1000)


####

CP.v3A <- function(dat, prd, prior) {
  lld1 <- function(x) {
    prod(dnbinom(dat, mu = x, size = k_um))
  }
  
  Pr3A <- (integrate(Vectorize(lld1), lower = prd, upper = Inf)$value * prior) /
    (integrate(Vectorize(lld1), lower = 0, upper = prd)$value * (1 - prior) + (integrate(Vectorize(lld1), lower = prd, upper = Inf)$value * prior))
  Pr3A
}


proced_STCH_cont <- function(d, prior){
  acu <- rep(NA, 100)
  fi <- rnbinom(mu = d, size = ka_es_r(d), n = 6000)
  mean_Re <- mean(fi)
  probs <- c(prior, rep(NA, 99))
  for(i in 1:100){
    acu[i] <- sample(fi, size = 1, replace = F)
    fi <- fi[-match(acu[i], fi)]
    probs[i + 1] <- CP.v3A(dat = acu[i], prd = 9, prior = probs[i])
    if ((length(!is.na(acu)) > 3) & ((probs[i + 1] < 0.01) | (probs[i + 1] > 0.99))) break
  }
  if (probs[i + 1] < 0.01) {
    resp <- 1
  } else {
    resp <- 0
  }
  return(list(Probabilities = probs, samples = i, recommendation = resp, col_data = acu,
              mean = mean_Re))
}


STCHA <- c(sum(replicate(1000, proced_STCH_cont(1, 0.1)$recommendation))/1000,
           sum(replicate(1000, proced_STCH_cont(2, 0.1)$recommendation))/1000,
           sum(replicate(1000, proced_STCH_cont(3, 0.1)$recommendation))/1000,
           sum(replicate(1000, proced_STCH_cont(4, 0.1)$recommendation))/1000,
           sum(replicate(1000, proced_STCH_cont(5, 0.1)$recommendation))/1000,
           sum(replicate(1000, proced_STCH_cont(6, 0.1)$recommendation))/1000,
           sum(replicate(1000, proced_STCH_cont(7, 0.1)$recommendation))/1000,
           sum(replicate(1000, proced_STCH_cont(8, 0.1)$recommendation))/1000,
           sum(replicate(1000, proced_STCH_cont(9, 0.5)$recommendation))/1000,
           sum(replicate(1000, proced_STCH_cont(10, 0.9)$recommendation))/1000,
           sum(replicate(1000, proced_STCH_cont(11, 0.9)$recommendation))/1000,
           sum(replicate(1000, proced_STCH_cont(12, 0.9)$recommendation))/1000,
           sum(replicate(1000, proced_STCH_cont(13, 0.9)$recommendation))/1000)


STCHAa <- c(sum(replicate(1000, proced_STCH_cont(1, 0.5)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(2, 0.5)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(3, 0.5)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(4, 0.5)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(5, 0.5)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(6, 0.5)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(7, 0.5)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(8, 0.5)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(9, 0.5)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(10, 0.5)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(11, 0.5)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(12, 0.5)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(13, 0.5)$recommendation))/1000)


STCHAb <- c(sum(replicate(1000, proced_STCH_cont(1, 0.9)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(2, 0.9)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(3, 0.9)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(4, 0.9)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(5, 0.9)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(6, 0.9)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(7, 0.9)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(8, 0.9)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(9, 0.1)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(10, 0.1)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(11, 0.1)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(12, 0.1)$recommendation))/1000,
            sum(replicate(1000, proced_STCH_cont(13, 0.1)$recommendation))/1000)

correct1 <- c(rep(1, 8), rep(0, 5))
