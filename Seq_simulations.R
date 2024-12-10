# Sequential simulations

# This code runs all simulations with conventional sequential computing 
# processing. This takes longer than parallel computing. Code that runs 
# simulations with parallel computing is provided in the corresponding 
# cases' files.


#---------------------------------------------------------------------------
# Case 1
#---------------------------------------------------------------------------

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


#---------------------------------------------------------------------------
# Metrics case 1
#---------------------------------------------------------------------------

# Overall error rate
sum(1 - (1 - abs(correct1 - SPRTA))) / 13 # for SPRT
sum(1 - (1 - abs(correct1 - STCHA))) / 13 # STBP, correct init priors
sum(1 - (1 - abs(correct1 - STCHAa))) / 13 # STBP, naive init priors
sum(1 - (1 - abs(correct1 - STCHAb))) / 13 # STBP, incorrect init priors

# Type II error
sum(1 - (1 - abs(correct1[9:13] - SPRTA[9:13]))) / 5 # for SPRT
sum(1 - (1 - abs(correct1[9:13] - STCHA[9:13]))) / 5 # STBP, correct init priors
sum(1 - (1 - abs(correct1[9:13] - STCHAa[9:13]))) / 5 # STBP, naive init priors
sum(1 - (1 - abs(correct1[9:13] - STCHAb[9:13]))) / 5 # STBP, incorrect init priors

# Type I error
sum(1 - (1 - abs(correct1[1:8] - SPRTA[1:8]))) / 8 # for SPRT
sum(1 - (1 - abs(correct1[1:8] - STCHA[1:8]))) / 8 # STBP, correct init priors
sum(1 - (1 - abs(correct1[1:8] - STCHAa[1:8]))) / 8 # STBP, naive init priors
sum(1 - (1 - abs(correct1[1:8] - STCHAb[1:8]))) / 8 # STBP, incorrect init priors

# Mean sample sizes required
mean(SPRTAs)
mean(STCHAs)
mean(STCHAas)
mean(STCHAbs)




#---------------------------------------------------------------------------
# Case 2
#---------------------------------------------------------------------------

# Decisions

repl_SPRT <- function(ns){
  levels <- seq(0.1, 1, 0.1)
  result <- rep(NA, 10)
  
  
  for(i in 1: 10) {
    result[i] <- replicate(1000, simu_SPRT(levels[i], ns)$result) |>
      sum() / 1000
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
    result[i] <- replicate(1000,
                           STBP_case2(levels[i], ns, prior1 = priors[i])$result
    ) |>
      sum() / 1000
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
    result[i] <- replicate(1000, STBP_case2(levels[i], ns)$result) |>
      sum() / 1000
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
    result[i] <- replicate(1000,
                           STBP_case2(levels[i], ns, prior1 = priors[i])$result
    ) |>
      sum() / 1000
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
    result[i] <- replicate(1000, simu_SPRT(levels[i], ns)$bouts) |>
      sum() / 1000
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
    result[i] <- replicate(1000,
                           STBP_case2(levels[i], ns, prior1 = priors[i])$bouts
    ) |>
      sum() / 1000
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
    result[i] <- replicate(1000,
                           STBP_case2(levels[i], ns)$bouts) |>
      sum() / 1000
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
    result[i] <- replicate(1000,
                           STBP_case2(levels[i], ns, prior1 = priors[i])$bouts
    ) |>
      sum() / 1000
  }
  result
}

result5CPAs2 <- repl_SCPTAs2(5)
result10CPAs2 <- repl_SCPTAs2(10)
result20CPAs2 <- repl_SCPTAs2(20)
result30CPAs2 <- repl_SCPTAs2(30)

correct2 <- c(rep(1, 5), rep(0, 5))


#---------------------------------------------------------------------------
# Metrics case 2
#---------------------------------------------------------------------------

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






#---------------------------------------------------------------------------
# Case 3
#---------------------------------------------------------------------------

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