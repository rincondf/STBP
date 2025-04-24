# COMPARISON BETWEEN EXPLICIT SPECIFICATION OF PRIOR DISTRIBUTIONS 
# AND THE USE OF IMPROPER PRIORS FOR DIFFERENT SAMPLE SIZES

#---------------------------------------------------------------------------
# DATA MATRICES WITH DIFFERENT SAMPLE SIZES FOR x = 0:10
#---------------------------------------------------------------------------

# n = 1

d1 <- matrix(rep(seq(0, 10), each = 1), 1, 11)

# n = 2

d2 <- matrix(rep(seq(0, 10), each = 2), 2, 11)

# n = 5

d5 <- matrix(rep(seq(0, 10), each = 5), 5, 11)

# n = 10

d10 <- matrix(rep(seq(0, 10), each = 10), 10, 11)

#---------------------------------------------------------------------------
# FUNCTIONS OBTAIN POSTERIOR PROBABILITIES FOR H1: MU > 5
#---------------------------------------------------------------------------

# Configuration 1 of Gamma conjugate prior

conjG1 <- function(data){
  probs <- rep(NA, 11)
  
  for(i in 1: length(data[1, ])) {
    probs[i] <- (1 - pgamma(5, shape = 1 + sum(data[, i]), 
                            rate = 0.139 + length(data[, i])))
  }
  
  probs
}

# Configuration 2 of Gamma conjugate prior

conjG2 <- function(data){
  probs <- rep(NA, 11)
  
  for(i in 1: length(data[1, ])) {
    probs[i] <- (1 - (pgamma(5, shape = 10 + sum(data[, i]), 
                             rate = 1.934 + length(data[, i]))))
  }
  
  probs
}

# Sequential test of Bayesian posterior probabilities

STBP_P <- function(data) {
  
  posterior <- rep(NA, length(data[1, ]))
  
  for(i in 1: length(data[1, ])) {
    likelihood_func <- function(data, lambda) {
      dpois(data, lambda)
    }
    
    likelihood <- function(x) {
      prod(likelihood_func(data = data[, i], lambda = x))
    }
    
    # Notice that if this likehood function is instead defined as:
    
    #likelihood <- function(x) {
    #  prod(likelihood_func(data = data[, i], lambda = x)) * 
    #    dgamma(x, shape, rate)
    #}
    
    # then STBP_P is equivalent to conjG1 or conjG2 with appropiate parameters 
    # for the gamma function
    
    H1 <- 0.5 *
      integrate(
        Vectorize(likelihood),
        lower = 5,
        upper = Inf
      )$value
    H0 <- (1 - 0.5) *
      integrate(
        Vectorize(likelihood),
        lower = 0,
        upper = 5
      )$value
    posterior[i] <- (H1) / ((H0) + (H1))
  }
  
  
  
  posterior
}



#---------------------------------------------------------------------------
# FIGURE
#---------------------------------------------------------------------------

jpeg("Blog versions/Fig2Post.jpeg", res = 96, width = 930, height = 590)

par(mfrow = c(2, 2))

par(mar = c(2, 5, 4, 2))
plot(seq(0, 10, 1), conjG1(d1), 
     ylim = c(0, 1), type = "l", col = "red", lwd = 2, cex.axis = 1.5,
     xlab = "", ylab = expression(paste("P(", mu  >= psi, "|X=x)")), 
     cex.lab = 2, xaxt = "n")

axis(1, at = seq(0, 10, 2), labels = FALSE)

lines(seq(0, 10, 1), conjG2(d1),
      lwd = 2)

lines(seq(0, 10), STBP_P(d1), ylim = c(0, 1), col = "blue", lwd = 2)

text(1, 0.8, "n = 1", cex = 2)

abline(v = seq(0, 10, 2), lty = 3, col = "grey")
abline(h = seq(0, 1, 0.2), lty = 3, col = "grey")




par(mar = c(2, 3, 4, 4))
plot(seq(0, 10), conjG1(d2), ylim = c(0, 1), type = "l", col = "red", lwd = 2, cex.axis = 1.5,
     xlab = "", ylab = "", cex.lab = 2, yaxt = "n", xaxt = "n")
lines(seq(0, 10), conjG2(d2), lwd = 2)
lines(seq(0, 10), STBP_P(d2), col = "blue", lwd = 2)

axis(1, at = seq(0, 10, 2), labels = FALSE)
axis(2, at = seq(0, 1, 0.2), labels = FALSE)

abline(v = seq(0, 10, 2), lty = 3, col = "grey")
abline(h = seq(0, 1, 0.2), lty = 3, col = "grey")

text(1, 0.8, "n = 2", cex = 2)



par(mar = c(4, 5, 2, 2))
plot(seq(0, 10), conjG1(d5), ylim = c(0, 1), type = "l", col = "red", lwd = 2, cex.axis = 1.5,
     xlab = "x", cex.lab = 2, ylab = expression(paste("P(", mu  >= psi, "|X=x)")), 
     cex.lab = 2)
lines(seq(0, 10), conjG2(d5), lwd = 2)
lines(seq(0, 10), STBP_P(d5), col = "blue", lwd = 2)

abline(v = seq(0, 10, 2), lty = 3, col = "grey")
abline(h = seq(0, 1, 0.2), lty = 3, col = "grey")

text(1, 0.8, "n = 5", cex = 2)



par(mar = c(4, 3, 2, 4))
plot(seq(0, 10), conjG1(d10), ylim = c(0, 1), type = "l", col = "red", lwd = 2, cex.axis = 1.5,
     xlab = "x", ylab = "", cex.lab = 2, yaxt = "n")
lines(seq(0, 10), conjG2(d10), lwd = 2)
lines(seq(0, 10), STBP_P(d10), col = "blue", lwd = 2)

axis(2, at = seq(0, 1, 0.2), labels = FALSE)

abline(v = seq(0, 10, 2), lty = 3, col = "grey")
abline(h = seq(0, 1, 0.2), lty = 3, col = "grey")

text(1, 0.8, "n = 10", cex = 2)

legend(6, 0.8, legend=c("Gamma(1, 0.139)", "Gamma(10, 1.934)", "No proper prior\n(our approach)"),  
       fill = c("red","black", "blue"), bty = "n", title = "Prior", title.adj = 0.2, title.cex = 1.2,
       y.intersp = 1.2)

dev.off()