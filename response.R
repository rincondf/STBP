par(mfrow = c(2, 2))

par(mar = c(2, 5, 4, 2))
plot(seq(0, 10, 0.1), (1-pgamma(5, shape = 1 + seq(0, 10, 0.1), rate = 0.139 +1)), 
     ylim = c(0, 1), type = "l", col = "red", lwd = 2, cex.axis = 1.5,
     xlab = "", ylab = expression(paste("P(", mu  >= psi, "|X=x)")), cex.lab = 2,
     xaxt = "n")

axis(1, at = seq(0, 10, 2), labels = FALSE)

lines(seq(0, 10, 0.1), (1-pgamma(5, shape = 10 + seq(0, 10, 0.1), rate = 1.934+1)),
      lwd = 2)

lines(seq(0, 10), test3(dd), ylim = c(0, 1), col = "blue", lwd = 2)

text(1, 0.8, "n = 1", cex = 2)

abline(v = seq(0, 10, 2), lty = 3, col = "grey")
abline(h = seq(0, 1, 0.2), lty = 3, col = "grey")




par(mar = c(2, 3, 4, 4))
plot(seq(0, 10), scn1, ylim = c(0, 1), type = "l", col = "red", lwd = 2, cex.axis = 1.5,
     xlab = "", ylab = "", cex.lab = 2, yaxt = "n", xaxt = "n")
lines(seq(0, 10), scn2, lwd = 2)
lines(seq(0, 10), test3(dd), col = "blue", lwd = 2)

axis(1, at = seq(0, 10, 2), labels = FALSE)
axis(2, at = seq(0, 1, 0.2), labels = FALSE)

abline(v = seq(0, 10, 2), lty = 3, col = "grey")
abline(h = seq(0, 1, 0.2), lty = 3, col = "grey")

text(1, 0.8, "n = 2", cex = 2)



par(mar = c(4, 5, 2, 2))
plot(seq(0, 10), scn1, ylim = c(0, 1), type = "l", col = "red", lwd = 2, cex.axis = 1.5,
     xlab = "x", cex.lab = 2, ylab = expression(paste("P(", mu  >= psi, "|X=x)")), 
     cex.lab = 2)
lines(seq(0, 10), scn2, lwd = 2)
lines(seq(0, 10), test3(dd), col = "blue", lwd = 2)

abline(v = seq(0, 10, 2), lty = 3, col = "grey")
abline(h = seq(0, 1, 0.2), lty = 3, col = "grey")

text(1, 0.8, "n = 5", cex = 2)



par(mar = c(4, 3, 2, 4))
plot(seq(0, 10), scn1, ylim = c(0, 1), type = "l", col = "red", lwd = 2, cex.axis = 1.5,
     xlab = "x", ylab = "", cex.lab = 2, yaxt = "n")
lines(seq(0, 10), scn2, lwd = 2)
lines(seq(0, 10), test3(dd), col = "blue", lwd = 2)

axis(2, at = seq(0, 1, 0.2), labels = FALSE)

abline(v = seq(0, 10, 2), lty = 3, col = "grey")
abline(h = seq(0, 1, 0.2), lty = 3, col = "grey")

text(1, 0.8, "n = 10", cex = 2)

legend(6, 0.8, legend=c("Gamma(1, 0.139)", "Gamma(10, 1.934)", "No proper prior\n(our approach)"),  
       fill = c("red","black", "blue"), bty = "n", title = "Prior", title.adj = 0.2, title.cex = 1.2,
       y.intersp = 1.2)









plot(seq(0, 10, 0.1), (1-pgamma(5, shape = 1 + seq(0, 10, 0.1), rate = 0.139 +1)), ylim = c(0, 1), type = "l", col = "red")


plot(seq(0, 10, 0.1), (1-pgamma(seq(0, 5), shape = 1 + seq(0, 10, 0.1), rate = 0.139 +1))*0.1 + 
       (1-pgamma(seq(5, 10), shape = 1 + seq(0, 10, 0.1), rate = 0.139 +1))*0.9)

points(seq(0, 10, 0.1), (1-pgamma(5, shape = 10 + seq(0, 10, 0.1), rate = 1.934 +1))*0.1 + 
       (1-pgamma(5, shape = 10 + seq(0, 10, 0.1), rate = 1.934 +1))*0.9)




test <- function(data) {
  
  for(i in 1: length(data[1, ])) {
    likelihood_func <- function(data, lambda) {
      dpois(data, lambda)
    }
    
    likelihood <- function(x) {
      prod(likelihood_func(data = data[, i], lambda = x)) * dgamma(x, shape = 1, rate = 0.139)
    }
    
    
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



test2 <- function(data) {
  
  for(i in 1: length(data[1, ])) {
    likelihood_func <- function(data, lambda) {
      dpois(data, lambda)
    }
    
    likelihood <- function(x) {
      prod(likelihood_func(data = data[, i], lambda = x)) * dgamma(x, shape = 10, rate = 1.934)
    }
    
    
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


test3 <- function(data) {
  
  for(i in 1: length(data[1, ])) {
    likelihood_func <- function(data, lambda) {
      dpois(data, lambda)
    }
    
    likelihood <- function(x) {
      prod(likelihood_func(data = data[, i], lambda = x))
    }
    
    
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

plot(seq(0, 10), test(dd), ylim = c(0, 1), type = "l", col = "red")
lines(seq(0, 10), test2(dd), ylim = c(0, 1))
lines(seq(0, 10), test3(dd), ylim = c(0, 1), col = "blue")
abline(h = 0.5)
abline(v = 5, lty = 3)


dd <- matrix(rep(seq(0, 10), each = 10), 10, 11)



scn1 <- rep(NA, 11)

for(i in 1: length(dd[1, ])) {
  scn1[i] <- (1-pgamma(5, shape = 1 + sum(dd[, i]), rate = 0.139 + length(dd[, i])))
}

scn2 <- rep(NA, 11)

for(i in 1: length(dd[1, ])) {
  scn2[i] <- (1 - (pgamma(5, shape = 10 + sum(dd[, i]), rate = 1.934 + length(dd[, i]))))
}



