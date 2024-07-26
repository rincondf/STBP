# FIGURE 1

require(khroma)

batlowW <- color("batlowW")
imola <- color("imola")

par(mfrow = c(2, 1))

par(mar = c(2, 9, 5, 2) + 0.1)

plot(seq(1, 13), 1 - abs(correct1 - SPRTA), ylim = c(0, 1), type = "o", 
     ylab = "Proportion of\ncorrect decisions\n", xlab = "", lwd = 2, cex.lab = 1.8, 
     xaxt = "n", yaxt = "n")

axis(1, at = seq(2, 12, 2), labels = FALSE)
axis(2, at = seq(0, 1, 0.2), labels = TRUE, cex.axis = 1.8, las = 2)

points(seq(1, 13), 1 - abs(correct1 - STCHA), col = imola(10)[1], type = "o", lwd = 2)
points(seq(1, 13), 1 - abs(correct1 - STCHAa), col = imola(10)[6], type = "o", lwd = 2)
points(seq(1, 13), 1 - abs(correct1 - STCHAb), col = imola(10)[8], type = "o", lwd = 2)

mtext("a", side = 3, cex = 2, outer = FALSE, line = 1, adj = 0, at = -4, xpd  = NA)

abline(v = 9, lty = 3, lwd = 2)

par(mar = c(5, 9, 2, 2) + 0.1)

plot(seq(1, 13), SPRTAs, type = "o", ylim = c(0, 50), ylab = "Average No. of\nsamples\n", 
     xlab = "Population size", lwd = 2, cex.lab = 1.8, cex.axis = 1.8,
     yaxt = "n")
points(seq(1, 13), STCHAs, col = imola(10)[1], type = "o", lwd = 2)
points(seq(1, 13), STCHAas, col = imola(10)[6], type = "o", lwd = 2)
points(seq(1, 13), STCHAbs, col = imola(10)[8], type = "o", lwd = 2)

axis(2, at = seq(0, 50, 10), labels = TRUE, cex.axis = 1.8, las = 2)

abline(v = 9, lty = 3, lwd = 2)

mtext("b", side = 3, cex = 2, outer = FALSE, line = 1, adj = 0, at = -4, xpd  = NA)




#  FIGURE 2

require(khroma)

batlowW <- color("batlowW")
imola <- color("imola")

par(mfrow = c(2, 3), oma = c(3, 5, 0, 0))


par(mar = c(0.5, 4, 4, 0) + 0.1)
plot(levels, 1 - abs(correct2 - result5), type = "o", xlab = "", ylab = "", lwd = 2, 
     cex.axis = 1.8, ylim = c(0, 1), yaxt = "n", xaxt = "n")

title(ylab = "Proportion of\ncorrect decisions\n", cex.lab = 2, xpd = NA)

axis(2, at = seq(0, 1, 0.2), labels = TRUE, cex.axis = 1.8, las = 2)
axis(1, at = seq(0, 1, 0.2), labels = FALSE)

points(levels, 1 - abs(correct2 - result10), type = "o", lwd = 2, lty = 2)
points(levels, 1 - abs(correct2 - result20), type = "o", lwd = 2, lty = 3)
points(levels, 1 - abs(correct2 - result30), type = "o", lwd = 2, lty = 4)

points(levels, 1 - abs(correct2 - result5CPA1), type = "o", col = imola(10)[1], lwd = 2)
points(levels, 1 - abs(correct2 - result10CPA1), type = "o", col = imola(10)[1], lwd = 2, lty = 2)
points(levels, 1 - abs(correct2 - result20CPA1), type = "o", col = imola(10)[1], lwd = 2, lty = 3)
points(levels, 1 - abs(correct2 - result30CPA1), type = "o", col = imola(10)[1], lwd = 2, lty = 4)

abline(v = 0.5, lty = 3, lwd = 2)

mtext("a", side = 3, cex = 1.8, line = 0.5, adj = 0, at = -0.02, xpd  = NA)


par(mar = c(0.5, 2, 4, 2) + 0.1)
plot(levels, 1 - abs(correct2 - result5), type = "o", ylab = "", xlab = "", lwd = 2, cex.lab = 1.8, 
     cex.axis = 1.8, ylim = c(0, 1), yaxt = "n", xaxt = "n")
axis(1, at = seq(0, 1, 0.2), labels = FALSE)
axis(2, at = seq(0, 1, 0.2), labels = FALSE)

points(levels, 1 - abs(correct2 - result10), type = "o", lwd = 2, lty = 2)
points(levels, 1 - abs(correct2 - result20), type = "o", lwd = 2, lty = 3)
points(levels, 1 - abs(correct2 - result30), type = "o", lwd = 2, lty = 4)


points(levels, 1 - abs(correct2 - result5CPA), type = "o", col = imola(10)[6], lwd = 2)
points(levels, 1 - abs(correct2 - result10CPA), type = "o", col = imola(10)[6], lwd = 2, lty = 2)
points(levels, 1 - abs(correct2 - result20CPA), type = "o", col = imola(10)[6], lwd = 2, lty = 3)
points(levels, 1 - abs(correct2 - result30CPA), type = "o", col = imola(10)[6], lwd = 2, lty = 4)

abline(v = 0.5, lty = 3, lwd = 2)

mtext("b", side = 3, cex = 1.8, line = 0.5, adj = 0, at = -0.02, xpd  = NA)

par(mar = c(0.5, 0, 4, 4) + 0.1)
plot(levels, 1 - abs(correct2 - result5), type = "o", ylab = "", xlab = "", lwd = 2, cex.lab = 1.8, 
     cex.axis = 1.8, ylim = c(0, 1), yaxt = "n", xaxt = "n")
axis(1, at = seq(0, 1, 0.2), labels = FALSE)
axis(2, at = seq(0, 1, 0.2), labels = FALSE)

points(levels, 1 - abs(correct2 - result10), type = "o", lwd = 2, lty = 2)
points(levels, 1 - abs(correct2 - result20), type = "o", lwd = 2, lty = 3)
points(levels, 1 - abs(correct2 - result30), type = "o", lwd = 2, lty = 4)

points(levels, 1 - abs(correct2 - result5CPA2), type = "o", col = imola(10)[8], lwd = 2)
points(levels, 1 - abs(correct2 - result10CPA2), type = "o", col = imola(10)[8], lwd = 2, lty = 2)
points(levels, 1 - abs(correct2 - result20CPA2), type = "o", col = imola(10)[8], lwd = 2, lty = 3)
points(levels, 1 - abs(correct2 - result30CPA2), type = "o", col = imola(10)[8], lwd = 2, lty = 4)

abline(v = 0.5, lty = 3, lwd = 2)

mtext("c", side = 3, cex = 1.8, line = 0.5, adj = 0, at = -0.02, xpd  = NA)


par(mar = c(2, 4, 2.5, 0) + 0.1)
plot(levels, result5s, type = "o", ylim = c(0, 9), ylab = "", 
     xlab = "", lwd = 2, cex.lab = 1.8, cex.axis = 1.8, yaxt = "n")

title(ylab = "Average No. of\nsampling bouts\n", cex.lab = 2, xpd = NA)

axis(2, at = seq(0, 8, 2), labels = TRUE, cex.axis = 1.8, las = 2)
axis(1, at = seq(0, 1, 0.2), labels = FALSE)

points(levels, result10s, type = "o", lty = 2, lwd = 2)
points(levels, result20s, type = "o", lty = 3, lwd = 2)
points(levels, result30s, type = "o", lty = 4, lwd = 2)

points(levels, result5CPAs1, type = "o", col = imola(10)[1], lwd = 2)
points(levels, result10CPAs1, type = "o", col = imola(10)[1], lwd = 2, lty = 2)
points(levels, result20CPAs1, type = "o", col = imola(10)[1], lwd = 2, lty = 3)
points(levels, result30CPAs1, type = "o", col = imola(10)[1], lwd = 2, lty = 4)

abline(v = 0.5, lty = 3, lwd = 2)

mtext("d", side = 3, cex = 1.8, line = 0.5, adj = 0, at = -0.02, xpd  = NA)

par(mar = c(2, 2, 2.5, 2) + 0.1)
plot(levels, result5s, type = "o", ylim = c(0, 9), ylab = "", 
     xlab = "", lwd = 2, cex.lab = 1.8, cex.axis = 1.8, yaxt = "n")
axis(1, at = seq(0, 1, 0.2), labels = FALSE)
axis(2, at = seq(0, 8, 2), labels = FALSE)

points(levels, result10s, type = "o", lty = 2, lwd = 2)
points(levels, result20s, type = "o", lty = 3, lwd = 2)
points(levels, result30s, type = "o", lty = 4, lwd = 2)

points(levels, result5CPAs, type = "o", col = imola(10)[6], lwd = 2)
points(levels, result10CPAs, type = "o", col = imola(10)[6], lwd = 2, lty = 2)
points(levels, result20CPAs, type = "o", col = imola(10)[6], lwd = 2, lty = 3)
points(levels, result30CPAs, type = "o", col = imola(10)[6], lwd = 2, lty = 4)

abline(v = 0.5, lty = 3, lwd = 2)

mtext("e", side = 3, cex = 1.8, line = 0.5, adj = 0, at = -0.02, xpd  = NA)


par(mar = c(2, 0, 2.5, 4) + 0.1)
plot(levels, result5s, type = "o", ylim = c(0, 9), ylab = "", 
     xlab = "", lwd = 2, cex.lab = 1.8, cex.axis = 1.8, yaxt = "n")
axis(1, at = seq(0, 1, 0.2), labels = FALSE)
axis(2, at = seq(0, 8, 2), labels = FALSE)

points(levels, result10s, type = "o", lty = 2, lwd = 2)
points(levels, result20s, type = "o", lty = 3, lwd = 2)
points(levels, result30s, type = "o", lty = 4, lwd = 2)

points(levels, result5CPAs2, type = "o", col = imola(10)[8], lwd = 2)
points(levels, result10CPAs2, type = "o", col = imola(10)[8], lwd = 2, lty = 2)
points(levels, result20CPAs2, type = "o", col = imola(10)[8], lwd = 2, lty = 3)
points(levels, result30CPAs2, type = "o", col = imola(10)[8], lwd = 2, lty = 4)

abline(v = 0.5, lty = 3, lwd = 2)

mtext("f", side = 3, cex = 1.8, line = 0.5, adj = 0, at = -0.02, xpd  = NA)

title(xlab = "Outbreak severity", cex.lab = 2.5, outer = TRUE, line = 1.8)



# FIGURE 3

require(khroma)

batlowW <- color("batlowW")
imola <- color("imola")

par(mfrow = c(2, 1))

par(mar = c(2, 9, 5, 4) + 0.1)

plot(means_det, 1-beta1, type = "o", ylab = "Beta\n", xlab = "", lwd = 2, 
     cex.lab = 1.8, ylim = c(0, 1), col = imola(10)[6], yaxt = "n", xaxt = "n")

axis(1, at = seq(0, 0.2, 0.05), labels = FALSE)
axis(2, at = seq(0, 1, 0.2), labels = TRUE, cex.axis = 1.8, las = 2)

axis(4, at = 0.05, labels = TRUE, cex.axis = 1.8, las = 2)

points(means_det, 1-beta3, type = "o", lwd = 2, lty = 2, col = imola(10)[6])
points(means_det, 1-beta5, type = "o", lwd = 2, lty = 3, col = imola(10)[6])
points(means_det, 1-beta10, type = "o", lwd = 2, lty = 4, col = imola(10)[6])


points(means_det, 1-betaGreen30, type = "o", lwd = 2, lty = 4)
points(means_det, 1-betaGreen20, type = "o", lwd = 2, lty = 3)
points(means_det, 1-betaGreen10, type = "o", lwd = 2, lty = 2)

lines(seq(0.01, 0.2, 0.001), beta_fun(seq(0.01, 0.2, 0.001), 30), lwd = 2)

abline(h = 0.05, col = "grey36", lwd = 2)

mtext("a", side = 3, cex = 2, line = 1, adj = 0, at = -0.05, xpd  = NA)

par(mar = c(5, 9, 2, 4) + 0.1)

plot(means_det, size1, type = "o", ylab = "Average No. of\nsampling bouts\n", xlab = "Mean individuals per sampling unit", 
     lwd = 2, cex.lab = 1.8, cex.axis = 1.8, ylim = c(0, 20), col = imola(10)[6], yaxt = "n")

axis(2, at = seq(0, 20, 5), labels = TRUE, cex.axis = 1.8, las = 2)

points(means_det, size3, type = "o", lwd = 2, lty = 2, col = imola(10)[6])
points(means_det, size5, type = "o", lwd = 2, lty = 3, col = imola(10)[6])
points(means_det, size10, type = "o", lwd = 2, lty = 4, col = imola(10)[6])

mtext("b", side = 3, cex = 2, line = 1, adj = 0, at = -0.05, xpd  = NA)



# FIGURE 4

require(khroma)

batlowW <- color("batlowW")
imola <- color("imola")

fig4 <- matrix(c(s1alt, size1_0, s1a, s3alt, size3_0, s3a, s5alt, size5_0, s5a, s10alt, size10_0, s10a), 3, 4)

par(mar = c(5, 9, 2, 2) + 0.1)
barplot(fig4, names.arg = c(1, 3, 5, 10),  yaxt = "n", 
        cex.names = 2, col = c(imola(10)[8], imola(10)[6], imola(10)[1]), border = NA, beside = TRUE)

axis(2, at = seq(0, 20, 5), labels = TRUE, cex.axis = 1.8, las = 2)
title(ylab = "No. of sampling bouts", cex.lab = 2, xpd = NA, line = 4)
title(xlab = "Subsample size", cex.lab = 2, xpd = NA, line = 3.5)
