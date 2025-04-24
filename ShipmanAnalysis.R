# ANALYSIS OF DEATH CERTIFICATES SIGNED BY HAROLD SHIPMAN BETWEEN 1977 AND 1998

#---------------------------------------------------------------------------
# Load datasets (Baker & Donaldson, 2001)
#---------------------------------------------------------------------------

female_vict <- read.csv("ShipmanData/female_vict.csv")
male_vict <- read.csv("ShipmanData/male_vict.csv")

#---------------------------------------------------------------------------
# Run STBP
#---------------------------------------------------------------------------

# Obtain cumulative deaths to compare with cumulative expected mortality

female_vict$CSdeaths <- cumsum(female_vict$deaths)
male_vict$CSdeaths <- cumsum(male_vict$deaths)

source("STBP.R")

# Obtain sequential posterior probabilities for female and male deaths for the 
# H1 that CSdeaths is 33.3% greater than CSexpected_est

fem_probs <- stbp(data = female_vict$CSdeaths, 
                  hypotheses = 1.333 * female_vict$CSexpected_est,
                  likelihood_func = function(data, x) {
                    prod(dpois(data, lambda = x))
                  },
                  prior = 0.5,
                  lower_criterion = 0,
                  upper_criterion = 1)$probabilities

male_probs <- stbp(data = male_vict$CSdeaths, 
                   hypotheses = 1.333 * male_vict$CSexpected_est,
                   likelihood_func = function(data, x) {
                     prod(dpois(data, lambda = x))
                   },
                   prior = 0.5,
                   lower_criterion = 0,
                   upper_criterion = 1)$probabilities


#---------------------------------------------------------------------------
# Figure
#---------------------------------------------------------------------------

jpeg("Blog_versions/Fig3Post.jpeg", res = 96, width = 930, height = 670)

par(mfrow = c(2, 1))

par(mar = c(2, 6, 4, 2) + 0.1)
plot(female_vict$Year, (female_vict$CSdeaths - female_vict$CSexpected_est), 
     type  = "o", ylim = c(-20, 180), ylab = "Cumulative\nexcess mortality", 
     xlab = "", lwd = 2, cex.lab = 1.8, cex.axis = 1.8, col = "red", xaxt = "n",
     xlim = c(1977, 1999))
text(1975, 210, "A", xpd = TRUE, cex = 2)

axis(1, at = seq(1980, 2000, 5), labels = FALSE)
points(male_vict$Year, (male_vict$CSdeaths - male_vict$CSexpected_est), 
       type  = "o", pch =5, lwd = 2, col = "blue")
abline(h = 0, lty = 3)

par(mar = c(5, 6, 1, 2) + 0.1)
plot(c(1976.5, female_vict$Year), fem_probs, type = "o", ylim = c(0, 1), 
     ylab = expression(paste("Prob(H"["1"], "|x)")), 
     xlab = "Year", lwd = 2, cex.lab = 1.8, cex.axis = 1.8, 
     col  = "red", xlim = c(1977, 1999))
abline(h = 0.99, lty = 2)
abline(h = 0.95, lty = 2)
abline(h = 0, lty = 3)
text(1975, 1.1, "B", xpd = TRUE, cex = 2)

points(c(1976.5, male_vict$Year), male_probs, pch =5, lwd = 2, type = "o", 
       col  = "blue")

dev.off()
