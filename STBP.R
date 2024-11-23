stbp_posterior <- function(data,
                           hypothesis,
                           likelihood_func,
                           prior,
                           lower_bnd = 0,
                           upper_bnd = Inf) {
  likelihood <- function(x) {
    prod(likelihood_func(data, x))
  }

  H1 <- prior *
            integrate(
              Vectorize(likelihood),
              lower = hypothesis,
              upper = upper_bnd
            )$value
  H0 <- (1 - prior) *
            integrate(
               Vectorize(likelihood),
               lower = lower_bnd,
               upper = hypothesis
            )$value
  posterior <- H1 / (H0 + H1)
  posterior
}

stbp <- function(data,
                 hypotheses,
                 likelihood_func,
                 prior = 0.5,
                 lower_bnd = 0,
                 upper_bnd = Inf,
                 lower_criterion = 0.01,
                 upper_criterion = 0.99) {
  # If hypothesis is just a single repeated value,
  # make a vector of that value repeated as many times as there are bouts.
  # This makes it so that the user can input either a single hypothesis,
  # or a trajectory of hypotheses.
  if(length(hypotheses) == 1) hypotheses <- rep(hypotheses, length(data))

  # Init vector with length equal to number of sampling bouts
  # and with initial prior as its first value
  posteriors <- c(prior, rep(NA, length(data) - 1))
  for(i in seq_along(data)) {
    bout = data[i]
    posteriors[i + 1] = stbp_posterior(bout,
                                     hypotheses[i],
                                     likelihood_func,
                                     prior = posteriors[i],
                                     lower_bnd = lower_bnd,
                                     upper_bnd = upper_bnd)
    # Break from iteration early if early_return is true,
    # the minimum iterations have been reached, and
    # if either of the decision criteria have been reached
    if(((posteriors[i + 1] < lower_criterion) ||
        (posteriors[i + 1] > upper_criterion))
      ) break
  }

  response <- if(posteriors[i + 1] < lower_criterion) 1 else 0
  indices_above <- which(posteriors > upper_criterion)
  indices_below <- which(posteriors < lower_criterion)

  return(list(
    probabilities = posteriors,
    recommendation = response,
    num_iterations = i,
    decision_indices = list(
      above = indices_above,
      below = indices_below
    )
  ))
}
