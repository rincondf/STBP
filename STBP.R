stbp_posterior <- function(data,
                           hypothesis,
                           likelihood_func,
                           prior=0.5,
                           lower_bnd=-Inf,
                           upper_bnd=Inf) {
  likelihood <- function(x) {
    prod(likelihood_func(data, x))
  }

  null <- prior *
            integrate(
              Vectorize(likelihood),
              lower = hypothesis,
              upper = upper_bnd
            )$value
  alt <- (1 - prior) *
            integrate(
               Vectorize(likelihood),
               lower = lower_bnd,
               upper = hypothesis
            )$value
  posterior <- null / (alt + null)
  posterior
}

stbp <- function(data,
                 hypotheses,
                 likelihood_func,
                 prior=0.5,
                 lower_bnd=-Inf,
                 upper_bnd=Inf,
                 lower_criterion=0.01,
                 upper_criterion=0.99,
                 early_return=FALSE,
                 min_iterations=4,
                 posterior_bnd_lower = 0,
                 posterior_bnd_upper = 1,
                 custom_posterior_func = NA) {
  # If hypothesis is just a single repeated value,
  # make a vector of that value repeated as many times as there are bouts.
  # This makes it so that the user can input either a single hypothesis,
  # or a trajectory of hypotheses.
  if(length(hypotheses) == 1) hypotheses <- rep(hypotheses, length(data))

  # If a custom function to calculate the posterior is supplied,
  # shadow the usual `stbp_posterior` function with that
  if(!is.na(custom_posterior_func)) stbp_posterior <- custom_posterior_func

  # Init vector with length equal to number of sampling bouts
  # and with inital prior as its first value
  posteriors <- c(prior, rep(NA,length(data)-1))
  for(i in seq_along(data)) {
    bout = data[i]
    # option to bind posterior values to set a maximum confidence level
    # this allows users to make sure that the system stays agile enough to
    # change responses when data changes,
    # even if it has been confident for a long time.
    if(posteriors[i] < posterior_bnd_lower) posteriors[i] <- posterior_bnd_lower
    if(posteriors[i] > posterior_bnd_upper) posteriors[i] <- posterior_bnd_upper
    posteriors[i+1] = stbp_posterior(bout,
                                     hypotheses[i],
                                     likelihood_func,
                                     prior,
                                     lower_bnd,
                                     upper_bnd)
    # Break from iteration early if early_return is true,
    # the minimum iterations have been reached, and
    # if either of the decision criteria have been reached
    if(early_return &&
       i >= min_iterations &&
       ((posteriors[i+1] < lower_criterion) ||
        (posteriors[i+1] > upper_criterion))
      ) break
  }

  response <- if(posteriors[i+1] < lower_criterion) 1 else 0
  indices_above <- which(posteriors > upper_criterion)
  indices_below <- which(posteriors < lower_criterion)

  return(list(
    probabilities = posteriors,
    recommendation = response,
    num_iterations = i,
    decision_indices = list(
      above = indices_above,
      below = indices_below,
      first = min(indices_above, indices_below)
    )
  ))
}
