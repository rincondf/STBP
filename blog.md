# Simple, intuitive procedure for sequential inferences
Post provided by Izzy McCabe

My colleague Diego Rincon and I have just published "Sequential testing of complementary hypotheses about population density" in MEE.
This paper started as an offshoot of our research into improving site-specific predictions of Codling moth, a prolific tree fruit pest.
The state of the industry for integrated pest management (IPM) typically involves a phenology model, which relates heat accumulation to
population dynamics, in combination with on-site sampling to inform decisions about controlling pest populations. A drawback of existing
methods is that phenology and sampling remain unlinked; growers use their intuition and experiences to connect the dots. The
intuition of experienced growers is not to be underestimated, and they are skilled at knowing how pest counts by certain periods in
the season are related to seasons where management interventions are needed. We wanted growers to be able to set maximum tolerable thresholds,
and for us to use our site-specific profiles in conjunction with early-season data to predict the probability of exceeding that threshold.
After surveying the state-of-the-art in sequential statistical procedures, we found existing methods disappointing. The most widely represented
in the literature is the Sequential Probability Ratio Test (SPRT), introduced in 1945 by Abraham Wald, but it requires fixed sample sizes and
the specification of two non-complementary models (i.e. the case where average captures are 8 vs 10, rather than the cases of being greater than
or less than 9). Furthermore, the original methodology only accounts for one-at-a-time sequential sampling, and requires complicated modifications
to support group sequential designs that are inherent to pest sampling networks. I believed that the problem of iteratively updating our knowledge
and predictions based on new data aligned with the philosophy of Bayesian statistics, and sure enough we found methods utilizing Bayes' Rule to
the domain of sequential hypothesis testing in Morgan & Cressie's Variable Probability Ratio Test (VPRT). Their method, while providing compelling
optimality properties, still requires non-composite hypothesis ratios and comes with the computational difficulties of Bayesian methods. Our new method,
the Sequential Test of Bayesian Posterior Probabilities (STBP) can be seen as building on the VPRT, and aims to provide an ergonomic and intuitive procedure
that explicitly supports variable sample sizes over time, threshold-style hypotheses, and even dynamic trajectories of hypotheses.

It works by making a crucial simplification of typical Bayesian methods: instead of the specification of proper prior probability distributions, it
assigns a single probability to the possibility of your parameter of inference being greater than the threshold and the complementary probability for the
possibility of it being less than that threshold. This is comparison to the typical procedure of giving each possible value of the parameter its own
probability density. This reduces the resulting posterior probability from an entire distribution to a single number, which we can then plug back into
the equation as the prior probability for the next iteration. It also eliminates the difficult computational exercise of integrating over constantly changing
joint probability density functions which forms the major limitation for typical Bayesian methods. Integrals only need be calculated over your desired likelihood
function, for which all commonly used distributions already have closed-form solutions in the form of their cumulative density functions! While our method had a
greater rate of false-positives compared to the SPRT in purely sequential designs, it outperformed it in group sequential designs (FIGURE HERE), requiring less sampling and
errors scaled more favorably with increasing sample sizes than the SPRT. We also tested our method as a procedure for invasive species monitoring, where it
demonstrated higher power for equal sample sizes than conventional approaches, with further reductions in required sampling when prior belief in the absence of
the invasive species is already high (FIGURE HERE). We also found that our simplified model quickly converges on the results produced by more complex Bayesian techniques
requiring difficult integrations. (INSERT FIGURE THAT WE SENT TO THE REVIEWER OF THE MODELS CONVERGING) 

Because our method allows for both variable sample sizes and sequential comparisons for changing-but-related hypotheses, it also has potential applications
in continuous monitoring for signal processing, fraud detection, forensics, and healthcare contexts (like that of changes in glucose levels or heart rate).
For example, by continuously comparing the average patient mortality rates for certain surgical procedures to
that of individual surgeons, we found that Harold Shipman's serial killings could have been detected with great confidence as early as (INSERT DATE AND FIGURE).

Our method represents innovation in a section of the ecological and agricultural literature which has remained largely stagnant for nearly many decades.
The SPRT, introduced in 1945, is still the most widely used method for sequential analysis, and the fixed-sample size method used in invasive species
monitoring has been used unchanged since its introduction in 1993. We plan on working to integrate the STBP into integrated pest management strategies
in the Pacific Northwest, and hope to see others find places in their corners of science and industry where it can improve their processes and research.
To facilitate this, we are working on an R package that will be easy to use and extensible.
The source code can be found in [this GitHub repository](https://github.com/PacificBird/STBP), and we will be working to release an R package to the CRAN
repository in the coming weeks.

