# nord_thyroid

This Repository shows R code and source data used in the publication *EARLY MORTALITY CRITICALLY IMPEDES IMPROVEMENTS IN THYROID CANCER SURVIVAL THROUGH A HALF CENTURY*, which is currently under review in the [*European Journal of Endocrinology*](https://academic.oup.com/ejendo)

All data used were downloaded from the [NORDCAN website](https://nordcan.iarc.fr/en). 

For the estimation of non-linear trends in survival, we used the NORDCAN data that include the estimates of 1-year and 5-years survival and their confidence intervals, based on the Pohar Perme estimator. This code utilizes the pre-calculated NORDCAN data (including the pre-calculated confidence intervals) to get continuous survival estimates and estimation of the magnitude of the change (slope) across a continuous timescale. The code also estimates 95% credible intervals for both the survival rate and the magnitude of change, taking into account all data. This is done via Generalized Additive Models in the Bayesian framework  with [brms package](https://cran.r-project.org/web/packages/brms/index.html)

Please, cite the [original article when using this R code

### Full citation:

>[will be added]^

### Statistical methods

Statistical modelling and all data visualizations were performed using [R statistical software](https://www.r-project.org/) in the [R studio](https://posit.co/download/rstudio-desktop/) environment

Time trends of 1-year and 5-years relative survival (in %; obtained from NORDCAN for each of the 5-year periods) were modelled via Gaussian generalized additive models (GAM) with thin plate splines (5 knots) and identity link. The GAM model included the effect of group (combination of sex and country) and group-specific non-linear effect of time (timepoint = middle year of each 5-years period) as predictors, allowing estimation of the relative survival across a continuous time scale despite the discrete distribution of data points. As the input data (estimates of the 1-year and 5-years survival in each of the 5-year periods) were variously uncertain (as expressed with confidence intervals calculated by NORDCAN), standard errors of estimates were accommodated into models. Models were run in the Bayesian framework using the [‘brms’ R package](https://cran.r-project.org/web/packages/brms/index.html) which employs [‘Stan’ software](https://mc-stan.org/) for probabilistic sampling. Separate models were used for 1-year and 5-years survival.

The prior distribution for the effect of the country was explicitly defined to have Gaussian distribution with zero mean and sigma of 30. Default brms priors were used for other parameters. We used Hamiltonian Monte Carlo sampling (2 chains, each of 7,000 samples including 2,000 warm-ups). All models were checked in terms of convergence, effective sample sizes and posterior predictive check. 

For 5/1-year survival ratio estimation, we divided posterior draws from the 5-year survival model by posterior draws from the 1-year model to get the posterior distribution of the conditional survival and its estimated annual changes over time. 

For all survival measures (relative 1-year and 5-years survival and 5/1-year ratio), we evaluated when the survival was changing over time with at least 95% plausibility (95% credible interval [CI] of the 1st derivation of given survival measure did not cross zero for at least 5 years). We also aimed to identify ‘breaking points’, i.e. times when the annual change of survival changed with at least 95% plausibility. This was assessed by calculation of the 2nd derivation of the given survival measure and its 95% CI; the ‘breaking point’ was defined as a peak value within at least a 3-year interval where 95% CI for the 2nd derivation did not cross zero.
