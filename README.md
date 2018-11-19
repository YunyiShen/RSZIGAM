This is a working repo for Yunyi's ecological modeling idea: Repeated Sampling Zero Inflated Generalized Additive Models. The idea is when there exist imperfect detection of the data, say counting number of birds, meanwhile, not all spots are occupied by certain birds, which cause zero inflate, we need to model theses factors at the same time.

The main function will follow R package [COZIGAM](https://www.jstatsoft.org/article/view/v035i11/0), by adding one more EM step to the main distribution (follow N-mixture models), to deal with the imperfect detection (count).

It equals to a zero inflated N-mixture model, but instead of using traditional linear regression strategy, here we use GAM version of all process including occupancy, latent lambda and detection probabilities. Though the point estimation of parameters comes from MPLE, the CI will adopt a Bayesian frame work, follow the Laplacian approximation suggested by COZIGAM paper.

 