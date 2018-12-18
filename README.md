# RAZIGAM
## Introduction
This is a working repo for Yunyi's ecological modeling idea: Repeated Sampling Zero Inflated Generalized Additive Models. The idea is when there exist imperfect detection of the data, say counting number of birds/presence of tigers, meanwhile, not all spots are occupied by certain birds, which cause zero inflate, we need to model theses factors at the same time. BBS data will be used to test the model.

## The Model
The main function ZIGAM.R consists of 2 functional parts. A GAM version Zero-Inflated N-mixture Model and a GAM version of single species Occupancy model. The main algorithm generally follows the EM and modified PIRLS method suggested by R package [COZIGAM](https://www.jstatsoft.org/article/view/v035i11/0) , (Hai Liu and Kung-Sik Chan 2010). The EM-PIRLS will iterate until converge.

For the zero inflated N-mixture model, instead of using traditional linear regression strategy, here we use GAM version of all process including occupancy, latent lambda and detection probabilities. By adding one more EM step for the main Poisson-Binomial distribution we can deal with the missing observation of occupancy and latent n at the same time. Three GAMs for psi, lambda and p are fitted with "pseudo-data" generated from Estep of the EM algorithm  and coupled usint the modified PIRLS. Though the point estimation of parameters comes from MPLE, the CI adopts a Bayesian frame work, follow the Laplacian approximation suggested by COZIGAM paper.

For the GAM-occupancy model, since the process is much simpler (no latent lambda), the algorithm of COZIGAM can be directly used to obtain the point estimation of parameters, but CI also adopts the Bayesian-Laplacian approximation approach.

## Simulation Study
### GAM Zero-Inflated N-Mixture
Simulation have been done for a simple additive case without interaction terms, two environmental factors and one detection factor were simulated. Psi were set to be linear with env.1 and env.2, lambda has a bimodal relation to env.1 and unimodal relation with env.2, detection rate p is linear with env.1 and det.1 while no relation with env.2. Same setting with 100, 400 and 900 sites, 10, 20 and 40 repeats were tested.  Difference between fit and real value was evaluated using relative 1-norm difference. Reaction curves were drew to compare the fit curves and real ones.


## Primary Results
### GAM Zero-Inflated N-Mixture
Results showed that the model is a data-hungry method which ask for 400 or more sites and 10 repeats. For 100 sites EM-PIRLS even did not converge. Relative difference according to 1-norm was around 10% for psi, 15% for lambda and 3% for p in the 900 sites simulation. Detailed "landscape" evaluation can be found in Fig.4. Compare of reaction curve can be seen in Fig.1~3. There can be a scaling factor different for lambda and constant different for psi and p.

![Fig.1 Reaction cure of psi, left is fitted](https://raw.githubusercontent.com/YunyiShen/RSZIGAM/master/Doc/figs/psi.png)
Fig.1 psi
![Fig.2 Reaction cure of lambda, left is fitted](https://raw.githubusercontent.com/YunyiShen/RSZIGAM/master/Doc/figs/lambda.png)
Fig.2 lambda
![Fig.3 Reaction cure of p, left is fitted](https://raw.githubusercontent.com/YunyiShen/RSZIGAM/master/Doc/figs/p.png)
Fig.3 p
![Fig.4 Lambda and Psi across the "landscape"](https://raw.githubusercontent.com/YunyiShen/RSZIGAM/master/Doc/figs/p.png)
Fit.4 landscape



