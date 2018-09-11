# mu_mu_lambda-ES
Implement a (mu/mu,lambda)-ES 

The algorithm is initialized with mu parents, choose the centroid of parents to generate lambda offsprings and then pick the mu best offspring as the parents of next iteration. 

**Note**: lambda >= mu

## Different noise level
Plot the performance of mml-ES under different noise-to-signal ratio defined as v = sigma_ep_star/sigma_star. 

The undesired pattern for green dots is that the convergence rate is nan.

Under high v and large sigma_star the algorithm does not converge, which gives the negative convergence rate.