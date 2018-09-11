# mu_mu_lambda-ES
Implement a $(\mu /\mu,\lambda )-ES$ 
It is initialized with $\mu$ parents, choose the centroid of parents to generate $\lambda$ offsprings and then pick the $\mu$ best offspring as the parents of next iteration. 
**Note**: $\lambda \geq \mu$ 

## Different noise level
Plot the performance of mml-ES under different noise-to-signal ratio defined as $v = \frac{\sigma_\epsilon^star}{\sigma^star} $ 
The undesired pattern for green dots is that the convergence rate is nan.
Under high $v$ and large $\sigma^*$ the algorithm does not converge, which gives the negative convergence rate.