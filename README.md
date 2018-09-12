# mu_mu_lambda-ES

Implement a (mu/mu,lambda)-ES

The algorithm is initialized with mu parents, choose the centroid of parents to generate lambda offsprings and then pick the mu best offspring as the parents of next iteration.

**Note**: lambda >= mu

## Different noise level

Plot the performance of mml-ES under different noise-to-signal ratio defined as v = sigma_ep_star/sigma_star.

The undesired pattern for green dots is that the convergence rate is nan.

Under high v (E.G. v=4) and large sigma_star the algorithm does not converge, which gives the negative convergence rate.



## Use CSA to adapt step size

Replace the previous step size adaptation using sigma = sigma_star/n*dist (denoted appraoch 0.) with Cumulative Step Length Adaptation (CSA) from https://pdfs.semanticscholar.org/2556/e4670234fc74232a1a419260fe8c1a29d0ba.pdf.
```
c = 1/sqrt(N)
D = sqrt(N)
s = 0
while not terminate() do{
	for j = 1,...,lambda{
		z(j) = randn(N,1)
		y(j) = x + sigma*z(j)
		f(j) = evaluate(y(i))
	}
	z = sorted(z)
	z = 1/mu*sum(z(1),...,z(mu)))
	x = x + sigma*z
	s = (1-c)*s + sqrt(mu*c*(2-c))*z
	sigma = sigma*exp((norm(s)-N)/(2*D*N)
}
```
Three approaches are used for CSA: 
1. Average each dim of N-dimensional Z
   z(j) = randn(N,1)
   s = (1-c)*s + sqrt(mu*c*(2-c))*mean(z)

2. N-dimensional Z dealt with by L2-norm 
   z(j) = randn(N,1)
   s = (1-c)*s + sqrt(mu*c*(2-c))*z

3. Use 1-dimensional Z lacking randomness 
   z(j) = randn()
   s = (1-c)*s + sqrt(mu*c*(2-c))*z

##### Comparison
CMA does not seem to work in the context on the quadratic sphere. Compared with 0), CSA presumably decreases the step size too fast. So that the algorithm is trapped in some local point and therefore cannot move further.

##### CMA is sensitive to initial step size

The performance is also very sensitive to the initial step size sigma0. 

If we set a large sigma0, the improvement in the first few steps is significant where the objective function value can be reduce to e^-2 in 10 runs. But the number of GP estimation in each iteration grows significantly. An example with sigma0 = 2, mu=3, lambda=10 gives the GP estimate in each iteration [0 0 0 0 460 1215 1706 4393 10751 755804 138400157] where the first 4 iterations uses the true objective function. 

For a small sigma0, the algorithm is trapped in some local point and therefore cannot move further as is in the case for 3) no matter small or large sigma0.

If we tune sigma0 to a proper value 0.5-0.6 we could achieve e^-4 within 20 iterations but the number of GP estimate also increases significantly as the number of iteration goes large. 

