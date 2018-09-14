# mu_mu_lambda-ES

Implement a (mu/mu,lambda)-ES

The algorithm is initialized with mu parents, choose the centroid of parents to generate lambda offsprings and then pick the mu best offspring as the parents of next iteration.

**Note**: lambda >= mu

## Different noise level

Plot the performance of mml-ES under different noise-to-signal ratio defined as v = sigma_ep_star/sigma_star.

The undesired pattern for green dots is that the convergence rate is nan.

Under high v (E.G. v=4) and large sigma_star the algorithm does not converge, which gives the negative convergence rate.



## Use CSA to adapt step size
### Algorithm
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

### Comparison
CMA does not seem to work in the context on the quadratic sphere. Compared with 0), CSA presumably decreases the step size too fast. So that the algorithm is trapped in some local point and therefore cannot move further.

### CMA is sensitive to initial step size

The performance is also very sensitive to the initial step size sigma0. 

If we set a large sigma0, the improvement in the first few steps is significant where the objective function value can be reduce to e^-2 in 10 runs. But the number of GP estimation in each iteration grows significantly. An example with sigma_ep_star = 4, sigma0 = 2, mu=3, lambda=10 gives the GP estimate in each iteration [0 0 0 0 460 1215 1706 4393 10751 755804 138400157] where the first 4 iterations uses the true objective function. 

For a small sigma0, the algorithm is trapped in some local point and therefore cannot move further as is in the case for 3) no matter small or large sigma0.

If we tune sigma0 to a proper value 0.5-0.6 we could achieve e^-4 within 20 iterations but the number of GP estimate also increases significantly as the number of iteration goes large. 


### Change of offspring selection criteria

Initially, we adapt the simple criteria by choose offsprings iff. f(centroid) > fep(candidate solution) where fep is the GP estimate of the candidate solution.

Using such criteria lead to the problem that even millions of sampling cannot obtain an offspring that matches the centroid. So we tried to loose the restriction on the offspring selected by adding a scalar sigma_GP i.e. the maximum tolerance that we can accpet for getting an offspring inferior to the centroid.

Another interpretation would be we have a centroid as the parent then we gennerate offspring in arbitrary direction. There must be a direction that leads to the optimal. What the offspring selection criteria does is to choose the area of that region that a candidate solution can be selected as an offspring. Once we increase the level of tolerance, the area that offspring can be obatined will increase and therefore the probabiliy of obtaining a candidate solution as an offspring where less sampling is needed.  

Now we choose an offspring iff. f(centroid)*sigma_GP > fep(candidate solution).

More tolerance (larger sigma_GP) would allow more progress (more steps), less sampling but less progress-to-step ratio i.e. less fitness gain for each step.

Presumably, there is a relation among sigma_GP and sigma_ep_star. It is understandable when there is more noise (variance) added to fy we need to put more selective on the offsprings. i.e. larger sigma_ep_star needs smaller sigma_GP.

#### Experiment
x0 = 
-0.203102904230873	-0.0856172663925321	1.30705750617960
-0.150834774366184	1.32683079014856	1.04826284432273
0.224821585068546	0.833055309940219	-0.763427819121924
-0.0629202204887581	-0.557425428225165	0.564949072093604
0.106706826487156	-0.599811245264293	0.0782480156340188
-1.54840125247811	-0.146604917055879	-0.340021761681992
-0.0408603991447191	1.87413826421482	-0.0702531036583463
0.287987538717003	-1.09111070759422	-0.358833142760595
-1.56273576980328	-0.531376791279362	1.39094457594521
-1.50993837819604	1.72474702674583	-0.702997115026487

| sigma0    | sigma_GP | sigma_ep_star| t(# of iterations)| # of GP estimates |
| :--------:|:--------:| :------:     | :------:          | :------:          |
| 1         |   3.5    |  1           | 55,73,64,55       | 231894,621,9789,65221588|


There is much randomness in the algorithm. Using the same parameter setting, the algorithm sometimes finishes with very limited sampling, sometimes finishes with over millions sampleing and sometimes just does not converge over 10000 iterations (almost half of the time).

Sampling a good candidate solution becomes much harder as we approaches the optimal, the number of sampling increases drastically in the last few steps before reaching the stopping criteria.

## Test CSA parameter

Verified in CSA, s = (1-c)*s + sqrt(mu*c*(2-c))*z which is n-dimensional. 

We plot |s|^2-n as the iteration goes. 

For noise-to-signal ratio v = 4 (quite high noise level) we obatin a similar pattern for sigma_star = 1,2 where the value goes up to around -9.7 then jump back to -10 and converge to -10. The scalar (exp((norm(s)-N)/(2*D*N)) mutiplied to sigma each time goes up to 0.216 then jump to 0.206. When sigma_star > 3, the algorithm does not converge and the both two terms are zero almost all time.

v = 2 similar pattern as stated.

---
## Change of Algorithm 

Usual (mu/mu,lambda)-ES replace the true objective function evaluation to candidate solutions generated by centroid with a Gaussian Process estimate, pick up the best mu candidate solutions as offsprings. Average the offsprings and get the new centroid (parent) for offspring generation in next iteration. Pseudocode is below. 

```
while not terminate() do{
	for j = 1,...,lambda{
		z(j) = randn(N,1)
		y(j) = x + sigma*z(j)
		fep(j) = evaluateUseGP(y(i))
	}
	z = sorted(z)
	z = 1/mu*sum(z(1),...,z(mu)))
	x = x + sigma*z
	update step size sigma
}
```

### 1. Model GP estimate using Gaussin distributed noise 

Assume the GP estimate can be modelled by the true objective function value with some random Gaussian noise terms.

fep(y) = f(y) + sigma_ep_star*randn(N,1)

### 2. Step size adapted by sigma = sigma_star/n*dist

The step size is proportional to the distance to the optimal point (dist).

### 3. Plot convergence rate over sigma* and noise-to-signal ratio

Plot is attached. With noise-to-signal ratio = 0,0.25,1,4. The negative convergence rate means the algorithm does not converge.

### 4. Step size adapted by cumulative step-size adaptation

The number of iteratins needed coincides that with the convergence rate plot, where a high convergence rate gives a small number of iterations to reach stopping criteria. 

### 5. Build GP model by adding training set 





