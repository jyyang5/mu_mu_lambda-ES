## STAGE 2
- Faster speed-up

- Simple and elegant step size adaptation mechanism 

### Schedule 
- 20181218
	- [x] Normalized step size (what we can achieve, as upper bound)
	- [x] If the opt. step size is independent of lambda
	- [] \*Can try weighted recombination but weights do not necessarily sum to 1 (ref Foundations_of_Genetic_Algorithms) 

- 20181219
	- [x] Corresponding success rate 
		- Success rate for mml-ES is significantly higher than that of GP-mml-ES
		One explanation is that the steps are much cheaper (one function evaluation rather than lambda per generation)

- [20181220](https://github.com/jyyang5/mu_mu_lambda-ES/tree/STAGE2/STAGE2/scatter_iteration_success)
	- [x] Scatter plot (\# of iterations vs. success rate, one dot per run using one $\sigma^* $)

- [20181223](https://github.com/jyyang5/mu_mu_lambda-ES/tree/STAGE2/STAGE2/mml-different_success_rate)
	- [x] Use success-rate-based step size adaptation for mml-ES with S=0.5:0.1:0.9 

- [20181230](https://github.com/jyyang5/mu_mu_lambda-ES/tree/STAGE2/STAGE2/(1%2B1)_different_success_rate) 
	- [x] Use success-rate-based step size adaptation for (1+1)-ES with S=0.1:0.1:0.6 





### 1. Step size adaptation 

#### 1.1 Normalized step size 

Use fixed normalized step size to update step size.	

#### 1.2 Success rate 

Adapt step size depending on the success rate.

- [mml-ES](https://github.com/jyyang5/mu_mu_lambda-ES/tree/STAGE2/STAGE2/mml-different_success_rate)
	- Increase S (large success rate 0.6-0.9)
		- better performance for the last two test functions 
		- worse on sphere functions


- [(1+1)-ES](https://github.com/jyyang5/mu_mu_lambda-ES/tree/STAGE2/STAGE2/(1%2B1)_different_success_rate) 
	- reduce S (small success rate 0.6-0.9)
		- better performance for sphere functions
		- worse for the other two


### Strategies used

- [cross-ES]
	- Evaluate centroid using true objective function in each generation  

- [plus-selection (bestOfTwoIterations))]
	- Plus-selection over two generations (the centroid can survive another generation) i.e., selection based on new centroid and the last centroid from previous generation 

- [plus-selection (bestSoFar)]
	- Plus-selection over all generations

- [variants (fourCases)]
	- Compare GP estimate with previous centroid (similar to surrogate model assisted (1+1)-ES)


#### 1.3 Prob. of four cases 

fep_centroid_temp: GP estimate of the new centroid 
f_centroid_temp: true objective function value of the new centroid 
f_centroid: true objective function value of best centroid so far


- fep_centroid_temp > f_centroid && f_centroid_temp > f_centroid

- fep_centroid_temp > f_centroid && f_centroid_temp < f_centroid

- fep_centroid_temp < f_centroid && f_centroid_temp > f_centroid

- fep_centroid_temp < f_centroid && f_centroid_temp < f_centroid




