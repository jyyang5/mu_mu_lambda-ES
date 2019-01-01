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

- 20181220
	- [x] Scatter plot (\# of iterations vs. success rate, one dot per run using one $\sigma^* $)

- 20181223
	- [x] Use success-rate-based step size adaptation for mml-ES with S=0.5:0.1:0.9 

- 20181230
	- [x] Use success-rate-based step size adaptation for (1+1)-ES with S=0.1:0.1:0.6 





### 1. Step size adaptation 

#### 1.1 Normalized step size 

Use fixed normalized step size to update step size.	

#### 1.2 Success rate 

Adapt step size depending on the success rate.

- mml-ES

	- Increase S (large success rate 0.6-0.9)
		- better performance for the last two test functions 
		- worse on sphere functions


- (1+1)-ES

	- reduce S (small success rate 0.6-0.9)
		- better performance for sphere functions
		- worse for the other two