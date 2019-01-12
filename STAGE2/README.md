## STAGE 2
- Faster speed-up

- Simple and elegant step size adaptation mechanism 


## Files

- [convergencePlot_success]
	- Step size: adapted by normalized step size (sigmaStar)
	- Plot success rate

- [scatter_iteration_success, scatter_iteration_success_large_scale]
	- Scatter plot: number of iteration vs. success rate (for each lambda and each test function)
	- Step size: adapted by normalized step size (sigmaStar)
	- Large_scale: more replicates 

---

- [mml-different_success_rate,step_size_adapatation_success_rate]
	- First attempt mml-GP (Comma-selection )
	- Step size: adapted based on success rate S

- [(1+1)_different_success_rate]
	- Strategy
	    - Train: mml-ES 
    	- After: bestSoFar (over 2 generations)

- [cross-different_success_rate]: 
	- Strategy
	    - Train: mml-ES 
    	- After: GP-mml-ES  (One objective function call per generation for centroid (obtained from progress vector))
	- Comma-selection 
	- Step size: adapted based on success rate S

- [bestOfTwoIterations_plus-selection]: plus-selection (1+1)-ES each centroid lifespan = 2 generations 
	- Strategy
	    - Train: mml-ES 
    	- After: BestSoFar (One objective function call per generation for centroid (obtained from progress vector))
	- Plus-selection (over two generations)
	- Step size: adapted based on success rate S	

- [bestSoFar_plus-selection]: plus-selection (1+1)-ES variant
	- Strategy
    	- Train: mml-ES 
    	- After: BestSoFar (One objective function call per generation for centroid (obtained from progress vector))
	- Plus-selection 
	- Step size: adapted based on success rate S

---
- [two_level_adaptation]: attempt to update step using two levels of success rate 

- [fourProb_sigmaStar]
	- Plot four probs of GP: TP,TN,FP,FN
	- Step size: adapted by **normalized step size (sigmaStar)**

- [fourProb_prob]
	- Plot four probs of GP: TP,TN,FP,FN
	- Step size: adapted based on **success rate S**


- [arash_variant]
    - Strategy
    	- Train: mml-ES 
    	- After: BestSoFar 
	- Step size: adapted based on S1,S2,S3

- [Dirk]
	- Strategy
    	- Train: (1+1)-ES 
    	- After: BestSoFar 
	- Step size: adapted based on S1,S2,S3 (S1=S2=0.8,S3=0.1)

- [replicate_1+1_arash_variant]
	-- Strategy
    	- Train: (1+1)-ES 
    	- After: BestSoFar 
	- Step size: adapted based on S1,S2,S3 
	- Objectives 


		1. plot change S3: invariant to changes from S1=S2=0.8, S3=0.05*2.^(0:1:3) [fun_multiPlus_change_C3_FUNCALLS.m AND run_mutiPlus_C3.m]
		2. plot change over S1 and S2 based on success rate S [fun_multiPlus_change_successRate_FUNCALLS.m AND run_mutiPlus_successRate.m]
			- S1 = S, S2 = 1-S
			- Where S1/(S1+S2) \approx success rate [fun_multiPlus_change_successRate_FUNCALLS.m AND run_mutiPlus_DF.m]
		3. plot change over DF (damping factor) 
			- S1 = S*DF, S2 = (1-S)*DF
	- Figures plots
		- histogram of objFunCalls
		- convergence plot 
		- step size 
		- normalized step size 
		- success rate and evaluation rate 
		- four probs TP,TN,FP,FN




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


- [20190110]
	- [x] Variant of what Arash did. Use C1 = 0.05, C2 = 0.6, C3 = 0.6
	   - [x] Length scale = 8, does not match that of what Dirk did 

	   	Because the first 40 generations uses (1+1)-ES
	   - Improvements when length scale = 16




- [20190111]
	- [x] First 40 generations using (1+1)-ES

	- [x] Replicate the results 
		- Test functions: condense defs (compare results)
	- [] C1, C2, C3

		Success C1/C2
		equal 1/2
		
		C1 = 0.2, C2 = 0.8 (1/5-success rate)

		C3: safeguard -> accuracy of the model (always get bad solutions)

		if C1 = C2 (1/2-success rate)



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


#### 1.3 Using prob. with four cases 

- [Arash's variant]
```matlab
if fep_y < fx
	sigma = sigma*-C3/D
else
	if fy < fx
		sigma = sigma*-C1/D
	else
		sigma = sigma*C2/D
	end 
end

```


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






