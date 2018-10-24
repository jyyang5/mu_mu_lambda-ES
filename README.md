# A mu_mu_lambda-ES using comparision-based surrogate(GP) 

## Overview
- Folders
    - [CSA](https://github.com/jyyang5/mu_mu_lambda-ES/tree/master/CSA)

   	    Implement a CSA for small lambda
   	    - Step size adaptation: CSA
   	    - Offspring evaluation: GP
   	    - Parent selection strategy: centroid of best mu offspring
    - [CSA_largeLambda](https://github.com/jyyang5/mu_mu_lambda-ES/tree/master/CSA_largeLambda)

	    Implement a CSA for large lambda **to be continued**
   	    - Step size adaptation: CSA
   	    - Offspring evaluation: GP
   	    - Parent selection strategy: centroid of best mu offspring
    - [GP](https://github.com/jyyang5/mu_mu_lambda-ES/tree/master/GP)

        Plot fep(x) and f(x) over number of iterations
   	    - Step size adaptation: CSA
   	    - Offspring evaluation: GP
   	    - Parent selection strategy: centroid of best mu offspring

    - [choose_centroid_bestSoFar](https://github.com/jyyang5/mu_mu_lambda-ES/tree/master/choose_centroid_bestSoFar)

        Plot fep(x) and f(x) over number of iterations
   	    - Step size adaptation: CSA
   	    - Offspring evaluation: GP
   	    - Parent selection strategy: best centroid of all iterations

    - [compare_result](https://github.com/jyyang5/mu_mu_lambda-ES/tree/master/compare_result)
        - Plot the following over objective function evaluations 
          - objective function value
          - step size
          - normalized step size
          - relative model error

        - Plot histogram for the following 
          - objective function evaluations
          - convergence rate
          - success rate 

        - Previous attempts
          
          Use a combination of mu and lambda plot f(x), sigma and sigmaStar for 5 test functions 
   	      - Step size adaptation: CSA
   	      - Offspring evaluation: GP
   	      - Parent selection strategy: centroid of best mu offspring

    - [convergence rate_noise](https://github.com/jyyang5/mu_mu_lambda-ES/tree/master/convergence%20rate_noise)

   	    Plot convergence rate for different noise to signal ratio and dimensionality with expected convergence rate plotted
   	    Code to intergrate the expected convergence rate
        - Step size adaptation: sigmaStar
   	    - Offspring evaluation: noise estimate 
   	    - Parent selection strategy: centroid of best mu offspring

    - [convergence_sigmaStar_GP](https://github.com/jyyang5/mu_mu_lambda-ES/tree/master/convergence_sigmaStar_GP)

        Use different strategies and plot convergence rate over sigmaStar
        **Conclusion: Choose centroid is the bset esp. with large lambda**
        - Step size adaptation: sigmaStar
   	    - Offspring evaluation: GP
   	    - Parent selection strategy: 
   		    1. centroid (centroid of the best mu offsprings ranked by GP)
   		    2. bestOfMu (best in the mu offspring ranked by GP)
   		    3. bestSoFar (bestOfMu in each iteration, replace y iff. f(current) < f(previous))
   		    4. bestCentroid (replace previous centroid iff. f(current) < f(previous))
   		    5. centroidAddTrain=1 (centroid as parent, add the best offspring ranked by GP to training set per iteration)
   		    6. centroidQuadratic_withGP (fit a quadratic model for z_centroid, choose lowest point, pos and neg evaluated by GP)
   		    7. centroidQuadratic_trueObjFun (fit a quadratic model for z_centroid, choose lowest point, pos and neg evaluated by true objective function)

   		Change the distance norm in GP use l-0.5 norm and l-3 norm **failed**

    - [mu_sigma_experiment](https://github.com/jyyang5/mu_mu_lambda-ES/tree/master/mu_sigma_experiment)

   	    Tuned TRAINING_SIZE with sigmaStar/CSA and did a couple of result.
   	    **Possible solutions for better performance: 1. update GP model 2. tune training size of GP model**
   	    **Best result using (3/3,10)-ES with sigmaStar = 1,2 relative speed-up 1.7,1.3,1.11 for linear, quadratic and cubic sphere respectively (sigmaStar = 1).**
   	    Performance degrades as sigmaStar grows, but large sigame appears to be more invariant to the growth of sigmaStar
   	    - tuneTrainSize_sigmaStar_GP
   	    - Step size adaptation: sigmaStar
   	    - tuneTrainSize_CSA_GP 
   	    - Step size adaptation: CSA

   - [convergence_lengthScale_GP](https://github.com/jyyang5/mu_mu_lambda-ES/tree/master/convergence_lengthScale_GP)

   	    Plot convergence rate over theta (fun_c_over_theta.m and run_c_over_theta.m) for sigmaStar=1:0.5:3 lambda = [30,50,70]
   	    **max convergence rate relative the same for (7/7,30)-ES, relative 0.1 increment for (12/12,50)-ES and (17/17,70)-ES**

| max convergence rate  |theta   | linear(sigmaStar)  | quadratic(sigmaStar)  | cubic(sigmaStar)  |
| :---------------------|-------:| ------------------:|--------------:|------------:|
| **(7/7,30)**  		|**8**   |0.48                |1.02		 |0.811|  
|   (7/7,30)    		|  8     |0.52(1.5)           |          |     |
|   (7/7,30)    		|  32    |                    |1.18(2.0) |     |
|   (7/7,30)    		|  2     |                    |          |0.89(1.5,2.0)   |  
| **(12/12,50)**		|**8**   | 0.53               |1.10		 |0.808|
|   (12/12,50)    		|  32    | 0.57(1.5) |||
|   (12/12,50)    		|  26    | |1.2(2.0,2.5) ||
|   (12/12,50)    		|  2     | ||0.96 (2.5)|
| **(17/17,70)**		|**8**   | 0.54   |1.07		 |0.799
|   (12/12,50)		    | 16     | 0.62(1.5)   |||
|   (12/12,50)		    | 28     | |1.25(2.0)   ||
|   (12/12,50)		    | 2      | ||0.96(3.0)   |

## 	Schedule 

- 20181011

	 - [x] Use GP sigmaStar for (3/3,10),(5/5,20), (10,40) for the plots below 
	    1. expected fitness gain eta over normalied step size sigmaStar
		2. opt. normalized step size over noise-to-signal ratio 
		3. opt. expected fitness gain over noise-to-signal ratio
     - Thought 
     	- Use sigmaStar test one sphere functions for optimal sigmaStar (largest convergence rate)
     	- Adapt step size s.t. the normalized step size close to opt. sigmaStar previously 
     	- Achieve appriximate similar convergence rate
     - Comment
     	1. scatter too close to each other 
     	2. strage n=100 has larger fitness gain compared with n=10 (in opt. normalized step size)
     	3. y-axies not correct
     	4. problem if 


- 20181012
	
	- [x] Get better approximate of progress rate using formula in paper when n = infty 

		- Use formula (11) from [paper](https://core.ac.uk/download/pdf/81976199.pdf) denoted n \approx \infty


	- [] Replace Gaussin distributed noise with GP 

- 20181016

	- [] Gap between the expected progress rate and the one using GP

		- TRAINING_SIZE
		- length scale

	- [x] Twice almost twice for the the opt. step size for experimental result and expected value

		- Fixed problem lies in the sigmaStar used for ploting

	- [] Performance for large dimensionality
		- Expect at least linear to the increase of dim

- 20181018/20181019

	- [] Thesis write-up
                3. Plots 
                        - Expected progress rate over \sigmaStar
                        - opt. normalized step size and opt. expected fitness gain over noise-to-signal ratio 
                        - step size 
                4. Step size adaptation (use GP)
                        - Normalized step size over iteration 
                        - Noise to signal ratio (var(fep-fy)/var(fy) over iteration) should be a straight line
                                
                all lambda offspring evaluated by GP (fep) and true objective function value (fy)
                
                calculate the noise-to signal ratio per iteration as var(fep,fy)/var(fy)
  - [] Draft for thesis 
    - Figure and tables are fine 
    - Little text 

  - [] Step size adaptation 
    - Use the standard [Niko's](http://www.cmap.polytechnique.fr/~nikolaus.hansen/handout2006.pdf)
    - Compare with opt. normalized step size using sigmaStar

	- [] TRAINING_SIZE (minor)

        3. Plots 
            - Expected progress rate over \sigmaStar (lambda=10,20,40 n=10 include expected and (1+1)-ES)
            - opt. normalized step size and opt. expected fitness gain over noise-to-signal ratio 
            - step size 
        4. Step size adaptation (use GP)
            - Normalized step size over iteration (plot)
            - Noise to signal ratio (var(fep-fy)/var(fy) over iteration) should be a straight line
                                
                <all lambda offspring evaluated by GP (fep) and true objective function value (fy)>
                
                <calculate the noise-to signal ratio per iteration as var(fep,fy)/var(fy)>
                
                <20,40,80 plot same -> 40 is the best (stick to n =10) too much information> 

                <focus on the same plot include the (1+1)-ES with GP>

                <variance (noise-to-signal ratio) over iteration> 

                <var(fye-fy)/var(fy) fitness in the population> 



	- [] TRAINING_SIZE 
    step_size_adaptation

		- Try TRAINING_SIZE = 4*n (round up to k*lambda)???

- 20181022
    - [] Step size adaptation 
       - Use the standard [Niko's](http://www.cmap.polytechnique.fr/~nikolaus.hansen/handout2006.pdf)

    - [x] compare result for plot 
       - Histogram 
       -  
        
    - Obserations (Use CSA from Dirk's lecture slides **wrong**)
      1. mml-ES seems to help avoid the potential poor direction resulted from an inaccurate GP model. Since the success rate is always around 0.48 as is expected while the success rate for  (1+1)-ES varies among test functions (from hist_xxx_sphere.fig)

      2. it is likely that mml-ES is more sensitive to relative model error which has a direct impact on the convergence rate. We can achieve similar performance to (1+1)-ES if the relative model of mml-ES is at least half smaller than that of (1+1)-ES (see xxx_sphere_funCalls.fig in linear sphere the model error for mml-ES is larger which potentially leads to the small convergence rate and therefore more objective evaluations).  

      3. one potential way out is to adapt the length scale factor (probably a length scale adaptation mechanism).The length scale factor 8 seems to be highly biased on quadratic sphere.  But from previous experiment I have no idea what the relation might be. I thought it could be related to the complexity of the test function. 

		
- 20181023
    - Success rate (48% for mml-ES)
        - Gnenerate one direction
        - (1+1) generate one direction (take step only when good)
        - mml generate the gradient (always make the step)
        - Step is bad every other iteration 
        - Reject the step
    - Rejecting steps (cannot be slower how much gives)

        - Every other step is bad
        - Could try what Arash did in step size adaptation

    - Impact of population(lambda) on noise-to-signal ratio(model error)

        how error depends on lambda (how model error depend on lambda)
        - Plot model error over lambda

    - Success rate use plot Histogram improvement(pdf) 
        Difference of current centroid compared with previous centroid 
        divide by R* derivative (R is the distance to optimal)
          - linear: R 
          - quadratic: 2R^2 
          - cubic: 3*R^3 

    - CSA
        - Just change the coefficients 
        - Try the latest one see if improvements






    








		
 





---
# mu_mu_lambda-ES

Implement a (mu/mu,lambda)-ES

The algorithm is initialized with mu parents, choose the centroid of parents to generate lambda offsprings and then pick the mu best offspring as the parents of next iteration.

**Note**: lambda >= mu

## Schedule

### 1. Model GP estimate using Gaussin distributed noise 

Assume the GP estimate can be modelled by the true objective function value with some random Gaussian noise terms.

fep(y) = f(y) + sigma_ep_star*randn(N,1)

### 2. Step size adapted by sigma = sigma_star/n*dist

The step size is proportional to the distance to the optimal point (dist).

### 3. Plot convergence rate over sigma* and noise-to-signal ratio

Plot is attached. With noise-to-signal ratio = 0,0.25,1,4. The negative convergence rate means the algorithm does not converge.

### 4. Step size adapted by cumulative step-size adaptation

The number of iteratins needed coincides that with the convergence rate plot, where a high convergence rate gives a small number of iterations to reach stopping criteria. 

- Convergence rate and normalized step size over noise ratio


### 5. Build GP model by adding training set 

The objective function evaluation value of offsprings is evaluated using the GP model after the model is built and therefore, reduce the objective function evaluation time to 1 in each iteration.

The improvement compared with a normal (mu/mu,lambda)-ES is the saving of objective function evaluation in each itertaion. The effectiveness of the improvement lies in the precision of the approximation using Gaussian Process model.


### 6. Compare the mml-ES with GP and without and similarly with (1+1)-ES 


The objective function evaluation obtained by averaging 400 runs is shown below. The plot of step size $\sigma$ and objective function value f(x) over # of objective function calls for 5 different test functions are attached. 

	(1+1)-speed-up = (# of objective function calls of (1+1)-ES withGP)/(# of … of (1+1)-ES noGP)
	mml-speed-up = (# of objective function calls of mml withGP)/(# of … of mml noGP)
	(1+1)-mml-ratio = (# of … of (1+1)-ES noGP)/(# of objective function calls of mml noGP)
	(1+1)-mml-GP-ratio = (# of … of (1+1)-ES withGP)/(# of objective function calls of mml withGP)
	
### 7. The precision of GP estimate

#### Compare funCalls (2/2,5)-ES lengthScale=16

Number of runs = 500
Length scale theta = **16** *sigma*sqrt(n)
(2/2,5)-ES
Evaluate and choose centroid as parent for each itertation
**Note**: plot and data stored in folder plot_cmp_result3

| Test function       |mml-ES withGP|mml-ES noGP|(1+1)-ES withGP|(1+1)-ES noGP|mml-speed-up|(1+1)-speed-up|mml-(1+1)-GP-ratio| mml-(1+1)-ratio|
| :-------------------|------------:| ---------:|--------------:|------------:|-----------:|-------------:|-----------------:|---------------:|
| linear sphere      |1128|4460|505|1279|      
| quadratic sphere   | 479|2350|214|676|   
| cubic sphere       | 386|1640|203|483|
| Schwefel’s function|1833|13825|1246|4233|  
| quartic function   |3356|11320|1502|2385|  

**Observation: the length scale factor does not seem to make a difference.**

#### Compare funCalls (2/2,5)-ES 

Number of runs = 500
Length scale theta = 8*sigma*sqrt(n)
(2/2,5)-ES
Evaluate and choose centroid as parent for each itertation
**Note**: plot and data stored in folder plot_cmp_result2

| Test function       |mml-ES withGP|mml-ES noGP|(1+1)-ES withGP|(1+1)-ES noGP|mml-speed-up|(1+1)-speed-up|mml-(1+1)-GP-ratio| mml-(1+1)-ratio|
| :-------------------|------------:| ---------:|--------------:|------------:|-----------:|-------------:|-----------------:|---------------:|
| linear sphere      |1097|4470|507|1271|      
| quadratic sphere   | 485|2340|213|678|   
| cubic sphere       | 385|1660|202|480|
| Schwefel’s function|2681|13825|1246|4233|  
| quartic function   |3283|11285|1498|2378|  

#### Compare funCalls (3/3,10)-ES

Number of runs = 400
Length scale: theta = 8*sigma*sqrt(n)
(3/3,10)-ES
Evaluate and choose centroid as parent for each itertation 
**Note**: plot and data stored in folder plot_cmp_result1

| Test function       |mml-ES withGP|mml-ES noGP|(1+1)-ES withGP|(1+1)-ES noGP|mml-speed-up|(1+1)-speed-up|mml-(1+1)-GP-ratio| mml-(1+1)-ratio|
| :-------------------|------------:| ---------:|--------------:|------------:|-----------:|-------------:|-----------------:|---------------:|
| linear sphere       |734|2700|505|1275| 3.68|2.52|0.688|0.472|
| quadratic sphere    |301|1400|214|681|  4.65|3.18|0.711|0.486|  
| cubic sphere        |267|960 |202|481|  3.60|2.38|0.756|0.501|
| Schwefel’s function|2599|5490|1490|2381|2.11|1.60|0.573|0.434|
| quartic function   |1024|6250|1259|4218|6.10|3.35|1.229|0.675| 


#### Compare funCalls (4/4,15)-ES

Number of runs = 500
Length scale theta = 8 *sigma*sqrt(n)
(4/4,15)-ES
Evaluate and choose centroid as parent for each itertation
**Note**: plot and data stored in folder plot_cmp_result3

| Test function       |mml-ES withGP|mml-ES noGP|(1+1)-ES withGP|(1+1)-ES noGP|mml-speed-up|(1+1)-speed-up|mml-(1+1)-GP-ratio| mml-(1+1)-ratio|
| :-------------------|------------:| ---------:|--------------:|------------:|-----------:|-------------:|-----------------:|---------------:|
| linear sphere      |670.5|2310|504|1270|      
| quadratic sphere   |238|1180|213|680.5|   
| cubic sphere       |276|800|202|476|
| Schwefel’s function|2683|13825|1246|4233|  
| quartic function   |841|4610|1250|4191|  

#### Next

- [x] Format the code s.t. the output is the same for 4 different strategies. 
	1. t:                  # of objective function calls                    
	2. x_last:             last parent x(centroid)
	3. fx:                 last objective function value for parent(centroid)
	4. sigma_array:        simage arrary over # of objective function calls  
	5. x_array:            parent set for x
	6. fx_array:           objective function values for parents
	7. convergence_rate:   rate of convergence
	8. GP_error:           if no GP error = -1
	9. sigma_star_array:   normalized step size

- [x] Modify the plot fucntion (fun_graph_funCall_merged.m,fun_graph_iterate_merged.m) script s.t. 
	- the data and plot is automatically saved.
	- plot for step size (sigma), objective function value (fx) and normalized step size (sigma*) over # objective function calls and # iterations respectively.
	- plot over # objective function call: run_graph_funCall_merged.m 
	- plot over # iterations call: run_graph_iterate_merged.m

- [x] Modify the training set s.t. the training data is added by a first come first out basis.
	See mml_GP.m

- [x] Double the GP length scale to see what happens.
	- Use the (2/2,5)-ES to test and the result for length scale mutiplied by 8 and 16 are stored in folder plot_cmp_result1 and plot_cmp_result2 respectively.
	- Result: the length scale does not matter in the context.

- [x] Run Arash's (1+1)-ES on the three sphere functions and plot the normalized step size against the iteration number.
	- Plot attached as the subplot combined with f(x) and sigma over number of objective function calls.
	- We need large normalized step size to beat the (1+1)-ES with GP essentially because of the step-size adaptation proposed. It could produce a large normalized step size s.t. the mml-ES with GP cannot beat that.
	- The normalized step size (sigma*) increases as we increase lambda. 
		- sigma* for (2/2,5)-ES < (1+1)-ES
		- sigma* for (3/3,10)-ES > (1+1)-ES
		- The difference grows even larger for (4/4,15)-ES (much larger normalized step size)


- [x] Plot the relative error of GP using the following and plot the GP error of (1+1)-ES with GP and mml-ES with GP accordingly.
	- GP_error = |f(y)-fep(y)|/(f(y)-f(x)| where x is the parent and y offspring.
	- Function for merged plot of the 5 test functions.
	- Below is by taking the median relative error and the number of objective function calls for each iteration by running 50 replicates.

	|strategy\median relative error| linear sphere | quadratic sphere | cubic  sphere | Schwefel’s function | quartic function
 	| :-------- | --------:| ------: |--------:| ------: |------: |
	| (4/4,15) withGP| 1.302|0.9086|2.0487|1.7776|1.2200|
	|(1+1) withGP|1.2042|1.2085|1.9871|1.3764|1.1586|

	|strategy\median funCall| linear sphere | quadratic sphere | cubic  sphere | Schwefel’s function | quartic function
	| :-------- | --------:| ------: |--------:| ------: |------: |
	|(4/4,15) withGP | 633|203|228.5|3005|817.5|
	|(1+1) withGP|488|216|206|1536|1236|

- [ ] Can you extend the GP code to compute variances and plot those in relation to errors and step size?
	GP estimate: fTest = mu + Ks'*(K\(fTrain'-mu))
	variances = Ks'*(K\(fTrain'-mu))
	error = |f(y)-fep(y)|/|f(y)-f(x)|
	Find the relation between variance & error, variance & step size.
	**Intepretation**: make the decision on variance rather than evaluate evaluate one per iteration. 

- [ ]In the (\mu/\mu,\;ambda)$-ES, would evaluating the best candidate solution rather than the centroid make a difference?

- [ ]Avoid the increase in objective function value by comparing the objective function value of the current centroid and the last one.
	- **Problem**: may easily trapped in a local minima E.G. a saddle point. 
	
- [ ] Weight using the GP model and the true objective function
	- For large lambda parameter setting: https://arxiv.org/pdf/1604.00772.pdf 
	- Ensure the surrogate model can be trusted: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.448.9371&rep=rep1&type=pdf
	- Balance the surroagte model and the true objective function: https://arxiv.org/pdf/1204.2356.pdf

### 8. Alternative strategy 

#### a. Original 
	
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
- Result for averaging 5 runs

	|mu     | lambda|  linear   |  quadratic|    cubic  |   Schwefel |  quartic  |
	| :---  | -----:| --------: |----------:| ---------:|----------: |--------:  |
	|2.00 	| 10.00 |	 697.00 |	 308.00 |	 267.00 |	 2848.00 |	 1084.00 |	
	|2.00 	| 20.00 |	 608.00 |	 219.00 |	 285.00 |	 3040.00 |	 748.00	 |
	|2.00 	| 30.00 |	 569.00 |	 208.00 |	 281.00 |	 3040.00 |	 762.00  |
	|2.00 	| 40.00 |	 580.00 |	 198.00 |	 364.00 |	 3040.00 |	 846.00  |
	|2.00 	| 50.00 |	 587.00 |	 195.00 |	 318.00 |	 3040.00 |	 3040.00 |
	|3.00 	| 10.00 |	 743.00 |	 272.00 |	 296.00 |	 2844.00 |	 969.00  | 
	|3.00 	| 20.00 |	 606.00 |	 216.00 |	 250.00 |	 3040.00 |	 715.00  |
	|3.00 	| 30.00 |	 606.00 |	 198.00 |	 269.00 |	 3040.00 |	 623.00  |
	|3.00 	| 40.00 |	 592.00 |	 197.00 |	 285.00 |	 3040.00 |	 3040.00 |
	|3.00 	| 50.00 |	 607.00 |	 195.00 |	 3040.00|	 3040.00 |	 3040.00 | 
	|4.00 	| 10.00 |	 785.00 |	 308.00 |	 270.00 |	 2654.00 |	 1100.00 |
	|4.00 	| 20.00 |	 641.00 |	 211.00 |	 264.00 |	 3040.00 |	 674.00  | 
	|4.00 	| 30.00 |	 609.00 |	 196.00 |	 269.00 |	 3040.00 |	 686.00  |
	|4.00 	| 40.00 |	 650.00 |	 191.00 |	 291.00 |	 3040.00 |	 3040.00 |
	|4.00 	| 50.00 |	 677.00 |	 190.00 |	 408.00 |	 3040.00 |	 3040.00 |
	|5.00 	| 10.00 |	 860.00 |	 329.00 |	 327.00 |	 2429.00 |	 1352.00 | 
	|5.00 	| 20.00 |	 654.00 |	 207.00 |	 274.00 |	 3040.00 |	 736.00  | 
	|5.00 	| 30.00 |	 656.00 |	 204.00 |	 273.00 |	 3040.00 |	 607.00  |  
	|5.00 	| 40.00 |	 695.00 |	 186.00 |	 315.00 |	 3040.00 |	 3040.00 | 

#### b. choose the best centrod so far as parent 

```
while not terminate() do{
	for j = 1,...,lambda{
		z(j) = randn(N,1)
		y(j) = x + sigma*z(j)
		fep(j) = evaluateUseGP(y(i))
	}
	z = sorted(z)
	z = 1/mu*sum(z(1),...,z(mu)))
	x_temp = x + sigma*z
	// update iff. centroid is better
	if(f(x) > f(x_temp)){
		x = x_temp
	}
}
	update step size sigma
}
```
- Result for averaging 5 runs

	|mu     | lambda|  linear   |  quadratic|    cubic  |   Schwefel |  quartic  |
	| :---  | -----:| --------: |----------:| ---------:|----------: |--------:  |
	|2.00   | 8.00 	|	791.00 	| 	 342.00 |	 399.00 |	 1794.00 |	 1147.00 |
	|2.00   | 9.00  |	722.00 	| 	 342.00 |	 340.00 |	 1420.00 |	 1022.00 |
	|2.00   | 10.00 |	 756.00 |	 306.00 |	 383.00 |	 1286.00 |	 1132.00 |
	|2.00   | 11.00 |	 696.00 |	 324.00 |	 463.00 |	 5040.00 |	 1048.00 |
	|2.00   | 12.00 |	 703.00 |	 310.00 |	 816.00 |	 1405.00 |	 960.00  |
	|2.00   | 13.00 |	 721.00 |	 310.00 |	 828.00 |	 5040.00 |	 967.00  |
	|2.00   | 14.00 | 	 719.00 |	 271.00 |	 809.00 |	 5040.00 |	 952.00	 |
	|2.00 	| 15.00 |	 771.00 |	 283.00 |	 913.00 |	 1461.00 |	 915.00  |
	|3.00 	| 8.00 	|    784.00 | 	 370.00 |	 378.00 |	 1507.00 |	 1279.00 | 
	|3.00 	| 9.00 	| 	 742.00 | 	 364.00 |	 450.00 |	 1299.00 |	 995.00  |
	|3.00 	| 10.00 |	 734.00 |	 329.00 |	 375.00 |	 1327.00 |	 1099.00 |
	|3.00 	| 11.00 |	 768.00 |	 326.00 |	 531.00 |	 1319.00 |	 1118.00 |
	|3.00 	| 12.00 |	 714.00 |	 315.00 |	 865.00 |	 1571.00 |	 1037.00 |
	|3.00 	| 13.00 |	 767.00 |	 304.00 |	 646.00 |	 1512.00 |	 995.00  |
	|3.00 	| 14.00 |	 745.00 |	 279.00 |	 895.00 |	 1431.00 |	 956.00  |
	|3.00 	| 15.00 |	 731.00 |	 279.00 |	 625.00 |	 2124.00 |	 957.00  |
	|4.00 	| 8.00 	|    903.00 | 	 400.00 |	 426.00 |	 1675.00 |	 1480.00 |
	|4.00 	| 9.00 	|    858.00 | 	 336.00 |	 401.00 |	 1611.00 |	 1242.00 |
	|4.00 	| 10.00 |	 814.00 |	 338.00 |	 420.00 |	 1402.00 |	 1119.00 |
	|4.00 	| 11.00 |	 737.00 |	 368.00 |	 619.00 |	 1349.00 |	 998.00  |
	|4.00 	| 12.00 |	 796.00 |	 308.00 |	 674.00 |	 1345.00 |	 1092.00 |
	|4.00 	| 13.00 |	 763.00 |	 310.00 |	 831.00 |	 1263.00 |	 1098.00 |
	|4.00 	| 14.00 |	 805.00 |	 312.00 |	 945.00 |	 1315.00 |	 1028.00 |
	|4.00 	| 15.00 |	 724.00 |	 300.00 |	 810.00 |	 1316.00 |	 1000.00 |


---

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






