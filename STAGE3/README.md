## STAGE 3
- Round experiments 
	- Test functions
		- Sphere: varying exp 
		- Quartic: varying \alpha 
		- ellipsoids & sewechfel: varying \beta 

	- Dimensions
		- n = 4,8,10,16

- Plot the speed-ups 
	- Base: (1+1)-ES
	- Compare
		- GP-(1+1)-ES
		- GP-(3/3,10)-ES
		- GP-(5/5,20)-ES
		- GP-(10/10,40)-ES

## 1. Files

### A. Strategies 
- (1+1)-ES [onePlusOne.m]
	- Step size adaptation: 1/5-rule
	- Fitness evaluated by: true objective function 
	- Selection: plus-selection 
- GP-(1+1)-ES[withGP.m]
	- Step size adaptation: 1/5-rule + c1,c2,c3
	- Fitness evaluated by: true objective function(train) + GP
	- Selection: plus-selection 
- GP-mml-ES [bestSoFar_arashVariant.m]
	- Step size adaptation: 1/5-rule + c1,c2,c3
	- Fitness evaluated by: true objective function(train) + GP
	- Selection: plus-selection 

### B. Graphing Function 
- Plot speed-ups for all test functions over dimension [fun_multi_over_dim.M]
	- Input:
		- f6_range: exponents in sphere function 
			(x'*x)^(exp/2)
		- f7_range: $\beta$ in quartic function
			```
			val = 0;
			for i = 1:1:length(x)-1
    			val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
			end
			``` 
		- f8_range: $\alpha$ in ellipsoids
			```
			if length(x) == 1
			    val = x'.*alpha.*x;
			elseif length(x) >= 2
			    val = x'*diag([alpha, ones(1,length(x)-1)])*x;
			end
			``` 
	- Return
		- Save fig 
		- Save data [med_dim=X.mat] 
			- T_med_f6[strategyName, parameterUsedInf6]: objective function calls in median run
			- eval_rate_med_f6[strategyName, parameterUsedInf6]: evaluation rate in median run
			- f_x_med_f6[strategyName, parameterUsedInf6],:]: fx array in median run
- Runner [run_multi_over_dim.m]
	- Input:
		- n: dimension of data 
		- range of parameters
	- Obtain:
		- Each call 

### C. Test strategies files
- Test all strategies once [test_one_run_merged.m]
	- Convergence plots, sigma, sigmaStar, four probs.
- Test (1+1)-ES [test_one_run_noGP.m]
- Test strategies with GP [test_one_run_GP.m]

## 2. Folders

### A. Base setting 
- Setting 
	- C1=1, C2=1, C3=0.1, LS=20 
- *range=5_11replicates*
- *rangeLen=10_31replicates*

### [*kappa*]
- Weights do not sum to 1, but $\kappa$
- Weights use settiing from niko's CMA Tutorial  
```
weights = weights = log((lambda+1)/2)-log(1:mu);
weights = kappa*weights/sum(weights);
```
- **Results**
	- Improve sphere functions, a loss in other functions 
	- Not worth doing 

### [*lengthScale*]


### B. Save all run data
- Objective 
	- Save data of all runs 
    - Test if all runs converge 
- Observation
	- Diverge 
		- sphere (1+1)-ES, GP-(1+1)-ES

#### a. [*save_all_data*]
#### b. [*saveAllData_rangeLen=10_31replicates*]
- Setting 
	- C1=1, C2=1, C3=0.2, LS=20, replicates=31
- Variables 
	- Median 
	    - T_med_fX: [strategyIndex, f6ParaIndex]
	    - f_x_med_fX: [strategyIndex, f6ParaIndex, 50000]
	- AllData
		- T_fX: [strategyIndex, f6ParaIndex, runIndex]
		- f_x_fX: [strategyIndex, f6ParaIndex, runIndex, 50000]

## 3. Analysis [*C=1,1,0.2_LS=20_save_AllData_analysis*]

Analyze the results obtained using C1=1.0,C2=1.0,C3=0.2, length scale = 20

### A. Some GP-mml-ES does not converge 
- Not all converge on quartic function 
	- N=12: $\alpha >= 1.3335$, 
	- Ok
		- a several local minima in the quartic functions
		- Each $\lapha$ has plateaus at a same objective function value 


### B. Some GP-mml-ES does not converge [modified]
- Previous version problem with **training size** (all fixed to 40)
	- Make TRAIN_SIZE = 4*N 
- Experiment with dimension N = 2,4,8,16
	



## Schedule

- [x] 20190116
	- Speed-ups for all strategies with GP over
		- test functions 
		1. Sphere: varying exp 
		2. Quartic: varying \alpha 
		3. ellipsoids & sewechfel: varying \beta 
		- Dimension of data
		n = 4,8,10,16

- [x] 20190117 
	- Curious 
		- weighted recombination 
	    	- mml: 1/4, mu=1, close 
	    	- middle in between 
	    	- should be in the middle  

		- Spheres 
			- Whether large values for exp (speed-up < 1 does not like it)
			- Tune length scale 
		- Quartic
			- Length scale may be independent  
			- 
		- Ellipsoids 
			- Dim: N=4 -> goes up 
				N=10 -> goes up again would be nice 
		- Swefel
			-  
	- Parameter setting 
		- Tune a little bit 
		- Sphere: End N = 8
		- s.t. gain does not lose 
	- Dimension 
		- N = 4,8,16
	- Matrix in GP usually a sign -> diverge (maybe with initialization)

- [x] 20190118
	- Kappa 
		- Weighted recombination 
		- All weights not sums to 1 (but to $\kappa$) i.e. weights = weights * $\kappa$
	- Experiments 
		- Seems larger $\kappa$ improves performance on sphere functions 
		- Not significantly different, GP-(1+1)-ES still outperform all when dim n grows 

- [x] 20190119
	- Save all data 
	

- [x] 20190129 
	- Table
		- Fine with more $\beta$s if have space
	- Speed-up Fig.
		- One test function/row (same scale) log||linear??
		- Get rid of schwefel 
	- Single step behaviour Fig.
		- N = 8
		- 4 cols:
		    - Sphere: quartic 
		    - Quartic: $\beta=1$
		    - Ellipsoids: $\beta = 0.1,1$   
		- 4 rows
			- hist of objective function calls
			- convergence plot 
			- 


