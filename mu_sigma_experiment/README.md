## Overview

- with_CSA_GP
Use CSA and GP estimate, experiment with a combination of mu and sigma 

- with_sigmaStar_GP
Use sigma_star fixed and GP estimate, experiment with a combination of sigma_star and sigma with mu = sigma/4

- tuneTrainSize_sigmaStar_GP
Use sigma_star fixed and GP estimate add TRAINING_FACTOR (for GP training set size), experiment with a combination of sigma_star, sigma and TRAINING_FACTOR 


## Trade-off surrogate (GP) model accurancy and objective function call saving

Two approaches are considered to fit an accurate GP model.

### 1. GP model updated if not accurate

Paper Source:  http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.448.9371&rep=rep1&type=pdf

##### 1.1 Description

The GP model is evaluated using the best mu points evaluated using GP model. If the ranking of the three points remains unchanged using the GP model after they are added to GP training set, then we can say the GP model is accurate.

Otherwise, we take few unevaluated (by objective function) but considered good points (by GP model) in the lambda candidate solutions gennerated. Evaluate them using the true objective function, put them in the GP traning set. This iterates until no points are left or the GP model is considered accurate.


##### 1.3 Experiment

Normalized step size is used to model the optimal step size. i.e. we need the ideal performance to decide whether preceed or not. 

sigma* = 1



##### 1.2 Problem

The number of iterations for the GP_update_assisted_mml-ES is very close or even better than the iteration needed for the mml-ES.

But the number of obejctive functions needed to fit an accurate GP model is tremendous. 

E.G. Quadratic Sphere

| Test function       | iteration number|objective functions|
| :-------------------|-------------------:| ---------:|
|  mml-ES withoutGP   |163|6520|     
| mml-ES accurateGP| 52|1887| 


##### 1.3 Solution 


1. Try different mu and lambda using sigma* for step update (**folder with_sigmaStar_GP**)

mu = lambda/4

convergence rate per function evaluations for differnet combination of lambda and sigma_star.

lambda 3 different values (things flaten out)

f(x) linear first iterations skip -> slope = convergence rate 

**Dynamicly show lengend when plot in a loop**

```
legend('-DynamicLegend');
for 
	d = sprintf('%d',x);
	plot(x,y,'DisplayName',d);hold on;
	legend('-DynamicLegend');
    legend('show');
    drawnow;
end
hold off;
```

- Result
	Not gain anything the reuslt is worse than baseline, cannot even solve the Quartic function

- Problem 
	
	1. The number of objective function call and iobjective function value is not corespondent 
	2. Cannot run testFun 5


2. Fix n_init 

Bound the number of iterations in 



### 2. Training size for GP model (folder tuneTrainSize_sigmaStar_GP)

Add a `TRAINING_FACTOR` parameter in mml_sigmaStar_GP_TrainSize.m file to control the training size.

Training set size is set to `TRAINING_SIZE = lambda*TRAINING_FACTOR;` 

The result of using training size can beat (1+1)-ES with GP for test function 1-4.

**(folder tuneTrainSize_sigmaStar_GP)**

1. Fix dynamically add legend when plot a curve
2. Add legend title
3. Fix ploting function value over number of objective function calls using iteration number 
4. Add TRAINING_FACTOR to model the GP training set size
5. Save the median number of objective function calls to a txt file 


### 3. Method II

Paper Source:  https://arxiv.org/pdf/1204.2356.pdf





