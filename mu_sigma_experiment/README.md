## Overview

- with_CSA_GP
Use CSA and GP estimate, experiment with a combination of mu and sigma 

- with_sigmaStar_GP
Use sigma_star fixed and GP estimate, experiment with a combination of sigma_star and sigma with mu = sigma/4

- tuneTrainSize_sigmaStar_GP
Use sigma_star fixed and GP estimate add TRAINING_FACTOR (for GP training set size), experiment with a combination of sigma_star, sigma and TRAINING_FACTOR 

- tuneTrainSize_CSA_GP
Use CSA and GP estimate add TRAINING_FACTOR (for GP training set size), experiment if the result obtained by tuneTrainSize_sigmaStar_GP holds.

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



### 2. Training size for GP model with sigma_star fixed (folder tuneTrainSize_sigmaStar_GP)

Add a `TRAINING_FACTOR` parameter in mml_sigmaStar_GP_TrainSize.m file to control the training size.

Training set size is set to `TRAINING_SIZE = lambda*TRAINING_FACTOR;` 

The result of using training size can beat (1+1)-ES with GP for test function 1-4.

**(folder tuneTrainSize_sigmaStar_GP)**

1. Fix dynamically add legend when plot a curve
2. Add legend title
3. Fix ploting function value over number of objective function calls using iteration number 
4. Add TRAINING_FACTOR to model the GP training set size
5. Save the median number of objective function calls to a txt file 

- Result 

Obtained by taking median of 10 runs.

| (mu/mu,lambda)-ES |sigma_star|TRAINING_FACTOR F |linear sphere|quadratic sphere|cubic sphere |
| :-----------------|---------:| ----------------:| -----------:|---------------:|------------:|
|Baseline (Arash's) |          |                  |503			|214		     |198		   |
|(3/3,10)  			|1		   |3		          |296		    |164		     |179          |
|(3/3,10)  			|2		   |3		          |348		    |166		     |162          |
|(3/3,10)  			|4		   |3		          |220		    |745		     |866          |
|(3/3,10)  			|5		   |3		          |269		    |infty		     |infty        |
|(3/3,10)  			|1		   |4		          |305		    |180		     |196          |
|(3/3,10)  			|2		   |4		          |356		    |168		     |197          |
|(3/3,10)  			|4		   |4		          |1334		    |339		     |495          |
|(3/3,10)  			|5		   |4		          |2626		    |infty		     |infty        |
|(9/9,30)  			|1		   |3		          |297		    |207		     |205          |
|(9/9,30)  			|2		   |3		          |291		    |174		     |249          |
|(9/9,30)  			|4		   |3		          |390		    |187		     |359          |
|(9/9,30)  			|5		   |3		          |497		    |196		     |435          |
|(9/9,30)  			|1		   |4		          |327		    |235		     |226          |
|(9/9,30)  			|2		   |4		          |321		    |205		     |267          |
|(9/9,30)  			|4		   |4		          |421		    |217		     |508          |
|(9/9,30)  			|5		   |4		          |515		    |217		     |infty        |
|(14/14,50)			|1		   |3		          |338		    |259		     |262          |
|(14/14,50)			|2		   |3		          |319		    |226		     |276          |
|(14/14,50)			|4		   |3		          |412		    |233		     |651          |
|(14/14,50)			|5		   |3		          |464		    |243		     |infty        |
|(14/14,50)			|1		   |4		          |388		    |305		     |304          |
|(14/14,50)			|2		   |4		          |365		    |279		     |338          |
|(14/14,50)			|4		   |4		          |447		    |284		     |infty        |
|(14/14,50)			|5		   |4		          |519		    |290		     |999          |


|(19/19,70)			|1		   |3		          |393		    |314		     |313          |
|(19/19,70)			|2		   |3		          |375		    |290		     |347          |
|(19/19,70)			|4		   |3		          |445		    |285		     |infty        |
|(19/19,70)			|5		   |3		          |502		    |294		     |infty        |

|(19/19,70)			|1		   |4		          |462		    |384		     |381          |
|(19/19,70)			|2		   |4		          |436		    |350		     |404          |
|(19/19,70)			|4		   |4		          |559		    |355		     |infty        |







### 3. Training size for GP model with CSA (folder tuneTrainSize_CSA_GP)

Add a `TRAINING_FACTOR` parameter in mml_CSA_GP_TrainSize.m file to control the training size.
First try TRAINING_FACTOR F=3 with sigma = 30




### 4. Method II

Paper Source:  https://arxiv.org/pdf/1204.2356.pdf





