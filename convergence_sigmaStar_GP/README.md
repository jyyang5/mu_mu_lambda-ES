## Overview

### Different Strategies
All below use sigma_star fixed and GP estimate with GP TRAINING_SIZE=40

- mml_sigmaStarGP_centroid
Evalauate the centroid of the mu best offsprings ranked by GP as parent.

- mml_sigmaStarGP_bestOfMu
Evalauate the best offspring in the best mu offsprings estimated by GP as parent.

- mml_sigmaStarGP_bestSoFar
Evalauate the best offspring so far as parent (evaluate the best one offspring generation ranked by GP. Evaluate use f(x) replace parent iff f(bestOfMu) < f(parent)). 

- mml_sigmaStar_GP_centroidBest
Evalauate the best centroid of the best mu offsprings ranked by GP so far as parent (Evaluate centroid use f(x) replace parent iff f(centroid) < f(parent)). 

- mml_sigmaStarGP_centroid_addTrain
Evalauate the centroid of the best mu offsprings ranked by GP as parent, add ADD_TRAIN_POINTS points evaluated by f(x) to GP training set in each iteration. 
**Then convegence should be divided by number of objective function call per iteration**

<<<<<<< HEAD
=======
- mml_sigmaStarGP_centroidQuadratic
Evaluate the centroid optimized by a quadratic model as parent for offspring generation. Choose the best direction vector by averaging the best mu offspring ranked by GP, then pick two points taking the positive and negative of the best direction (evaluated by GP) together with the previous centroid to form a quadratic model picking the lowest point referred to as the optimized centroid. 

### Plots and data for convergence rate over normalized step size 

- centroid_TRAINING_SIZE=40
Use mml_sigmaStarGP_centroid.m

- bestOfmu_TRAIN_SIZE=40
Use mml_sigmaStarGP_bestOfMu.m

- bestSoFar_TRAIN_SIZE=40
Use mml_sigmaStarGP_bestSoFar.m

- centroidBest_TRAIN_SIZE=40
Use mml_sigmaStar_GP_centroidBest.m

- centroidAddTrain=1_TRAIN_SIZE=40
Use mml_sigmaStarGP_centroid_addTrain.m
**NOTE**: two objective function calls (centroid and the best offspring ranked by GP) per iteration after GP model is built. 

- centroidQuadratic_withGP
Use mml_sigmaStarGP_centroidQuadratic.m

- centroidQuadratic_trueObjFun
Use mml_sigmaStarGP_centroidQuadratic.m
**NOTE**: three objective function calls (optimized centroid, points take the positive and negative direction) per iteration after GP model is built. 



---
>>>>>>> master

### Conclusion 

Taking the centroid is the best.


### What to do next

- Build a quadratic model for the centroid. 

- Adaptation of TRAINING_SIZE should scale with dimension n.

- Adapt the training set according to the model accurancy E.G. (f(x) -fep(x))/std(y) 


## Optimize centroid quadratic model (failed)

- Use GP estimate to evaluate the points (positive and negative directions) to fit the quadratic model
	- The convergence rate is relative small compared to just take the centroid as parent for offspring generation.
	- Very sensitive to GP estimation. Cases with a large lambda (lambda>40), the strategy is likely to diverge given a completely wrong lowerest point. Possibility would be:
		1. The GP estimate for the two points sampled using positive and negative directions are completely inaccurate.
		2. Large lambda may give offsprings that is beyond the trust region of GP model. So that the model cannot be trusted anymore.  
	
- Use true objective function to evaluate the points (positive and negative directions) to fit the quadratic model
	- The Gaussian noise is even larger then the one using GP estimates for fitting the quadratic model. 
    - More interestingly, even if the centroid is optimized by fitting a quadratic model with true objective function evaluation, the strategy converges with a large lambda(E.G. lambda > 25).

<<<<<<< HEAD
## Adapt length scale in GP (in progress)
=======
## Kernel and length scale for GP  

- Kernel function
- Length scale 
	- The square exponential works best for the quadratic sphere
	- Try change length scale (two smllaer two larger) for linear and cubic sphere see if it could match the performance of quadratic sphere. 
>>>>>>> master

The length scale factor 8 is the most appropriate for the quadratic sphere where the large convergence rate for quadratic sphere with two relative convergence rate for linear and cubic sphere can be expalained.

We want to find the optimzal length sacle factor for the three sphere functions and think about adapting the length scale factor accordingly for a relative large convergence rate.

Plot convergence rate over a range of theta (2,4,8,16,32) for different sigmaStar and lambda over the three sphere functions.



