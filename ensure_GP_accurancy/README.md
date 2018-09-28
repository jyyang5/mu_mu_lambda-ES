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


### 2. Method II

Paper Source:  https://arxiv.org/pdf/1204.2356.pdf

