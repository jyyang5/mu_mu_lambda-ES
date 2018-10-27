## Dart art 

In this folder we tried to modify step size adpatation by using a range of hybried approaches in terms of CSA or some idea borrowed from 1/5-rule

The most successful one is the one using CSA with a emergency tigger 

**Remaining Question**: why is CSA perform better than the CSA ???
See MATLAB.fig (cmp_delta_fitGain.fig, cmp_fx_sigma_modelError_sigmaStar.fig)
The below result uses lambda = 40, mu = ceil(lambda/4), sigma0=1

| Step size adaptation  |  test function|convergence rate|  success rate| objFunCalls|number of emergencey|
| :---------------|-------:| ------------------:|--------------:|------------:|--------:|
| CSA|               quadratic sphere|0.4907|0.5045 | 221|0
|CSA with emergency| quadratic sphere|0.6496|0.2235 | 180|100/140

- Is a larger convergence that importnat???
- Is there a better way to balance the trade off between convergence rate and success rate 
- Is 0.714 emergency cheating in the context of sphere funtions 
- How is emergency rate affected by the DECREASE_FACTOR  

### Files 
- mml_GP_final_emergency.m

  Final version of the CSA with emergency (described in the MATLAB pesudo code below)
  - Fix problem with **negative** normalized fitness gain/iteration 
  - Get rid of for-loops calculating the array of n  

- mml_GP_CSA_Niko_hybrid.m
```matlab
set DECREASE_FACTOR          % a constnat < 1
y = x + sigma *z;            % x parent y offspring
fy = f(y);  
xTrain(n,T) = y;
fTrain(T) = fy;
% Training havent yet done
if(t<=TRAINING_SIZE)
  update sigma using CSA;
  x = y;
  fx = fy;
  
% Offspring inferior 
elseif(f(y) > f(x))
  sigma=sigma*DECREASE_FACTOR;

% Offspring superior
else
  update sigma using CSA;
  x = y
end


```

- mml_GP_CSA_Niko.m

Just use CSA from Niko's CMA tutorial paper.
```

```

Just use Niko's CSA from CMA tutorial paper 2016

### Folders

- attempt1
    - Use mml_GP_CSA_v1.m
    - Add a l-2 norm in Niko's CSA in CMA tutorial slides
    - else similar to mml_GP_CSA_Niko_hybrid.m

- niko_CSA_attempt_decreaseFactor
    - Use mml_GP_CSA_Niko_hybrid.m

- niko_CSA_pure
    - Use mml_GP_CSA_Niko.m

