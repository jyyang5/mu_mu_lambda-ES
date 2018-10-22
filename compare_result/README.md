## Compare result of different ES
#### Files
- Strategies
    - (1+1)-ES (noGP.m) **CSA from Dirk's slides**
    - Surrogate-assisted (1+1)-ES (withGP.m)
    - (mu/mu,lambda)-ES (mml_CSA.m)
    - Surrogate-assisted (mu/mu,lambda)-ES (mml_GP_CSA.m) **CSA from Dirk's slides**
- Plots
    - fun_graph_funCall_NEW.m

        Plot the following over objective function evaluations 
        - objective function value
        - step size
        - normalized step size
        - relative model error

        Plot histogram for the following 
        - objective function evaluations
        - convergence rate
        - success rate

        **NOTE:** plot for normalized step size and relative model error fails (very flat no obvious pattern)

    - run_funCall_NEW.m
        - runner for fun_graph_funCall_NEW.m
        - Run 5 test functions
        - Output the median of objective function evaluations, convergence rate, success rate

 
#### Result (Use CSA from Dirk's slides)

The (10/10,40)-ES achieve better performance in a sense the sucess rate is higher than that of (1+1)-ES

It is suprising that the relative model error in quadratic sphere for (10/10,40)-ES is almost 10 times smaller for (10/10,40)-ES.

Another interesting thing is that the performance is closely related to the relative model error

- Linear sphere
    - The relative model error is relatively similar but the mml-ES has larger variance
    - Gives 5% decrease in success rate (0.152, 0.208) for mml-ES and (1+1)-ES respectively
    - Although success rate for mml-ES is greater (0.525 vs. 0.315)
    - Alomst a quarter more objective function evaluations (651 vs. 494)

- Quadratic sphere
    - Relative model error for mml-ES is much smaller almost 1/8 of that of (1+1)-ES
    - Similar convergence rate (0.553 vs. 0.558) for mml-ES and (1+1)-ES
    - Higher success rate for mml-ES (0.481 vs. 0.428)
    - Similar objective function calls 203 vs. 215 for mml-ES and (1+1)-ES

- Cubic sphere
    - Relative model error for mml-ES relative better 1/4 of that of (1+1)-ES
    - Convergence rate is inferior (0.533 vs. 0.636)
    - Better success rate (0.481 vs. 0.334)
    - 10% inferior than (1+1)-ES (218 vs. 197)
    
- Schwefel's Problem 1.2
    - Fails
- Quartic Function 
    - Better performance 
 

**Conclude**: mml-ES is more sensitive to model error (invariant to test functions in terms of approximately 0.5 success rate), for improvement seek tuning for theta for different test function may be a way out.


---




## Previous
Compare the performance of different EAs namely 
    - (1+1)-ES
    - Surrogate-assisted (1+1)-ES
    - (mu/mu,lambda)-ES
    - Surrogate-assisted (mu/mu,lambda)-ES

#### Performance

The performance is evaluated 
1. Plot the following over number of objective function evaluations and iterations 
    - objective function values 
    - step size sigma
    - normalized step size
    
2. Do the similar plots for different test functions
    - Linear sphere 
    - Quadratic sphere 
    - Cubic sphere 
    - Schwefel function (problem 2.1)
    - Quartic function 

#### Files

- EAs

- Function for plots

- Function runner 

#### Folders
   
