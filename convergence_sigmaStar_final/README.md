## Convergence rate merged (with and without GP)

Use the Gaussian distributed noise to model GP and use GP experiment with progress rate
Use normalized step size where sigma is proportional to the distance to optimal.
Expected progress rate and the coresponding sigmaSatr for different n and lambda is obtained from equation (11) https://github.com/jyyang5/mu_mu_lambda-ES/blob/master/convergence_sigmaStar_merged/expected_progressRate.pdf


### Files

- Different ES
    - mml_normal.m
    - mml_noise.m
    - mml_GP.m
        
- Function runs ES and plot 
    - fun_GP_fitness_sigmaStar.m
        1. Plot expected progress rate for different configuration [1fitGain_mu_mu_lambda_ES.fig]
            - vartheta: different colour 
            - n: different type of line
            - n=10: dashed ('--')
            - n=100: dotted (':')
            - n=infty: solid ('-') 
    
        2. Plot experimental result [1fitGain_mu_mu_lambda_ES.fig]
            - vartheta: different colour 
            - n: different type of scatter 
            - n=10: cross ('x')
            - n=100: circle ('o')
        3. Plot opt. step size [1opt_normalized_step_size_mu_mu_lambda_ES.fig]
            - n: different colours 
            - n=10: blue ('b')
            - n=100: red ('r')
    - fun_precise_fitness_sigmaStar_multi.m
        Plot opt. progress rate [1opt_fitGain_mu_mu_lambda_ES.fig]
            - n: different colour
                - n=10: blue ('b')
                - n=100: red ('r')
                - n=infty: blacok ('k')
    -  fun_precise_optFitGain_over_v.m
        
- Function runner       
    - run_precise_multi_colour.m
        1. Plot opt. progress rate [1opt_fitGain_mu_mu_lambda_ES.fig]
            - n: different colour
            - n=10: blue ('b')
            - n=100: red ('r')
            - n=infty: blacok ('k')
        2. Run all the functions above
        


        
