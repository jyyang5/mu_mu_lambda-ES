## Convergence rate_noise

Model the GP estimate using the true objective function with a Gaussin distributed noise.
Use normalized step size where sigma is proportional to the distance to optimal.

- Folders
    - [noSolid](https://github.com/jyyang5/mu_mu_lambda-ES/tree/master/convergence%20rate_noise/noSolid)

	    Plot 
   	    - Convergence rate for different signal-to-noise ratio and sigmaStar (scatter)
        - v = 0,0.25,1,4
    - [add_solid_line](https://github.com/jyyang5/mu_mu_lambda-ES/tree/master/convergence%20rate_noise/add_solid_line)

        Add
        - Expected value of convergence (curve)
        
    - [solid_merged](hhttps://github.com/jyyang5/mu_mu_lambda-ES/tree/master/convergence%20rate_noise/solid_merged)

        Merge
   	    - Scatter and curve
   	    - Legends formatted
   	    - Into two functions

    

- Files
    
    - For n -> infty 
        - mml_noise.m

            Noise estimation of GP 

        - mml_normal.m

            Use exact function evaluation for offspring ranking

        - fun_fitness_sigmaStar_multi.m

            Run multiple times take the median and plot expected convergence rate over sigmaStar for 
            - noise-to-signal ratio = 0,0.25,1, 4 (different colour)
            - dimension n = 10,100 (different scatter o/x)

        - fun_optFitGain_over_v.m

            Plot the for usual mml-ES 

        - run_multi_colour.m

            Runner for all previous functions

    - For n -> 10,100,infty
        - fun_precise_fitness_sigmaStar_multi.m

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

        - fun_precise_optFitGain_over_v

            Plot opt. progress rate [1opt_fitGain_mu_mu_lambda_ES.fig]
                - n: different colour
                    - n=10: blue ('b')
                    - n=100: red ('r')
                    - n=infty: blacok ('k')

        - run_precise_multi_colour

            Run all the functions above