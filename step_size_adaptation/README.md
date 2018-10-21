## Step size adaptation

### CSA(Nico's method)

- Weighted recombination is used 




- Algorithm (without weighted recombination)
```matlab

```


- Algorithm (Don't work )


```matlab
% offspring generation
for i=1:lambda
  x(t+1,i) = m(t) + sigma(t)*N(n,1);
end
% parent for next iteration
m(t+1) = mean(x(t+1,:),2);                                    % m(t+1) = m(t) + mean(x(t+1,:)-m(t),2); 
% update of evolution path 
s(t+1) = (1-c)*s(t) + sqrt(c*(2-c)*mu)*(m(t+1)-m(t))/sigma(t)    % s(t+1) = (1-c)*s(t) + sqrt(c*(2-c)*mu)*z; where z is the mean of mu randn(n,1)


```

### Observation 

- Configurations
    - quadratic sphere
    - sigma0=1
    - (10/10,40)-ES
- Result    
  - Experiment with CSA (sample_plot_10_10_40_ES.fig)
    - relative error is around 0.5 while sigmaStar \in [4,15] 
    - number of iterations t = 169
    - convergence rate = 0.571
    - number of objective evaluations = 209

  - Compare that with the progress rate plot.
      - The dashed magenta line has a similar setting expect that it uses normalized step sie not CSA. 

      The performance at sigmaStar = 8 has a progress rate around 0.6 
      - The CSA did adapt the step size properly

- ToDoNext
  
  - Narrow the gap between the magenta line and the red dash line (or even larger noise-to-signal ratio)
