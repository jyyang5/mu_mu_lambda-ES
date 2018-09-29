##### Experiment

- Try different mu and lambda using sigma* for step update. Only evaluate the centroid after the training is done  

mu = lambda/4

sigma_star fixed to some constant

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
	Not gain anything the reuslt is worse than baseline, cannot even solve the Quartic function.

