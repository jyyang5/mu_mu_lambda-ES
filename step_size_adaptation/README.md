## Step size adaptation

### CSA(Nico's method)

- Weighted recombination is used 

- Algorithm (without weighted recombination)
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