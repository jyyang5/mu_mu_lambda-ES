function genData(c1,c2,c3,mu,lambda,LENGTH_SCALE,n)
N=21;    % number of runs
evals = @(x) x(end, 1);

% c1 = 0.8;
% c2 = 0.8;
% c3 = 0.1;
% LENGTH_SCALE = 20;
% lambda = 10;
% 
% mu = ceil(lambda/4);


f1 = @(x) sqrt(x'*x);
f2 = @(x) x'*x;
f3 = @(x) sqrt(x'*x)^3;
f4 = @(x) sum(cumsum(x).^2);                                    % Schwefel's function
f5 = @(x) sum((x(2:end)-x(1:end-1).^2).^2 + (x(1:end-1)-1).^2); % quartic function

rec = cell(N, 5);
for k=1:N
    [~, rec{k, 1}] = ES(f1, mu, lambda, 1.0e+00*randn(n, 1), 1.0, LENGTH_SCALE, 40, c1, c2, c3);
    [~, rec{k, 2}] = ES(f2, mu, lambda, 1.0e+00*randn(n, 1), 1.0, LENGTH_SCALE, 40, c1, c2, c3);
    [~, rec{k, 3}] = ES(f3, mu, lambda, 1.0e+00*randn(n, 1), 1.0, LENGTH_SCALE, 40, c1, c2, c3);
    [~, rec{k, 4}] = ES(f4, mu, lambda, 1.0e+00*randn(n, 1), 1.0, LENGTH_SCALE, 40, c1, c2, c3);
    [~, rec{k, 5}] = ES(f5, mu, lambda, 1.0e+00*randn(n, 1), 1.0, LENGTH_SCALE, 40, c1, c2, c3);
end

Nevals = cellfun(evals, rec);
[~, idx] = sort(Nevals);
for k=1:5
    r(k) = rec(idx((N+1)/2, k), k);
end

fprintf("lambda=%d, \theta=%d",lambda,LENGTH_SCALE);
fEvals = median(Nevals) 

T = [502,212,202,1503,1250];
T./fEvals

end



