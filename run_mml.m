%function val = mml(f,x0,sigma_star,sigma_ep_star,lambda,sigma0,NUM_OF_ITERATIONS)
% initialization
% f:                  objective function value
% x0:                 mu initial point
% mu:                 population size
% sigma0:             initial muttaion strength
% NUM_OF_ITERATIONS:  number of maximum iterations



% OPTIMAL:            global optima
% TARGET_DISTANCE:    target distance to global optima
% example input:      fun = @(x) x' * x
%                     noGP(fun, randn(10,1),1,1000) 
f = @(x) (x'*x);
n = 10;
mu = 2;
lambda = 5;
x0 = randn(n,mu);
sigma_star = 1;
sigma_ep_star = 0;
sigma0 = 1;
NUM_OF_ITERATIONS = 5000;
a = mml(f,x0,sigma_star,sigma_ep_star,lambda,sigma0,NUM_OF_ITERATIONS);
centroid = cell2mat(a(2));
fcentroid = cell2mat(a(3));
sigma_array = cell2mat(a(4));