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
mu = 3;
%x0 = randn(n,mu);
lambda = 10;
sigma0 = 5;
v = 4;
sigma_ep_star = v*sigma_star;
NUM_OF_ITERATIONS = 5000;
a = mml_csa(f,x0,sigma0,sigma_ep_star,lambda,NUM_OF_ITERATIONS);
t = cell2mat(a(1));
disp("number of iterations");
disp(t);

t = cell2mat(a(1));
centroid = cell2mat(a(2));
f_centroid = cell2mat(a(3));
sigma_array = cell2mat(a(4));
fcentroid_array = cell2mat(a(6));
convergence_rate = cell2mat(a(7));
s_array = cell2mat(a(8));



disp("convergence rate");
disp(convergence_rate);