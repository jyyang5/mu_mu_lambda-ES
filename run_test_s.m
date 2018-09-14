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
lambda = 10;
sigma0 = 1;
sigma_star = 1;
sigma_ep_star = 4*sigma_star;
NUM_OF_ITERATIONS = 5000;
a = mml_noise(f,x0,sigma0,sigma_star,sigma_ep_star,lambda,NUM_OF_ITERATIONS);
t = cell2mat(a(1));
disp("number of iterations");
disp(t);

centroid = cell2mat(a(5));
fcentroid = cell2mat(a(6));
sigma_array = cell2mat(a(4));
convergence_rate = cell2mat(a(7));
t_gp = cell2mat(a(8));
t = cell2mat(a(1));
s_array = cell2mat(a(9));



disp("convergence rate");
disp(convergence_rate);