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

% Test functions
f1 = @(x) (x'*x)^(1/2);  % linear sphere
f2 = @(x) (x'*x);        % quadratic sphere
f3 = @(x) (x'*x)^(3/2);  % cubic sphere


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(x) (x'*x);
n = 10;
mu = 3;
x0 = randn(n,mu);
%v = 4;
%sigma_ep_star = v*sigma_star;
lambda = 10;
sigma0 = 1;

NUM_OF_ITERATIONS = 10000;
a = mml_GP(f1,x0,sigma0,lambda,NUM_OF_ITERATIONS);
t = cell2mat(a(1));
centroid = cell2mat(a(2));
f_centroid = cell2mat(a(3));
sigma_array = cell2mat(a(4));
fcentroid_array = cell2mat(a(6));
convergence_rate = cell2mat(a(7));
fep_centroid_array = cell2mat(a(8));

b = mml(f1,x0,sigma0,lambda,NUM_OF_ITERATIONS);

t1 = cell2mat(b(1));
centroid1 = cell2mat(b(2));
f_centroid1 = cell2mat(b(3));
sigma_array1 = cell2mat(b(4));
fcentroid_array1 = cell2mat(b(6));
convergence_rate1 = cell2mat(b(7));
fep_centroid_array1 = cell2mat(b(8));

disp("number of iterations");
disp(t);
disp(t1);

disp("convergence rate");
disp(convergence_rate);
disp(convergence_rate1);

T = 300;
figure(10);
hold on;
x_axis = 41:1:T-1;
scatter(x_axis,fcentroid_array(41:T-1));
scatter(x_axis,fep_centroid_array(41:T-1));
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schwefel?s Problem 1.2
function val = f4(x)
    val = 0;
    for i = 1:1:length(x)
        temp = 0;
        for j = 1:1:i
            temp = temp + x(j);
        end
        val = val + temp^2;
    end
end
% quartic function
function val = f5(x)
    beta = 1;
    val = 0;
    for i = 1:1:n-1
        val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
    end
end
