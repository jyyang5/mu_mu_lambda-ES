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
% x0 for mml-ES
f = @(x) (x'*x);
n = 10;
mu = 4;
lambda = 30;
% x0 for mml-ES
x0 = randn(n,mu);
% x1 for (1+1)-ES
x1 = randn(n,1);


sigma0 = 1;
NUM_OF_RUNS = 10;

NUM_OF_ITERATIONS = 10000;
% mml-ES with GP
a = mml_GP(f1,x0,sigma0,lambda,NUM_OF_ITERATIONS);
t = cell2mat(a(1));
centroid = cell2mat(a(2));
f_centroid = cell2mat(a(3));
sigma_array = cell2mat(a(4));
fcentroid_array = cell2mat(a(6));
convergence_rate = cell2mat(a(7));
GP_error = cell2mat(a(8));
sigma_star_array = cell2mat(a(9));

% (1+1)-ES with GP
b = withGP(f1,x1,sigma0,NUM_OF_ITERATIONS);
t1 = cell2mat(b(1));
centroid1 = cell2mat(b(2));
f_centroid1 = cell2mat(b(3));
sigma_array1 = cell2mat(b(4));
fcentroid_array1 = cell2mat(b(6));
convergence_rate1 = cell2mat(b(7));
GP_error1 = cell2mat(b(8));
sigma_star_array1 = cell2mat(b(9));


% disp('GP error');
% disp(GP_error);
% disp(GP_error1);
disp('# of objective functions');
disp(t);
disp(t1);

figure(10);
hold on;
% plot(x_axis,fcentroid_array(1:t));
plot(GP_error);
plot(GP_error1);
ylabel('relative GP error','FontSize',15);%
xlabel('number of objective function calls','FontSize',15); 
set(gca, 'YScale', 'log');
legend({'mml-ES','(1+1)-ES'},'FontSize',10); %
hold off;

% b = mml(@f5,x0,sigma0,lambda,NUM_OF_ITERATIONS);
% 
% t1 = cell2mat(b(1));
% centroid1 = cell2mat(b(2));
% f_centroid1 = cell2mat(b(3));
% sigma_array1 = cell2mat(b(4));
% fcentroid_array1 = cell2mat(b(6));
% convergence_rate1 = cell2mat(b(7));
% fep_centroid_array1 = cell2mat(b(8));

% disp("number of iterations");
% disp(t);
% % disp(t1);

% disp("convergence rate");
% disp(convergence_rate);
% disp(convergence_rate1);

% T = t;
% figure(10);
% hold on;
% x_axis = 41:1:T-1;
% plot(x_axis,fcentroid_array(41:T-1));
% plot(x_axis,fepcentroid_array(41:T-1));
% set(gca, 'YScale', 'log');
% ylabel('logarithmic value','FontSize',15);%
% xlabel('number of objective function calls','FontSize',15); 
% legend({'f(x)','fep(x)'},'FontSize',10); %
% hold off;
% 
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
    for i = 1:1:length(x)-1
        val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
    end
end