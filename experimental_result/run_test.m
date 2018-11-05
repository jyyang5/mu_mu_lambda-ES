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
% Ini
fname = 1;
n = 10;
lambda = 20;
mu = ceil(lambda/4);
TRAINING_SIZE = 40;
LENGTH_SCALE = 8;
% x0 for mml-ES
x0 = randn(n,mu);
% x1 for (1+1)-ES
x1 = randn(n,1);


sigma0 = 1;
NUM_OF_RUNS = 10;

NUM_OF_ITERATIONS = 10000;
% mml-ES with GP
a = mml_GP_CSA_Niko(fname,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE);
t = cell2mat(a(1));
centroid = cell2mat(a(2));
f_centroid = cell2mat(a(3));
sigma_array = cell2mat(a(4));
T = cell2mat(a(5));
fcentroid_array = cell2mat(a(6));
convergence_rate = cell2mat(a(7));
GP_error = cell2mat(a(8));
sigma_star_array = cell2mat(a(9));
success_rate = cell2mat(a(10));
delta_array = cell2mat(a(11));

% (1+1)-ES with GP
b = withGP(fname,x1,sigma0,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE);
t1 = cell2mat(b(1));
centroid1 = cell2mat(b(2));
f_centroid1 = cell2mat(b(3));
sigma_array1 = cell2mat(b(4));
T1 = cell2mat(b(5));
fcentroid_array1 = cell2mat(b(6));
convergence_rate1 = cell2mat(b(7));
GP_error1 = cell2mat(b(8));
sigma_star_array1 = cell2mat(b(9));
success_rate1 = cell2mat(b(10));
delta_array1 = cell2mat(b(11));


% disp('GP error');
% disp(GP_error);
% disp(GP_error1);
disp('# of objective functions');
disp(T);
disp(T1);
disp('convergence rate');
disp(convergence_rate);
disp(convergence_rate1);
disp('success rate');
disp(success_rate);
disp(success_rate1);
close all;
figure(10);

% hist objective function calls
subplot(1,4,1);
t_start = ceil(TRAINING_SIZE/lambda);
fx_range1 = fcentroid_array(1:t_start);
fx_range2 = fcentroid_array(t_start+1:t);
t_range1 = 1:lambda:lambda*t_start;
t_range2 = lambda*t_start+1:lambda*t_start+length(fx_range2);
semilogy([t_range1 t_range2], [fx_range1 fx_range2]);hold on;% mml with GP
% plot([lambda lambda+1:T-lambda-1],[fcentroid_array(1) fcentroid_array(2:t)]);hold on;
plot(1:1:T1,fcentroid_array1(1:T1));
ylabel('objective function value f(y)','FontSize',15);%
xlabel('objective function evaluations','FontSize',15); 
set(gca, 'YScale', 'log');
legend({'mml-ES','(1+1)-ES'},'FontSize',10); %
title('objective function value','FontSize',20);
hold off;

% f(x) over objective funtion calls
subplot(1,4,2);
sigma_range1 = sigma_array(1:t_start);
sigma_range2 = sigma_array(t_start+1:t);
plot([t_range1 t_range2], [sigma_range1 sigma_range2]);hold on;
% plot(1:1:t,sigma_star_array(1:t));hold on;
plot(1:1:T1,sigma_array1(1:T1));
ylabel('step size \sigma','FontSize',15);%
xlabel('objective function evaluations','FontSize',15); 
set(gca, 'YScale', 'log');
legend({'mml-ES','(1+1)-ES'},'FontSize',10); %
title('step size \sigma','FontSize',20);
hold off;

% sigmaStar over objective funtion calls
subplot(1,4,3);
sigma_star_range1 = sigma_star_array(1:t_start);
sigma_star_range2 = sigma_star_array(t_start+1:t);
plot([t_range1 t_range2], [sigma_star_range1 sigma_star_range2]);hold on;
% plot(1:1:t,sigma_star_array(1:t));hold on;
plot(1:1:T1,sigma_star_array1(1:T1));
ylabel('normalized step size \sigma*','FontSize',15);%
xlabel('objective function evaluations','FontSize',15); 
set(gca, 'YScale', 'log');
legend({'mml-ES','(1+1)-ES'},'FontSize',10); %
title('normalized step size \sigma*','FontSize',20);
hold off;

% relative model error over objective funtion calls
subplot(1,4,4);
GP_error_range1 = GP_error(1:t_start);
GP_error_range2 = GP_error(t_start+1:t);
plot([t_range1 t_range2], [GP_error_range1 GP_error_range2]);hold on;
% plot(1:t,GP_error(1:t));hold on;
plot(TRAINING_SIZE:T1,GP_error1(TRAINING_SIZE:T1));hold on;
% GP smoothed
smoothed_GP_range1 = smoothdata(GP_error_range1,'gaussian',40);
smoothed_GP_range2 = smoothdata(GP_error_range2,'gaussian',40);

plot([t_range1 t_range2], [smoothed_GP_range1 smoothed_GP_range2],'LineWidth',2);hold on;
plot(TRAINING_SIZE+1:T1,smoothdata(GP_error1(TRAINING_SIZE+1:T1),'gaussian',40),'LineWidth',2);
ylabel('relative model error','FontSize',15);%
xlabel('objective function evaluations','FontSize',15); 
set(gca, 'YScale', 'log');
legend({'mml-ES','(1+1)-ES','smoothed mml','smoothed (1+1)'},'FontSize',10); %
title('relative model error','FontSize',20);
hold off;

figure(12);
histogram(delta_array,'Normalization','probability');hold on;
histogram(delta_array1,'Normalization','probability');hold on;
legend({'mml-ES','(1+1)-ES'},'FontSize',10); %
ylabel('Normalized pdf','FontSize',15);%
xlabel('Normalized delta','FontSize',15); 
set(gca, 'YScale', 'log');


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
        val = val + sum(x(1:i))^2;
    end
end
% function val = f4(x)
%     val = 0;
%     for i = 1:1:length(x)
%         temp = 0;
%         for j = 1:1:i
%             temp = temp + x(j);
%         end
%         val = val + temp^2;
%     end
% end

% quartic function
function val = f5(x)
    beta = 1;
    val = 0;
    for i = 1:1:length(x)-1
        val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
    end
end
