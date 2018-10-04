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
f = f1;
n = 10;

ADD_TRAIN_POINTS= 2;
lambda = 20;
mu = floor(lambda/4);
% x0 for mml-ES
%x0 = randn(n,mu);
% x1 for (1+1)-ES
x1 = randn(n,1);


sigma_satr = 1;
NUM_OF_RUNS = 10;

NUM_OF_ITERATIONS = 5000;
% mml-ES with GP
a = mml_sigmaStarGP_centroid_addTrain(f,x0,sigma_satr,lambda,NUM_OF_ITERATIONS,ADD_TRAIN_POINTS);
t = cell2mat(a(1));
centroid = cell2mat(a(2));
f_centroid = cell2mat(a(3));
sigma_array = cell2mat(a(4));
T = cell2mat(a(5));
fcentroid_array = cell2mat(a(6));
convergence_rate = cell2mat(a(7));
GP_error_array = cell2mat(a(12));
GP_error_med = cell2mat(a(8));
sigma_star_array = cell2mat(a(9));
success_rate = cell2mat(a(11));
fep_centroid_array = cell2mat(a(10));
% % (1+1)-ES with GP
% b = withGP(f1,x1,sigma0,NUM_OF_ITERATIONS);
% t1 = cell2mat(b(1));
% centroid1 = cell2mat(b(2));
% f_centroid1 = cell2mat(b(3));
% sigma_array1 = cell2mat(b(4));
% fcentroid_array1 = cell2mat(b(6));
% convergence_rate1 = cell2mat(b(7));
% GP_error1 = cell2mat(b(8));
% sigma_star_array1 = cell2mat(b(9));


disp('convergence_rate');
disp(success_rate);
disp('success_rate');
disp(convergence_rate)
% disp(GP_error);
% disp(GP_error1);
disp('# of iteratons');
disp(t);
disp('# of objective functions');
disp(T);
disp('# GP error');
disp(GP_error_med);

% 
% t = t -1;
% T = T-1;
% figure(11);
% hold on;
% t_range1 = 1:lambda+1:(lambda+1)*TRAINING_FACTOR+1;
% t_range2 = (lambda+1)*TRAINING_FACTOR+2:T;
% f_x_range1 = fcentroid_array(1:TRAINING_FACTOR+1);
% f_x_range2 = fcentroid_array(TRAINING_FACTOR+2:t);
% plot([t_range1 t_range2], [f_x_range1 f_x_range2]);
% fep_x_med_range2 = fep_centroid_array(TRAINING_FACTOR+2:t);
% plot( t_range2, fep_x_med_range2);
% ylabel('logarithmic value','FontSize',15);%
% xlabel('number of objective function calls','FontSize',15); 
% legend({'f(x)','fep(x)'},'FontSize',10); %
% set(gca, 'YScale', 'log')
% hold off;

% [add ADD_TRAIN_POINTS points to training set] convergence plot over number of function evaluations 
temp = 40/lambda;
figure(11);
hold on;
t_range1 = 1:lambda+1:(lambda+1)*temp+1;
t_range2 = (lambda+1)*temp+1+(1+ADD_TRAIN_POINTS):(1+ADD_TRAIN_POINTS):(lambda+1)*temp+1+(1+ADD_TRAIN_POINTS)*(t-temp-1);
f_x_range1 = fcentroid_array(1:temp+1);
f_x_range2 = fcentroid_array(temp+2:t);
plot([t_range1 t_range2], [f_x_range1 f_x_range2]);
fep_x_med_range2 = fep_centroid_array(temp+2:t);
plot( t_range2, fep_x_med_range2);
ylabel('logarithmic value','FontSize',15);%
xlabel('number of objective function calls','FontSize',15); 
legend({'f(x)','fep(x)'},'FontSize',10); %
set(gca, 'YScale', 'log')
hold off;

% % convergence plot over number of function evaluations
% hold on;
% t_range1 = 1:lambda+1:(lambda+1)*TRAINING_FACTOR+1;
% t_range2 = (lambda+1)*TRAINING_FACTOR+2:t+lambda*TRAINING_FACTOR;
% f_x_range1 = fcentroid_array(1:TRAINING_FACTOR+1);
% f_x_range2 = fcentroid_array(TRAINING_FACTOR+2:t);
% plot([t_range1 t_range2], [f_x_range1 f_x_range2]);
% fep_x_med_range2 = fep_centroid_array(TRAINING_FACTOR+2:t);
% plot( t_range2, fep_x_med_range2);
% ylabel('logarithmic value','FontSize',15);%
% xlabel('number of objective function calls','FontSize',15); 
% legend({'f(x)','fep(x)'},'FontSize',10); %
% set(gca, 'YScale', 'log')
% hold off;

% disp(t1);
% 
% figure(10);
% hold on;
% % plot(x_axis,fcentroid_array(1:t));
% plot(GP_error);
% plot(GP_error1);
% ylabel('relative GP error','FontSize',15);%
% xlabel('number of objective function calls','FontSize',15); 
% set(gca, 'YScale', 'log');
% legend({'mml-ES','(1+1)-ES'},'FontSize',10); %
% hold off;

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
% x_axis = 41:1:t-1;
% plot(x_axis,fcentroid_array(41:t-1));
% plot(x_axis,fep_centroid_array(41:t-1));
% set(gca, 'YScale', 'log');
% ylabel('logarithmic value','FontSize',15);%
% xlabel('number of objective function calls','FontSize',15); 
% legend({'f(x)','fep(x)'},'FontSize',10); %
% hold off;




% a_range1 = 1:lambda+1:(lambda+1)*(TRAINING_SIZE/lambda)+1;
% a_range2 = (lambda+1)*(TRAINING_SIZE/lambda)+2:T;
% b_range1 = fcentroid_array(1:(TRAINING_SIZE/lambda)+1);
% b_range2 = fcentroid_array((TRAINING_SIZE/lambda)+2:t);



% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schwefel's Problem 1.2
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
