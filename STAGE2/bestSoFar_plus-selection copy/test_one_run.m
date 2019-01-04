n=10;
lambda = 40;
mu=ceil(lambda/4);
NUM_OF_ITERATIONS = 2000;
x0 = randn(n,mu);
sigma0 = 1;
LENGTH_SCALE = 16;
TRAINING_SIZE = 40;
SUCCESS_RATE = 0.8;
fname = 2;

a = bestSoFar_GP_change_success_rate(fname,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,SUCCESS_RATE);
t_array = cell2mat(a(1));
sigma_matrix = cell2mat(a(4));
T_array = cell2mat(a(5));
f_x_matrix = cell2mat(a(6));
success_rate_array = cell2mat(a(10));
sigma_star_matrix = cell2mat(a(9));



t_start = ceil(TRAINING_SIZE/lambda);
T_range_1 = (0:1:t_start).*(lambda+1)+1;
T_range_2 = (t_start+2:t_array)+lambda*t_start;
T_range = [T_range_1 T_range_2];

plot(T_range,f_x_matrix(1:t_array));
set(gca, 'YScale', 'log');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = t_array;
T = T_array;
fcentroid_array = sigma_matrix;
sum(fcentroid_array(t_start:T-1)>fcentroid_array(t_start+1:T))/length(fcentroid_array(t_start:T-1))