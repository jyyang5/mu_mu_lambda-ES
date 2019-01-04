n=10;
mu=3;
lambda = 10;
NUM_OF_ITERATIONS = 2000;
x0 = randn(n,mu);
sigma0 = 1;
LENGTH_SCALE = 16;
TRAINING_SIZE = 40;
SUCCESS_RATE = 0.8;
fname = 2;

a = bestOfTwo_GP_change_success_rate(fname,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,SUCCESS_RATE);
t_array = cell2mat(a(1));
sigma_matrix = cell2mat(a(4));
T_array = cell2mat(a(5));
f_x_matrix = cell2mat(a(6));
success_rate_array = cell2mat(a(10));
sigma_star_matrix = cell2mat(a(9)); 

plot(f_x_matrix);
set(gca, 'YScale', 'log');