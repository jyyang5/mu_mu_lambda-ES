n=10;
lambda = 40;
mu=ceil(lambda/4);
NUM_OF_ITERATIONS = 2000;
x0 = randn(n,mu);
sigma0 = 1;
LENGTH_SCALE = 16;
TRAINING_SIZE = 40;
SUCCESS_RATE = 0.4;
fname = 3;
SIGMA_STAR = 4;

a = bestSoFar_fourProb_GP_SIGMA_STAR(fname,x0,SIGMA_STAR,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE);
t_array = cell2mat(a(1));
sigma_matrix = cell2mat(a(4));
T_array = cell2mat(a(5));
f_x_matrix = cell2mat(a(6));
success_rate_array = cell2mat(a(10));
sigma_star_matrix = cell2mat(a(9));

four_categories = cell2mat(a(12));
% four_categories_test = cell2mat(a(13));
disp(four_categories);
% disp(four_categories_test);

t_start = ceil(TRAINING_SIZE/lambda);
T_range_1 = (0:1:t_start).*(lambda+1)+1;
T_range_2 = (t_start+2:t_array)+lambda*t_start;
T_range = [T_range_1 T_range_2];

figure(10)
subplot(1,3,1)
plot(T_range,f_x_matrix(1:t_array));
set(gca, 'YScale', 'log');

subplot(1,3,2)
plot(T_range,sigma_matrix(1:t_array));
set(gca, 'YScale', 'log');

subplot(1,3,3)
bar(four_categories/sum(four_categories));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = t_array;
T = T_array;
fcentroid_array = sigma_matrix;
sum(fcentroid_array(t_start:T-1)>fcentroid_array(t_start+1:T))/length(fcentroid_array(t_start:T-1))