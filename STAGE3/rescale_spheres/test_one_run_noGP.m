n=8;
lambda = 10;
mu=ceil(lambda/4);
NUM_OF_ITERATIONS = 20000;
x0 = 0.1*randn(n,1);
sigma0 = 1;
LENGTH_SCALE = 20;
TRAINING_SIZE = 40;
SUCCESS_RATE = 0.4;
SIGMA_STAR = 1;


fname = 8;
para = 0.01;
a = onePlusOne(fname,para,x0,sigma0,NUM_OF_ITERATIONS);

para = 100;
b = onePlusOne(fname,para,x0,sigma0,NUM_OF_ITERATIONS);

t_array = cell2mat(a(1));
sigma_matrix = cell2mat(a(4));
T = cell2mat(a(5));
f_x_matrix = cell2mat(a(6));
success_rate_array = cell2mat(a(10));
sigma_star_matrix = cell2mat(a(9));

t_array1 = cell2mat(b(1));
sigma_matrix1 = cell2mat(b(4));
T1 = cell2mat(b(5));
f_x_matrix1 = cell2mat(b(6));
success_rate_array1 = cell2mat(b(10));
sigma_star_matrix1 = cell2mat(b(9));

% four_categories = cell2mat(a(12));
% % four_categories_test = cell2mat(a(13));
% disp(four_categories);
% % disp(four_categories_test);
% eval_ratio = cell2mat(a(13));
% fprintf('Evaluation rate = %.2f\n',eval_ratio);


T_range = 1:T;
T_range1 = 1:T1;

figure(10)
subplot(1,3,1)
plot(T_range,f_x_matrix(T_range));hold on;
plot(T_range1,f_x_matrix1(T_range1));
legend({'$\beta$=0.1','$\beta$=10'},'interpreter','latex')
xlabel('function calls','FontSize',15);
ylabel('function value','FontSize',15);
set(gca, 'YScale', 'log');

subplot(1,3,2)
plot(T_range,sigma_matrix(T_range));hold on;
plot(T_range1,sigma_matrix1(T_range1));
legend({'$\beta$=0.1','$\beta$=10'},'interpreter','latex')
xlabel('function calls','FontSize',15);
ylabel('step size','FontSize',15);
set(gca, 'YScale', 'log');

subplot(1,3,3)
plot(T_range,sigma_star_matrix(T_range));hold on;
plot(T_range1,sigma_star_matrix1(T_range1));
legend({'$\beta$=0.1','$\beta$=10'},'interpreter','latex')
xlabel('function calls','FontSize',15);
ylabel('normalized step size','FontSize',15);
set(gca, 'YScale', 'log');

% subplot(1,4,4)
% bar(four_categories/sum(four_categories));
% str_cell_SIGMA_STAR = {'TN','FP','FN','TP'};
% set(gca,'xticklabel',str_cell_SIGMA_STAR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = t_array;
T = T_array;
fcentroid_array = sigma_matrix;
% sum(fcentroid_array(t_start:T-1)>fcentroid_array(t_start+1:T))/length(fcentroid_array(t_start:T-1))