n=10;
lambda = 10;
mu=ceil(lambda/4);
NUM_OF_ITERATIONS = 10000;
sigma0 = 1;
LENGTH_SCALE = 20;
TRAINING_SIZE = 40;
SUCCESS_RATE = 0.4;
SIGMA_STAR = 1;

fname = 6;
para = 8;
C1 = 1.0;
C2 = 1.0;
C3 = 0.2;

x0 = randn(n,1);
a = withGP(fname,para,x0,sigma0,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE);
b = bestSoFar_arashVariant(fname,para,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,C1,C2,C3);

t = cell2mat(a(1));
sigma_matrix = cell2mat(a(4));
T = cell2mat(a(5));
f_x_matrix = cell2mat(a(6));
success_rate_array = cell2mat(a(10));
sigma_star_matrix = cell2mat(a(9));

t1 = cell2mat(b(1));
sigma_matrix1 = cell2mat(b(4));
T1 = cell2mat(b(5));
f_x_matrix1 = cell2mat(b(6));
success_rate_array1 = cell2mat(b(10));
sigma_star_matrix1 = cell2mat(b(9));

four_categories = cell2mat(a(12));
four_categories1 = cell2mat(b(12));
% disp(four_categories_test);
eval_ratio = cell2mat(a(13));
eval_ratio1 = cell2mat(b(13));
fprintf('Evaluation rate = %.2f,%.2f\n',eval_ratio,eval_ratio1);


T_range = 1:T;
T_range1 = 1:T1;

figure(10)
subplot(1,4,1)
plot(T_range,f_x_matrix(T_range));hold on;
plot(T_range1,f_x_matrix1(T_range1));
legend({'(1+1)-ES',sprintf('(%d/%d,%d)-ES',mu,mu,lambda)});
xlabel('function calls','FontSize',15);
ylabel('function value','FontSize',15);
set(gca, 'YScale', 'log');

subplot(1,4,2)
plot(T_range,sigma_matrix(T_range));hold on;
plot(T_range1,sigma_matrix1(T_range1));
legend({'(1+1)-ES',sprintf('(%d/%d,%d)-ES',mu,mu,lambda)});
xlabel('function calls','FontSize',15);
ylabel('step size','FontSize',15);
set(gca, 'YScale', 'log');

subplot(1,4,3)
plot(T_range,sigma_star_matrix(T_range));hold on;
plot(T_range1,sigma_star_matrix1(T_range1));
legend({'(1+1)-ES',sprintf('(%d/%d,%d)-ES',mu,mu,lambda)});
xlabel('function calls','FontSize',15);
ylabel('normalized step size','FontSize',15);
set(gca, 'YScale', 'log');

subplot(1,4,4)
% bar(four_categories/sum(four_categories));hold on;
% bar(four_categories1/sum(four_categories1));
bar([four_categories/sum(four_categories);four_categories1/sum(four_categories1)])
str_cell_SIGMA_STAR = {'(1+1)-ES',sprintf('(%d/%d,%d)-ES',mu,mu,lambda)};
legend({'TN','FP','FN','TP'});
set(gca,'xticklabel',str_cell_SIGMA_STAR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

