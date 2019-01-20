n=10;
mu=ceil(lambda/4);
NUM_OF_ITERATIONS = 20000;
sigma0 = 1;
TRAINING_SIZE = 40;
SUCCESS_RATE = 0.4;
SIGMA_STAR = 1;

close all;


fname = 8;
para = 0.6;


x0 = randn(n,1);
a = onePlusOne(fname,para,x0,sigma0,50000);
LENGTH_SCALE = 8;
b = withGP(fname,para,x0,sigma0,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE);

C1 = 1.0;
C2 = 1.0;
C3 = 0.2;
LENGTH_SCALE = 20;
kappa = 2;

lambda = 10;
c1 = bestSoFar_arashVariant(fname,para,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,C1,C2,C3);
d1 = bestSoFar_arashVariant_w(fname,para,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,C1,C2,C3,kappa);
lambda = 20;
c2 = bestSoFar_arashVariant(fname,para,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,C1,C2,C3);
d2 = bestSoFar_arashVariant_w(fname,para,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,C1,C2,C3,kappa);

lambda = 40;
c3 = bestSoFar_arashVariant(fname,para,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,C1,C2,C3);
d3 = bestSoFar_arashVariant_w(fname,para,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,C1,C2,C3,kappa);



t = cell2mat(a(1));
sigma_matrix = cell2mat(a(4));
T = cell2mat(a(5));
f_x_matrix = cell2mat(a(6));
success_rate_array = cell2mat(a(10));
sigma_star_matrix = cell2mat(a(9));


t_b = cell2mat(b(1));
sigma_matrix_b = cell2mat(b(4));
T_b = cell2mat(b(5));
f_x_matrix_b = cell2mat(b(6));
success_rate_array_b = cell2mat(b(10));
sigma_star_matrix_b = cell2mat(b(9));
four_categories_b = cell2mat(b(12));
% disp(four_categories_test);
eval_ratio_b = cell2mat(b(13));

t1 = cell2mat(c1(1));
sigma_matrix1 = cell2mat(c1(4));
T1 = cell2mat(c1(5));
f_x_matrix1 = cell2mat(c1(6));
success_rate_array1 = cell2mat(c1(10));
sigma_star_matrix1 = cell2mat(c1(9));
four_categories1 = cell2mat(c1(12));
eval_ratio1 = cell2mat(c1(13));

t2 = cell2mat(c2(1));
sigma_matrix2 = cell2mat(c2(4));
T2 = cell2mat(c2(5));
f_x_matrix2 = cell2mat(c2(6));
success_rate_array2 = cell2mat(c2(10));
sigma_star_matrix2 = cell2mat(c2(9));
four_categories2 = cell2mat(c2(12));
eval_ratio2 = cell2mat(c2(13));

t3 = cell2mat(c3(1));
sigma_matrix3 = cell2mat(c3(4));
T3 = cell2mat(c3(5));
f_x_matrix3 = cell2mat(c3(6));
success_rate_array3 = cell2mat(c3(10));
sigma_star_matrix3 = cell2mat(c3(9));
four_categories3 = cell2mat(c3(12));
eval_ratio3 = cell2mat(c3(13));

T4 = cell2mat(d1(5));
T5 = cell2mat(d2(5));
T6 = cell2mat(d3(5));
f_x_matrix4 = cell2mat(d1(6));
f_x_matrix5 = cell2mat(d1(6));
f_x_matrix6 = cell2mat(d1(6));




fprintf('Evaluation rate = %.2f,%.2f\n',eval_ratio1,eval_ratio2);

if T == 99999
    T_range = 50000;
else
    T_range = 1:T;
end
if T_b == 99999
    T_range_b = 50000;
else
    T_range_b = 1:T_b;
end
if T1 == 99999
    T_range1 = 50000;
else
    T_range1 = 1:T1;
end
if T2 == 99999
    T_range2 = 50000;
else
    T_range2 = 1:T2;
end
if T3 == 99999
    T_range3 = 50000;
else
    T_range3 = 1:T3;
end


figure(10)
subplot(1,4,1)
plot(T_range,f_x_matrix(T_range));hold on;
plot(T_range_b,f_x_matrix_b(T_range_b));hold on;
plot(T_range1,f_x_matrix1(T_range1));hold on;
plot(T_range2,f_x_matrix2(T_range2));hold on;
plot(T_range3,f_x_matrix3(T_range3));
legend({'(1+1)-ES','GP-(1+1)-ES','GP-(3/3,10)-ES','GP-(5/5,20)-ES','GP-(10/10,40)-ES'});
% legend({'(1+1)-ES','GP-(1+1)-ES',sprintf('GP-(%d/%d,%d)-ES',mu,mu,lambda)});
xlabel('function calls','FontSize',15);
ylabel('function value','FontSize',15);
set(gca, 'YScale', 'log');

subplot(1,4,2)
plot(T_range,sigma_matrix(T_range));hold on;
plot(T_range_b,sigma_matrix_b(T_range_b));hold on;
plot(T_range1,sigma_matrix1(T_range1));hold on;
plot(T_range2,sigma_matrix2(T_range2));hold on;
plot(T_range3,sigma_matrix3(T_range3));

legend({'(1+1)-ES','GP-(1+1)-ES','GP-(3/3,10)-ES','GP-(5/5,20)-ES','GP-(10/10,40)-ES'});
xlabel('function calls','FontSize',15);
ylabel('step size','FontSize',15);
set(gca, 'YScale', 'log');

subplot(1,4,3)
plot(T_range,sigma_star_matrix(T_range));hold on;
plot(T_range_b,sigma_star_matrix_b(T_range_b));hold on;
plot(T_range1,sigma_star_matrix1(T_range1));hold on;
plot(T_range2,sigma_star_matrix2(T_range2));hold on;
plot(T_range3,sigma_star_matrix3(T_range3));

legend({'(1+1)-ES','GP-(1+1)-ES','GP-(3/3,10)-ES','GP-(5/5,20)-ES','GP-(10/10,40)-ES'});
xlabel('function calls','FontSize',15);
ylabel('normalized step size','FontSize',15);
set(gca, 'YScale', 'log');

subplot(1,4,4)
% bar(four_categories/sum(four_categories));hold on;
% bar(four_categories1/sum(four_categories1));
bar([four_categories_b/sum(four_categories_b);four_categories1/sum(four_categories1);...
    four_categories2/sum(four_categories2);four_categories3/sum(four_categories3)]);
str_cell_SIGMA_STAR = {'GP-(1+1)-ES','GP-(3/3,10)-ES','GP-(5/5,20)-ES','GP-(10/10,40)-ES'};
legend({'TN','FP','FN','TP'});
set(gca,'xticklabel',str_cell_SIGMA_STAR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

