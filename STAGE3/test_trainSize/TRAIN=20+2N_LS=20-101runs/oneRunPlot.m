% (5,11,51)
s = 5;
p = 11;
r = 51;
f_x_array = squeeze(f_x_f6(s,p,r,:));
sigma_array = squeeze(sigma_f6(s,p,r,:));
sigma_star_array = squeeze(sigma_star_f6(s,p,r,:));
success_rate = squeeze(success_f6(s,p,r));
eval_rate = squeeze(eval_rate_f6(s,p,r));
four_probs = squeeze(four_prob_f6(s,p,r,:));

T_range = 1:9999;

figure(p)

subplot(1,3,1)
hold on;
plot(T_range,f_x_array(T_range));hold on;
xlabel('function calls','FontSize',15);
ylabel('function value','FontSize',15);
set(gca, 'YScale', 'log');

subplot(1,3,2)
hold on;
plot(T_range,sigma_array(T_range));hold on;
xlabel('function calls','FontSize',15);
ylabel('step size','FontSize',15);
set(gca, 'YScale', 'log');

subplot(1,3,3)
hold on;
plot(T_range,sigma_star_array(T_range));hold on;
xlabel('function calls','FontSize',15);
ylabel('normalized step size','FontSize',15);
set(gca, 'YScale', 'log');