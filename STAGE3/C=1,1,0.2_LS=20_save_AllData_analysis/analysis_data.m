

% % Count # of diverge 
% f_diverge_COUNT = [sum(sum(sum(T_f6 == 9999,1),2),3),...
%     sum(sum(sum(T_f7 == 9999,1),2),3),...
%     sum(sum(sum(T_f8 == 9999,1),2),3),...
%     sum(sum(sum(T_f9 == 9999,1),2),3)];
% 
% 
% 
% f_GPmml_diverge_COUNT_all = [sum(sum(sum(T_f6(3:5,:,:) == 9999,1),2),3),...
%     sum(sum(sum(T_f7(3:5,:,:) == 9999,1),2),3),...
%     sum(sum(sum(T_f8(3:5,:,:) == 9999,1),2),3),...
% 	sum(sum(sum(T_f9(3:5,:,:) == 9999,1),2),3)]
% 
% f6_GPmml_diverge_COUNT_REP = sum(T_f6(3:5,:,:),3);
% f7_GPmml_diverge_COUNT_REP = sum(T_f7(3:5,:,:),3);
% f8_GPmml_diverge_COUNT_REP = sum(T_f8(3:5,:,:),3);
% 
% 
% 
% f_GPmml_diverge_PEC = f_GPmml_diverge_COUNT_all./...
%     [sum(sum(sum(T_f6(3:5,:,:) ~= 9999,1),2),3),...
%     sum(sum(sum(T_f7(3:5,:,:) ~= 9999,1),2),3),...
%     sum(sum(sum(T_f8(3:5,:,:) ~= 9999,1),2),3),...
%     sum(sum(sum(T_f9(3:5,:,:) ~= 9999,1),2),3)]
% 
% 
% f_onPlusOne_diverge_COUNT = [sum(sum(sum(T_f6(1,:,:) == 9999,1),2),3),...
%     sum(sum(sum(T_f7(1,:,:) == 9999,1),2),3),...
%     sum(sum(sum(T_f8(1,:,:) == 9999,1),2),3),...
%     sum(sum(sum(T_f9(1,:,:) == 9999,1),2),3)]


% strategy index
s = 5;
% function parameter index
f = 21;
% replicate index 
r = 50;

FIG_NUM = f;


f_x_array = squeeze(f_x_f7(s,f,r,:));
sigma_array = squeeze(sigma_f7(s,f,r,:));
sigma_star_array = squeeze(sigma_star_f7(s,f,r,:));
success_rate = squeeze(success_f7(s,f,r));
eval_rate = squeeze(eval_rate_f7(s,f,r));
four_probs = squeeze(four_prob_f7(s,f,r,:));

T_range = 1:9999;

figure(FIG_NUM)

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






