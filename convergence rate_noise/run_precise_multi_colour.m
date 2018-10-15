% final version
% run fun_fitness_sigmaStar to plot convergence rate 
% at noise-to-signal ratio 0,.25,1,4 respectively
% also for data with dim = 10, 100

%fun_fitness_sigmaStar(ita,n,scatterColour,typeDot,LINE_OR_NOT)
%Input
%   mu                    # parents AND # of parents replaced by offspring
%   lambda                generated offspring size
%   ita:                  ita = sigma_ep_star /sigma_star
%   n:                    dim of data
%   scatterColour:        r  Red, g  Green, b   Blue, c   Cyan, m  Magenta, y  Yellow, k   Black,w  White
%   typeDot:              type of dot (+  Plus sign, o  Circle, *  Asterisk, x    Cross)
%                         *: n=10     s: square n=100
%   LINE_OR_NOT:          justScatter = 0, includeLine = 1, dottedLine = 2
%   c_mu_lambda:          parameter for ploting the solide line 

%Return
%   scatter or scatter+curve

%
NUM_OF_RUNS = 2;
mu = 3;
lambda = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For expected convergence rate 
step = 0.0000001;
x = -10:step:10;
% expected fitness gain for n->\infty  under mu and lambda
if(mu==3 && lambda==10)
    c_mu_lambda = 1.065389626877247; % expected convergence rate for (3/3,10)-ES
elseif(mu==5 && lambda==20)
    c_mu_lambda = 1.214478382788638; % expected convergence rate for (5/5,20)-ES
elseif(mu==10 && lambda==40)
    c_mu_lambda = 1.242204493664515; % expected convergence rate for (10/10,40)-ES
else 
    c_mu_lambda = (lambda-mu)/(2*pi)*nchoosek(lambda,mu)*sum(exp(-x.^2).*(normcdf(x)).^(lambda-mu-1).*(1-normcdf(x)).^(mu-1))*step;
end


figure(4);
legend('-DynamicLegend'); 
% green
ita = 4;
%fun_fitness_sigmaStar(4,100000,'g','.',1);        
fun_precise_fitness_sigmaStar_multi(NUM_OF_RUNS,mu,lambda,ita,100,'g','s',0,c_mu_lambda);
fun_precise_fitness_sigmaStar_multi(NUM_OF_RUNS,mu,lambda,ita,10,'g','*',1,c_mu_lambda);        
disp('green done');


% red
ita = 1;
fun_precise_fitness_sigmaStar_multi(NUM_OF_RUNS,mu,lambda,ita,100,'r','s',0,c_mu_lambda);
fun_precise_fitness_sigmaStar_multi(NUM_OF_RUNS,mu,lambda,ita,10,'r','*',1,c_mu_lambda);        


disp('red done');
% blue
ita = 0.25;
fun_precise_fitness_sigmaStar_multi(NUM_OF_RUNS,mu,lambda,ita,100,'b','s',0,c_mu_lambda);
fun_precise_fitness_sigmaStar_multi(NUM_OF_RUNS,mu,lambda,ita,10,'b','*',1,c_mu_lambda);        


disp('blue done');


% black
ita = 0;
fun_precise_fitness_sigmaStar_multi(NUM_OF_RUNS,mu,lambda,ita,100,'k','s',0,c_mu_lambda);
fun_precise_fitness_sigmaStar_multi(NUM_OF_RUNS,mu,lambda,ita,10,'k','*',1,c_mu_lambda);        


disp('black done');

% normal mml-ES
fun_precise_fitness_sigmaStar_multi(NUM_OF_RUNS,mu,lambda,ita,10,'k','.',0,c_mu_lambda);

disp('dotted black done')



% plot opt. step size for (n=10,100,infty)
v_array = [0.001 0.005 0.01 0.05 0.1 0.25 0.4 1 2 4 8 16 32 64];                                % plot n=10,100 dots
v_curve_array = 0.1:0.001:v_array(length(v_array));                       % plot curve n -> infty

FIG_NUM = 5;                
fun_precise_optFitGain_over_v(NUM_OF_RUNS,mu,lambda,v_array,10,'r','*',FIG_NUM,c_mu_lambda,v_curve_array);    % plot n = 10 (x)
disp('* done'); 
fun_precise_optFitGain_over_v(NUM_OF_RUNS,mu,lambda,v_array,100,'b','s',FIG_NUM,c_mu_lambda,v_curve_array);   % plot n = 100 (o)
disp('s done');




p2 = sprintf('exp_fitGain_%d_%d_%d_ES.fig',mu,mu,lambda);
ylim([0 inf]);
saveas(gcf,p2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % noise-to-signal ratio
% % v = 0.1:0.00001:10;
% v = 0.001*2.^(1:1:18);
% % solve the plot by taking derivative
% % opt. step size
% figure(5);
% plot(v,c_mu_lambda*mu./(sqrt(1+v.*v)));
% xlabel('noise-to-signal ratio \upsilon','FontSize',15);%
% ylabel('opt. normalized step size \sigma^*','FontSize',15); 
% set(gca, 'XScale', 'log');
% set(gca,'FontSize',15);
% p1 = sprintf('opt. normalized step size (%d/%d,%d)-ES',mu,mu,lambda);
% title(p1,'fontsize',20);
% p2 = sprintf('opt_step_size_%d_%d_%d_ES.fig',mu,mu,lambda);
% saveas(gcf,p2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opt. fitness gain given opt. step size n= 10,100,infty
sigma_star = 0.1:0.001:8;
v_range = 0.001*2.^(1:1:18);
v_range_trans = transpose(v_range);
% a matrix 
% row: different noise-to-signal ratio v
% col: different normalized step size sigmaStar
n=10;
expected_cure_10 = c_mu_lambda*sigma_star.*(1+sigma_star.^2/2/mu/n)./(sqrt(1+sigma_star.^2/mu/n).*sqrt(1+v_range_trans.^2+sigma_star.^2/2/n))-n*(sqrt(1+sigma_star.^2/mu/n)-1);
[max_fitness_array_10 s_index] = max(expected_cure_10,[],2);
n=100;
expected_cure_100 = c_mu_lambda*sigma_star.*(1+sigma_star.^2/2/mu/n)./(sqrt(1+sigma_star.^2/mu/n).*sqrt(1+v_range_trans.^2+sigma_star.^2/2/n))-n*(sqrt(1+sigma_star.^2/mu/n)-1);
[max_fitness_array_100 s_index] = max(expected_cure_100,[],2);




figure(6);
sigma_star = c_mu_lambda*mu./(sqrt(1+v_range.*v_range));
hold on
plot(v_range,sigma_star.*(c_mu_lambda)./sqrt(1+v_range.^2)-sigma_star.^2./(2*mu),'k');
plot(v_range,transpose(max_fitness_array_10),'-.x','Color','r');
plot(v_range,transpose(max_fitness_array_100),'-.o','Color','b');
hold off;
legend('n \rightarrow \infty','n=10','n=100');

ylabel('opt. expected fitness gain \eta','FontSize',15);
xlabel('noise-to-signal-ratio \upsilon','FontSize',15); 
set(gca, 'XScale', 'log');
set(gca,'FontSize',15);
p1 = sprintf('opt. expected fitness gain (%d/%d,%d)-ES',mu,mu,lambda);
title(p1,'fontsize',20);
p2 = sprintf('opt_fitGain_%d_%d_%d_ES.fig',mu,mu,lambda);
saveas(gcf,p2);
