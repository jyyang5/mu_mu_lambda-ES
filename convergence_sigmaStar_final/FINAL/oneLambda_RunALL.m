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
function oneLambda_RunALL(subplotNum,NUM_OF_RUNS)

f = @(x) (x'*x);

% subplotNum = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(subplotNum == 0) % (1+1)-ES
    lambda= 0;
    mu = 0;
    opt_plot_colour = 'c';
elseif (subplotNum == 1)
    lambda = 10;
    mu = ceil(lambda/4);
    opt_plot_colour = 'k';
elseif(subplotNum == 2)
    lambda = 20;
    mu = ceil(lambda/4);
    opt_plot_colour = 'b';
elseif(subplotNum == 3)
    lambda = 40;
    mu = ceil(lambda/4);
    opt_plot_colour = 'm';
else
end


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
elseif(mu ==0 && lambda == 0)        % for (1+1)-ES
    c_mu_lambda = 0;                               
else 
    c_mu_lambda = (lambda-mu)/(2*pi)*nchoosek(lambda,mu)*sum(exp(-x.^2).*(normcdf(x)).^(lambda-mu-1).*(1-normcdf(x)).^(mu-1))*step;
end

% If NOT (1+1)-ES
if(subplotNum~=0)
    figure(4);
    subplot(1,3,subplotNum)
    legend('-DynamicLegend'); 
    % green
    ita = 4;
    %fun_fitness_sigmaStar(4,100000,'g','.',1);        
    fun_precise_fitness_sigmaStar_multi(f,NUM_OF_RUNS,mu,lambda,ita,100,'g','o',0,c_mu_lambda);
    fun_precise_fitness_sigmaStar_multi(f,NUM_OF_RUNS,mu,lambda,ita,10,'g','x',1,c_mu_lambda);        
    disp('green done');


    % red
    ita = 1;
    fun_precise_fitness_sigmaStar_multi(f,NUM_OF_RUNS,mu,lambda,ita,100,'r','o',0,c_mu_lambda);
    fun_precise_fitness_sigmaStar_multi(f,NUM_OF_RUNS,mu,lambda,ita,10,'r','x',1,c_mu_lambda);        
    disp('red done');
    % blue
    ita = 0.25;
    fun_precise_fitness_sigmaStar_multi(f,NUM_OF_RUNS,mu,lambda,ita,100,'b','o',0,c_mu_lambda);
    fun_precise_fitness_sigmaStar_multi(f,NUM_OF_RUNS,mu,lambda,ita,10,'b','x',1,c_mu_lambda);        
    disp('blue done');


    % black
    ita = 0;
    fun_precise_fitness_sigmaStar_multi(f,NUM_OF_RUNS,mu,lambda,ita,100,'k','o',0,c_mu_lambda);
    fun_precise_fitness_sigmaStar_multi(f,NUM_OF_RUNS,mu,lambda,ita,10,'k','x',1,c_mu_lambda);       
    disp('black done');
    
    box on;
%     % normal mml-ES
%     fun_precise_fitness_sigmaStar_multi(f,NUM_OF_RUNS,mu,lambda,ita,10,'k','.',0,c_mu_lambda);
%     disp('dotted black done')
    % save fig. expected fitness gain over normalized step size
    if(subplotNum==3)
        box on;
        p2 = sprintf('expectedFitGain.fig');
        saveas(gcf,p2);
    end
end

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot opt. expected fitness gain over noise-to-signal ratio 
% % experimental (x o for n = 10,100) and solid n->infty
% 
% s_start = 0.0001;
% increment =1;
% s_end = 10+s_start;
% v_array = s_start:increment:s_end;                                          % For experimental result
% 
% % v_expedted_curve_array = 0.0001:0.0001:20;                                % plot curve n -> infty
% v_expedted_curve_array = exp(-2.302585092994046: 0.0461:2.302585092994046+0.01);
% 
% FIG_NUM = 5;
% figure(FIG_NUM);
% subplot(1,2,1);hold on;
% if(subplotNum~=0)
%     fun_precise_optFitGain_over_v(f,NUM_OF_RUNS,mu,lambda,v_array,10,opt_plot_colour,'x',c_mu_lambda,v_expedted_curve_array);    % plot n = 10 (x)
% else 
%     subplotNum_step = 2;
%     fun_precise_optFitGain_over_v_ONE(f,NUM_OF_RUNS,FIG_NUM,subplotNum_step,v_array,10,opt_plot_colour,'x',c_mu_lambda,v_expedted_curve_array);
% end
% disp('* done'); 
% figure(FIG_NUM);
% subplot(1,2,1);hold on;
% fun_precise_optFitGain_over_v(f,NUM_OF_RUNS,mu,lambda,v_array,100,opt_plot_colour,'o',c_mu_lambda,v_expedted_curve_array);   % plot n = 100 (o)
% disp('s done');
% xlim([0 10]);
% ylim([0 inf]);
% box on;
% 
% 
% % % save opt. expected fitGain
% % if(subplotNum==3)
% %     box on;
% %     p2 = sprintf('opt_fitGain.fig');
% %     saveas(gcf,p2);
% % end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % noise-to-signal ratio
% % % v = 0.1:0.00001:10;
% % v = 0.001*2.^(1:1:18);
% % % solve the plot by taking derivative
% % % opt. step size
% % figure(5);
% % plot(v,c_mu_lambda*mu./(sqrt(1+v.*v)));
% % xlabel('noise-to-signal ratio \upsilon','FontSize',15);%
% % ylabel('opt. normalized step size \sigma^*','FontSize',15); 
% % set(gca, 'XScale', 'log');
% % set(gca,'FontSize',15);
% % p1 = sprintf('opt. normalized step size (%d/%d,%d)-ES',mu,mu,lambda);
% % title(p1,'fontsize',20);
% % p2 = sprintf('opt_step_size_%d_%d_%d_ES.fig',mu,mu,lambda);
% % saveas(gcf,p2);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % opt. fitness gain given opt. step size n= 10,100,infty
% % sigma_star = 0.1:0.001:8;
% % v_range = 0.001/2*2.^(1:1:19);
% % v_range_trans = transpose(v_range);
% % % a matrix 
% % row: different noise-to-signal ratio v
% % col: different normalized step size sigmaStar
% % n=10;
% % expected_cure_10 = c_mu_lambda*sigma_star.*(1+sigma_star.^2/2/mu/n)./(sqrt(1+sigma_star.^2/mu/n).*sqrt(1+v_range_trans.^2+sigma_star.^2/2/n))-n*(sqrt(1+sigma_star.^2/mu/n)-1);
% % [max_fitness_array_10 s_index] = max(expected_cure_10,[],2);
% % n=100;
% % expected_cure_100 = c_mu_lambda*sigma_star.*(1+sigma_star.^2/2/mu/n)./(sqrt(1+sigma_star.^2/mu/n).*sqrt(1+v_range_trans.^2+sigma_star.^2/2/n))-n*(sqrt(1+sigma_star.^2/mu/n)-1);
% % [max_fitness_array_100 s_index] = max(expected_cure_100,[],2);
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot opt. expected step size gain over noise-to-signal ratio 
% % solid n->infty
% if(subplotNum~=0)
%     figure(FIG_NUM);hold on;
%     subplot(1,2,2);hold on;
%     v_range = v_expedted_curve_array;
%     sigma_star = c_mu_lambda*mu./(sqrt(1+v_range.*v_range));
%     plot(v_range,sigma_star,'Color',opt_plot_colour);
%     ylabel('opt. normalized step size \sigma^*_{opt}','FontSize',15);
%     xlabel('noise-to-signal-ratio \upsilon','FontSize',15); 
%     % set(gca, 'XScale', 'log');
%     set(gca,'FontSize',15);
%     box on;
%     xlim([0 10]);
%     ylim([0 inf]);
% end
% 
% % Add legend when last run(largest lambda)
% if(subplotNum==3)
%     leg = legend('(1+1)-ES','(3/3,10)-ES','(5/5,20)-ES','(10/10,40)-ES');
%     leg.FontSize = 10;
%     title(leg,'Model assisted');
%     box on;
%     p2 = sprintf('opt_stepSize_fitGain_FINAL.fig');
%     saveas(gcf,p2);
% end


end