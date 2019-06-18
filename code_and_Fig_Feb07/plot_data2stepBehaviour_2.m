% test file
% plot from data saved
% INPUTEL: n:  dim of saved data
% NOTE: DO NOT CHANGE THE NUMBER WHEN CALL fun_graph_merged
%       OR CANNOT SAVE ITERTAION DATA
% MODIFICATION 
%       take median run values (mediean objective function calls)
%       Plot 1
%           Row 1: histogram of objective function calls
%           Row 2: convergence plots
%           Row 3: step size
%       Plot 2
%           Row 1: normalized convergence rate
%           Row 2: normalized success rate [good step]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data
% load('all_dim=8.mat');

% Load data with dim n;
n = 8;
load(sprintf('data_dim=%d.mat',n));

% For compact subplots 
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.04], [0.06 0.04], [0.06 0.04]);
if ~make_it_tight,  clear subplot;  end


% GP smooth 
window_length = 40;
kernel = exp(-(-3*window_length:3*window_length).^2/window_length^2/2);
kernel = kernel/sum(kernel);        % Normalized   


para_array = [1 ,1, 1, 2];

subplot_ROW = 4;
subplot_COL = 4; 
STRATEGY_NUM = 5;
FIG_NUM = 100;
% n = 16;

for fname = 1:1:subplot_COL
    para = para_array(fname);
    if fname == 1
            f_x_med =  f_x_med_f6;
            sigma_med  =  sigma_med_f6;
            T_med_matrix = T_med_f6;
            t_med_matrix = t_med_f6;
            eval_rate_med = eval_rate_med_f6;
            error_array_med = error_array_med_f6;
            
            T_matrix =  T_f6;
            four_prob = four_prob_med_f6;
            eval_rate = eval_rate_med_f6;
            f_range = f6_range;
    elseif fname == 2
            f_x_med =  f_x_med_f7;
            sigma_med  =    sigma_med_f7;
            T_med_matrix = T_med_f7;
            t_med_matrix = t_med_f7;
            eval_rate_med = eval_rate_med_f7;
            error_array_med = error_array_med_f7;
            
            T_matrix =  T_f7;
            four_prob = four_prob_med_f7;
            eval_rate = eval_rate_med_f7;
            f_range = f7_range;
     elseif fname == 3 || fname == 4
            f_x_med =  f_x_med_f8;
            sigma_med  =    sigma_med_f8;
            T_med_matrix = T_med_f8;
            t_med_matrix = t_med_f8;
            eval_rate_med = eval_rate_med_f8;
            error_array_med = error_array_med_f8;
            
            T_matrix =  T_f8;
            four_prob = four_prob_med_f8;
            eval_rate = eval_rate_med_f8;
            f_range = f8_range;
%       elseif fname == 4
%             f_x_med =  f_x_med_f9;
%             sigma_med  =    sigma_med_f9;
%             T_matrix =  T_f9;
%             four_prob = four_prob_med_f9;
%             eval_rate = eval_rate_med_f9;
     end
    
    for i = 1:1:STRATEGY_NUM
        % Subplot row and col number of fig1
        
        %==========================
        % Assign reltiave model error in median runs
        %==========================
        T = squeeze(T_med_matrix(i,para,:));
        if i ~=  1      
            t = t_med_matrix(i,para);
            GP_error = squeeze(error_array_med(i,para,:));
            range_t = TRAINING_SIZE+2:1:t;
        end
        
        
        T_array = squeeze(T_matrix(i,para,:));
        range = 1:median(T_matrix(i,para,:));
        
        f_x_array = squeeze(f_x_med(i,para,:));
        sigma_array = squeeze(sigma_med(i,para,:));
        
       
        if i == 1
            lambda = 0;
        elseif i == 2
            lambda = 1;
        elseif i == 3
            lambda = 10;
        elseif i == 4
            lambda = 20;
        elseif i == 5
            lambda = 40;
        end
        
        figure(FIG_NUM);
        legend('-DynamicLegend'); 
        hold on;
        mu = ceil(lambda/4);
        if lambda==0
             d = sprintf('no model');
        elseif lambda == 1
            d = sprintf('(1/1,1)');
            ds = sprintf('(1/1,1) [S]');
        else
            d = sprintf('(%d/%d,%d)',mu,mu,lambda);
            ds = sprintf('(%d/%d,%d) [S]',mu,mu,lambda);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1.objective function [row 1]
        if fname == 2
            subplot(subplot_ROW,subplot_COL,fname+2);
        elseif fname >= 3
            subplot(subplot_ROW,subplot_COL,fname-1);
        else
            subplot(subplot_ROW,subplot_COL,fname);
        end
        h = histogram(T_array,'DisplayName',d);hold on;
        if(fname == 1)
            h.BinWidth = 45;
        % forth col    
        elseif(fname == 2)
            h.BinWidth = 40;
        % second col
        elseif(fname == 3)
            h.BinWidth = 60;
         % third col
        elseif(fname == 4)
            h.BinWidth = 35;
        end
        ylim([0 40]);

        if(lambda==40)
            legend({'no model','(1/1,1)','(3/3,10)','(5/5,20)','(10/10,40)'},'Fontsize',13,'Interpreter','latex');
            if(fname == 1)
                title(sprintf('spheres ($n=%d$, $\\alpha=%.2f$)',n,f_range(para)),'Fontsize',13,'interpreter','latex');
            elseif(fname == 2)
                title(sprintf('quartic function ($n=%d$, $\\gamma=%.2f$)',n,f_range(para)),'Fontsize',13,'interpreter','latex');
            elseif(fname == 3 || fname == 4)
                title(sprintf('ellipsoids ($n=%d$, $\\beta=%.2f$)',n,f_range(para)),'Fontsize',13,'interpreter','latex');
            end
        end
        hLegend = findobj(gcf, 'Type', 'Legend');
        set(hLegend,'Fontsize',12,'interpreter','latex');
        xlabel('objective function calls','Fontsize',13,'interpreter','latex');
        set(gca, 'Fontsize',13);       
        box on

%         ylim([0 0.3]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2.objective function [row 2]
        
        if fname == 2
            subplot(subplot_ROW,subplot_COL,fname+2+subplot_COL);
        elseif fname >= 3
            subplot(subplot_ROW,subplot_COL,fname-1+subplot_COL);
        else
            subplot(subplot_ROW,subplot_COL,fname+subplot_COL);
        end
        plot(range,f_x_array(range),'DisplayName',d);hold on;
        hLegend = findobj(gcf, 'Type', 'Legend');
        set(hLegend,'Fontsize',12,'interpreter','latex');
        if(fname==1)
            ylabel('objective function value','Fontsize',13,'interpreter','latex');
        end
        xlabel('objective function calls','Fontsize',13,'interpreter','latex');
        set(gca, 'YScale', 'log','Fontsize',13);
        ylim([10^(-8) 10^2]);
        yticks([10^-8 10^-6 10^-4 10^-2 10^2]);
        yticklabels({'10^{-8}','10^{-6}','10^{-4}','10^{-2}', '10^{2}'});
        legend('-DynamicLegend'); 
        legend('show');
        box on

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % 3.sigma [row 3]
        if fname == 2
            subplot(subplot_ROW,subplot_COL,fname+2+subplot_COL*2);
        elseif fname >= 3
            subplot(subplot_ROW,subplot_COL,fname-1+subplot_COL*2);
        else
            subplot(subplot_ROW,subplot_COL,fname+subplot_COL*2);
        end
        plot(range,sigma_array(range),'DisplayName',d);hold on;
        hLegend = findobj(gcf, 'Type', 'Legend');
        set(hLegend,'Fontsize',12,'interpreter','latex');
        if(fname==1)
            ylabel('step size','Fontsize',13,'interpreter','latex');
        end
        xlabel('objective function calls','Fontsize',13,'interpreter','latex');
        set(gca, 'YScale', 'log','Fontsize',13);
        ylim([10^(-36) 10^4]);
        yticks([10^(-36) 10^(-26) 10^(-16) 10^-6 10^0 10^4]);
        yticklabels({'10^{-36}','10^{-26}','10^{-16}','10^{-6}','10^{0}' , '10^{4}'});
%         yticks([10^-12 10^-8 10^-4 10^0 10^4]);
%         yticklabels({'10^{-12}','10^{-8}','10^{-4}','10^{0}', '10^{4}'});
        legend('-DynamicLegend'); 
        legend('show');
        box on


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4. GP error [row 4]
         if fname == 2
            subplot(subplot_ROW,subplot_COL,fname+2+subplot_COL*3);
        elseif fname >= 3
            subplot(subplot_ROW,subplot_COL,fname-1+subplot_COL*3);
        else
            subplot(subplot_ROW,subplot_COL,fname+subplot_COL*3);
        end
        %==========================
        % Plot relative model error 
        %==========================
        if lambda == 1 || lambda == 20
            plot(range_t,GP_error(range_t),'DisplayName',d);hold on;
            plot(range_t,exp(conv(log(GP_error(range_t)), kernel, 'same')),'LineWidth',1,'DisplayName',ds);hold on;
        end
        hLegend = findobj(gcf, 'Type', 'Legend');
        set(hLegend,'Fontsize',12,'interpreter','latex','NumColumns',2);
        if(fname==1)
            ylabel('relative model error','Fontsize',13,'interpreter','latex');
        end
        xlabel('iterations','Fontsize',13,'interpreter','latex');
        set(gca, 'YScale', 'log','Fontsize',13);
        ylim([10^(-5) 10^5]);
        legend('-DynamicLegend'); 
        legend('show');
        box on

%         saveas(gcf,sprintf('step_behaviour.fig')); 
        saveas(gcf,sprintf('step_behaviour2.fig')); 

        
    end   
end


