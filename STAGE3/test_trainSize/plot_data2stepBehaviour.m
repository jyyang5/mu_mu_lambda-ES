% test file
% using 5 different obejctive functions plot 5 seperate graphs
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



% For compact subplots 
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.08], [0.05 0.08], [0.05 0.08]);
if ~make_it_tight,  clear subplot;  end

% GP smooth 
window_length = 40;
kernel = exp(-(-3*window_length:3*window_length).^2/window_length^2/2);
kernel = kernel/sum(kernel);        % Normalized   


para = 1;

subplot_ROW = 3;
subplot_COL = 4; 
STRATEGY_NUM = 5;
FIG_NUM = 100;
n = 16;

for fname = 1:1:4
    if fname == 1
            f_x_med =  f_x_med_f6;
            sigma_med  =  sigma_med_f6;
            T_array = 
            T_matrix =  T_f6;
            four_prob = four_prob_med_f6;
            eval_rate = eval_rate_med_f6;
            f_range = f6_range;
    elseif fname == 2
            f_x_med =  f_x_med_f7;
            sigma_med  =    sigma_med_f7;
            T_matrix =  T_f7;
            four_prob = four_prob_med_f7;
            eval_rate = eval_rate_med_f7;
            f_range = f7_range;
     elseif fname == 3
            f_x_med =  f_x_med_f8;
            sigma_med  =    sigma_med_f8;
            T_matrix =  T_f8;
            four_prob = four_prob_med_f8;
            eval_rate = eval_rate_med_f8;
            f_range = f8_range;
      elseif fname == 4
            f_x_med =  f_x_med_f9;
            sigma_med  =    sigma_med_f9;
            T_matrix =  T_f9;
            four_prob = four_prob_med_f9;
            eval_rate = eval_rate_med_f9;
     end
    
    for i = 1:1:STRATEGY_NUM
        % Subplot row and col number of fig1
        

        
    
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
             d = sprintf('(1+1)-ES');
        elseif lambda == 1
            d = sprintf('GP-(1+1)-ES');
        else
            d = sprintf('GP-(1+1)_{(%d/%d,%d)}-ES',mu,mu,lambda);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1.objective function [row 1]
        subplot(subplot_ROW,subplot_COL,fname);
        h = histogram(T_array,'DisplayName',d);hold on;
        if(fname == 1)
            h.BinWidth = 150;
        elseif(fname == 2)
            h.BinWidth = 200;
        elseif(fname == 3)
            h.BinWidth = 600;
        elseif(fname == 4)
            h.BinWidth = 150;
        end
        ylim([0 80]);

        if(lambda==40)
            legend({'(1+1)-ES','(1,1)','(3/3,10)','(5/5,20)','(10/10,40)'},'interpreter','latex');
            if(fname == 1)
                title(sprintf('spheres (n=%d, $\beta$=%.2f)',n,f_range(para)),'fontsize',15,'interpreter','latex');
            elseif(fname == 2)
                title(sprintf('quartic function (n=%d, \\beta=%.2f)',n,f_range(para)),'fontsize',15,'interpreter','latex');
            elseif(fname == 3)
                title(sprintf('ellipsoids (n=%d, $\\beta$=%.2f)',n,f_range(para)),'fontsize',15,'interpreter','latex');
            elseif(fname == 4)
                title(sprintf("Schwefel's function (n=%d)",n),'fontsize',15,'interpreter','latex');
            end
        end
        hLegend = findobj(gcf, 'Type', 'Legend');
        set(hLegend, 'Fontsize', 8,'interpreter','latex');
        xlabel('number of objective function calls','FontSize',15,'interpreter','latex');
%         ylim([0 0.3]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2.objective function [row 2]
        subplot(subplot_ROW,subplot_COL,fname+subplot_COL);
        plot(range,f_x_array(range),'DisplayName',d);hold on;
        hLegend = findobj(gcf, 'Type', 'Legend');
        set(hLegend, 'Fontsize', 8,'interpreter','latex');
  
        if(fname==1)
            ylabel('objective function value','FontSize',15,'interpreter','latex');
        end
        xlabel('objective function calls','FontSize',15,'fontname','interpreter','latex');
        set(gca, 'YScale', 'log');
        ylim([10^(-10) 100]);
        legend('-DynamicLegend'); 
        legend('show');
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % 3.sigma [row 3]
        subplot(subplot_ROW,subplot_COL,fname+subplot_COL*2);
        plot(range,sigma_array(range),'DisplayName',d);hold on;
        hLegend = findobj(gcf, 'Type', 'Legend');
        set(hLegend, 'Fontsize', 8,'interpreter','latex');
        if(fname==1)
            ylabel('step size','FontSize',15,'interpreter','latex');
        end
        xlabel('objective function calls','FontSize',15,'interpreter','latex');
        set(gca, 'YScale', 'log');
        ylim([10^(-12) 1]);
        legend('-DynamicLegend'); 
        legend('show');
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4. GP error [row 4]
        subplot(subplot_ROW,subplot_COL,fname+subplot_COL*2);
        plot(range,sigma_array(range),'DisplayName',d);hold on;
        hLegend = findobj(gcf, 'Type', 'Legend');
        set(hLegend, 'Fontsize', 8,'interpreter','latex');
        if(fname==1)
            ylabel('step size','FontSize',15,'interpreter','latex');
        end
        xlabel('objective function calls','FontSize',15,'interpreter','latex');
        set(gca, 'YScale', 'log');
        ylim([10^(-12) 1]);
        legend('-DynamicLegend'); 
        legend('show');
        
%         if(fname<4)
%             subplot(subplot_ROW,subplot_COL,fname+15);
%             if lambda==0
%                 plot(1:1:T_med,sigma_star_matrix_med(1:T_med),'DisplayName',d);hold on;
%             else 
%                 sigma_star_range1 = sigma_star_matrix_med(1:t_start);
%                 sigma_star_range2 = sigma_star_matrix_med(t_start+1:t_med);
%                 plot([t_range1 t_range2], [sigma_star_range1 sigma_star_range2],'DisplayName',d);hold on;
%             end
%             if(fname==1)
%                 ylabel('normalized step size \sigma*','FontSize',15);%
%             end
%             xlabel('objective function calls','FontSize',15); 
%             set(gca, 'YScale', 'log');
%         %     title('normalized step size \sigma*','FontSize',20);
%              ylim([0.05 15]);
%             legend('-DynamicLegend'); 
%             legend('show');
%         end
        saveas(gcf,'merged_plot_v2.fig'); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % 5. histogram normalized fitness gain per iteration
%         figure(500);
%         subplot(1,subplot_COL,fname);
%         % No bar colour, bold binwidth
%         if(T_med == 5000) % early stopping
%             h2 = histogram(nonzeros(delta_matrix_temp),'Normalization','probability','DisplayName',d);hold on;
%         else
%             h2 = histogram(nonzeros(delta_matrix_temp),'Normalization','probability','DisplayName',d);hold on;
%         end
%         xlim([-3 3]);
%         if(lambda==0)
%             h2.EdgeColor= [0  0.4470 0.7410];
%         elseif(lambda == 10)
%             h2.EdgeColor= [0.8500  0.3250  0.0980];
%         elseif(lambda==20)
%             h2.EdgeColor= [0.9290  0.6940  0.1250];
%         elseif(lambda==40)
%             h2.EdgeColor= [0.4940  0.1840  0.5560];
%         end
%         h2.LineWidth=1.2;
%         h2.BinWidth = 0.2;
        
%         if(fname == 1 || fname == 2 || fname == 3)
%             figure(201);
%             subplot(2,subplot_COL,fname);
%             value = h2.Values;		% height of the bar
%             width = h2.BinWidth;				% width of the bar
%             range = h2.BinLimits;		% [startX endX]
%             disp(length(range(1)+width/2:width:range(2)-width/2)==length(value));
%             % Did not do the exact range right ends at [range(1)+width/2:width:range(2)-width/2] 
%             plot((range(1)+width/2):width:(range(2)),value,'DisplayName',d);hold on;
%             if(fname==1)
%                 ylabel('Probability','FontSize',15);%
%             end
%             xlabel('Normalized fitness gain','FontSize',15); 
%             % set(gca, 'YScale', 'log');
%             % title('step size \sigma','FontSize',20);
%             legend('-DynamicLegend'); 
%             legend('show');
%             xlim([-3 3]);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Subplot row and col number of fig2
%             subplot_ROW = 4;
%             subplot_COL = 5; 
% 
%             xNameSprintf = sprintf('onvergence rate');
%             xLimit = [0 1];
%             plot_pdf(convergence_array,T_med,201,subplot_ROW,subplot_COL,fname,lambda,xNameSprintf,xLimit);
%             xNameSprintf = sprintf('success rate');
%             xLimit = [0 1];
%             plot_pdf(success_med,T_med,201,subplot_ROW,subplot_COL,fname+subplot_COL*2,lambda,xNameSprintf,xLimit);
% %             xNameSprintf = sprintf('fitGain per iteration');
% %             xLimit = [0 1];
%     %         plot_pdf(delta_matrix_temp,T_med,300,3,5,fname+10,lambda,xNameSprintf,xLimit);
%             legend('-DynamicLegend'); 
%             legend('show');
%         end
        

        
    end   
end




function plot_pdf(data,T_med,figureName,fig_row,fig_col,fig_index,lambda,xNameSprintf,xLimit)
%  Input: 
%       handler for histogram
%  Plot:
    % No bar colour, bold binwidth
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot histogram 
    figure(figureName);
    subplot(fig_row,fig_col,fig_index);
    if lambda==0
        d = sprintf('(1+1)-ES');
    else
        mu = ceil(lambda/4);
        d = sprintf('(%d/%d,%d)-ES',mu,mu,lambda);
    end
    if(T_med == 2000) % early stopping
        h2 = histogram(nonzeros(data),'Normalization','probability','DisplayName',d);hold on;
    else
        h2 = histogram(nonzeros(data),'Normalization','probability','DisplayName',d);hold on;
    end
    
    h2.LineWidth=0.5;
    
    if(lambda==0)
        h2.EdgeColor= [0  0.4470 0.7410];
    elseif(lambda == 10)
        h2.EdgeColor= [0.8500  0.3250  0.0980];
    elseif(lambda==20)
        h2.EdgeColor= [0.9290  0.6940  0.1250];
    elseif(lambda==40)
        h2.EdgeColor= [0.4940  0.1840  0.5560];
    end
    if(rem(fig_index,fig_col)==1)
        ylabel('probability','FontSize',15);%
    end
    % Set titles for the first row 
    if(fig_index == 1)
        d3 =sprintf('linear sphere');
        title(d3,'fontsize',15);
    elseif(fig_index == 2)
        d3 =sprintf('quadratic sphere');
        title(d3,'fontsize',15);
    elseif(fig_index == 3)
        d3 =sprintf('cubic sphere');
        title(d3,'fontsize',15);
    elseif(fig_index == 4)
        d3 =sprintf('Schwefel function');
        title(d3,'fontsize',15);
    elseif(fig_index == 5)
        d3 =sprintf('quartic function');
        title(d3,'fontsize',15);
    end
    ylim([0 0.4]);
    
    legend('-DynamicLegend'); 
    legend('show');
    xlim(xLimit);
    xlabel(xNameSprintf,'FontSize',15); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot pdf curve
    figure(figureName);
    subplot(fig_row,fig_col,fig_index+fig_col);
    value = h2.Values;		% height of the bar
    width = h2.BinWidth;				% width of the bar
    range = h2.BinLimits;		% [startX endX]
    % Did not do the exact range right ends at [range(1)+width/2:width:range(2)-width/2] 
    plot((range(1)+width/2):width:(range(2)),value,'DisplayName',d);hold on;
    if(rem(fig_index,fig_col)==1)
        ylabel('probability','FontSize',15);%
    end
    xlim(xLimit);
    xlabel(xNameSprintf,'FontSize',15); 
    % set(gca, 'YScale', 'log');
    % title('step size \sigma','FontSize',20);
    legend('-DynamicLegend'); 
    legend('show');
    

end