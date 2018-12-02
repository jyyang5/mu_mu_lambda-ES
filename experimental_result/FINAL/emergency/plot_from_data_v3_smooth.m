% test file
% using 5 different obejctive functions plot 5 seperate graphs
% NOTE: DO NOT CHANGE THE NUMBER WHEN CALL fun_graph_merged
%       OR CANNOT SAVE ITERTAION DATA
% MODIFICATION 
%       take median run values (mediean objective function calls)
%       Plot 1
%           Row 1: histogram of objective function calls
%           Row 2: histogram of objective function calls
%           Row 3: histogram of objective function calls
%           Row 4: normalized step size [mediuam run]
%       Plot 2
%           Row 1: normalized convergence rate
%           Row 2: normalized success rate [good step]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = @(x) (x'*x)^(1/2);
f2 = @(x) (x'*x);
f3 = @(x) (x'*x)^(3/2);

close all;


% For compact subplots 
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.05 0.02], [0.05 0.02]);
if ~make_it_tight,  clear subplot;  end

% GP smooth 
window_length = 40;
kernel = exp(-(-3*window_length:3*window_length).^2/window_length^2/2);
kernel = kernel/sum(kernel);        % Normalized   

% NUM_OF_RUNS = 1000;
% % NUM_OF_RUNS = 2;
% TRAINING_SIZE = 40;
% LENGTH_SCALE = 8;


% 
% t_array = zeros(5,4,NUM_OF_RUNS,1);                                           % # of iterations for the stop creteria
% sigma_matrix = zeros(5,4,NUM_OF_RUNS,10000);                                  % store all sigma
% T_array = zeros(5,4,NUM_OF_RUNS,1);                                           % # of objective function evaluations for the stop creteria
% f_x_matrix = zeros(5,4,NUM_OF_RUNS,10000);                                    % store all fx
% convergence_rate_array = zeros(5,4,NUM_OF_RUNS,1);                            % convergence rate 
% GP_error_matrix = zeros(5,4,NUM_OF_RUNS,10000);                               % store similar to noise-to-signal ratio
% sigma_star_matrix = zeros(5,4,NUM_OF_RUNS,10000);                             % normalized step size 
% success_rate_array = zeros(5,4,NUM_OF_RUNS,1);                                % success rate 
% delta_matrix = zeros(5,4,NUM_OF_RUNS,10000);                                  % each [i,j] stores a delta array 

% t_med = median(t_array);
% t_max = max(t_array);
% t_min = min(t_array);
% sigma_matrix_med = squeeze(median(sigma_matrix));
% T_med = median(T_array);
% f_x_med = squeeze(median(f_x_matrix));
% convergence_med = median(convergence_rate_array);
% GP_error_matrix_med = squeeze(median(GP_error_matrix));
% sigma_star_matrix_med = squeeze(median(sigma_star_matrix));
% success_med = median(success_rate_array);
% delta_matrix_med = squeeze(median(delta_matrix));
lambda_array = [0,10,20,40];

close all;
for fname = 1:1:5
    for i = 1:1:length(lambda_array)
        % Subplot row and col number of fig1
        subplot_ROW = 4;
        subplot_COL = 5; 
%         t_med = median(t_array(fname,i,:));
        t_max = max(t_array(fname,i,:));
        t_min = min(t_array(fname,i,:));
%         sigma_matrix_med = median(squeeze(sigma_matrix(fname,i,:,:)));
        
        % Find the index of median T 
        T_array_temp = T_array(fname,i,:);
        sorted_T = sort(T_array_temp);
        med_index = find(T_array_temp == sorted_T(ceil(length(sorted_T)/2)));
        med_index = med_index(1);   % randomly pick one
        T_med = T_array(fname,i,med_index);
        
        t_med = squeeze(t_array(fname,i,med_index));
        sigma_matrix_med = transpose(squeeze(sigma_matrix(fname,i,med_index,:)));
        f_x_med = transpose(squeeze((f_x_matrix(fname,i,med_index,:))));
        
        convergence_array = squeeze((convergence_rate_array(fname,i,:)));
        GP_error_matrix_med = transpose(squeeze(GP_error_matrix(fname,i,med_index,:)));
        sigma_star_matrix_med = transpose(squeeze(sigma_star_matrix(fname,i,med_index,:)));
        success_med = success_rate_array(fname,i,:);
        delta_matrix_temp =squeeze(delta_matrix(fname,i,:,:));
%         f_x_med = median(squeeze((f_x_matrix(fname,i,:,:))));
%         convergence_array = squeeze((convergence_rate_array(fname,i,:)));
%         GP_error_matrix_med = median(squeeze(GP_error_matrix(fname,i,:,:)));
%         sigma_star_matrix_med = median(squeeze(sigma_star_matrix(fname,i,:,:)));
%         success_med = success_rate_array(fname,i,:);
%         delta_matrix_temp =squeeze(delta_matrix(fname,i,:,:));
        
        lambda =  lambda_array(i);
        figure(200);
        legend('-DynamicLegend'); 
        hold on;
        mu = ceil(lambda/4);
        if lambda==0
            d = sprintf('(1+1)-ES');
        else
            d = sprintf('(%d/%d,%d)-ES',mu,mu,lambda);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1.objective function [row 1]
        subplot(subplot_ROW,subplot_COL,fname);
        h = histogram(T_array(fname,i,:),'Normalization','probability','DisplayName',d);hold on;
        if(fname == 1)
            h.BinWidth = 10;
        elseif(fname == 2)
            h.BinWidth = 5;
        elseif(fname == 3)
            h.BinWidth = 5;
        elseif(fname == 4)
            h.BinWidth = 60;
        elseif(fname == 5)
            h.BinWidth = 35; 
        end
%         h.BinWidth = 20;          % Bin width
        if(lambda==40)
            legend({'(1+1)-ES','(3/3,10)-ES','(5/5,20)-ES','(10/10,40)-ES'})
            if(fname == 1)
                d3 =sprintf('linear sphere');
                title(d3,'fontsize',15);
            elseif(fname == 2)
                d3 =sprintf('quadratic sphere');
                title(d3,'fontsize',15);
            elseif(fname == 3)
                d3 =sprintf('cubic sphere');
                title(d3,'fontsize',15);
            elseif(fname == 4)
                d3 =sprintf('Schwefel;s function');
                title(d3,'fontsize',15);
            elseif(fname == 5)
                d3 =sprintf('quartic function');
                title(d3,'fontsize',15); 
            end
        end
        if(fname==1)
            ylabel('probability','FontSize',15);
        end
        xlabel('objective function calls','FontSize',15); 
        ylim([0 0.3]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2.objective function [row 2]
        subplot(subplot_ROW,subplot_COL,fname+5);
        if lambda==0
            t_start = TRAINING_SIZE+1+2;
            plot(1:1:t_med,f_x_med(1:T_med),'DisplayName',d);hold on;
        else 
            t_start = ceil(TRAINING_SIZE/lambda);
            fx_range1 = f_x_med(1:t_start);
            fx_range2 = f_x_med(t_start+1:t_med);
            t_range1 = 1:lambda:lambda*t_start;
            t_range2 = lambda*t_start+1:lambda*t_start+length(fx_range2);
            plot([t_range1 t_range2], [fx_range1 fx_range2],'DisplayName',d);hold on;% mml with GP
        end
        if(fname==1)
            ylabel('objective function value','FontSize',15);%
        end
        xlabel('objective function calls','FontSize',15); 
        set(gca, 'YScale', 'log');
        ylim([10^(-10) 50]);
        legend('-DynamicLegend'); 
        legend('show');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3.GP error [row 3]
        subplot(subplot_ROW,subplot_COL,fname+10);
        if lambda==0
            plot(t_start:1:T_med,GP_error_matrix_med(t_start:T_med),'DisplayName',d);hold on;
%             smoothed_GP = smoothdata(GP_error_matrix_med(t_start:T_med),'gaussian',40);
            smoothed_GP = exp(conv(log(GP_error_matrix_med(t_start:T_med)), kernel, 'same'));
            d1 = sprintf('(1+1)-ES[smoothed]');
            plot(t_start:1:T_med,smoothed_GP,'DisplayName',d1,'LineWidth',2);hold on;
        else 
%             GP_error_range1 = GP_error_matrix_med(1:t_start);
            GP_error_range2 = GP_error_matrix_med(t_start+1:t_med);
            plot(t_range2,GP_error_range2,'DisplayName',d);hold on;
%             plot([t_range1 t_range2], [GP_error_range1 GP_error_range2],'DisplayName',d);hold on;
            % GP smoothed
%             smoothed_GP_range1 = smoothdata(GP_error_range1,'gaussian',40);
%             smoothed_GP_range2 = smoothdata(GP_error_range2,'gaussian',40);
            smoothed_GP_range2 = exp(conv(log(GP_error_range2), kernel, 'same'));

            d1 = sprintf('(%d/%d,%d)-ES[smoothed]',mu,mu,lambda);
            plot(t_range2, smoothed_GP_range2,'DisplayName',d1,'LineWidth',2);hold on;
%             plot([t_range1 t_range2],[smoothed_GP_range1 smoothed_GP_range2],'DisplayName',d1,'LineWidth',2);hold on;
        end
        if(fname==1)
            ylabel('relative model error','FontSize',15);%
        end
        xlabel('objective function calls','FontSize',15); 
        set(gca, 'YScale', 'log');
        ylim([0.01 10^4]);
        legend('-DynamicLegend'); 
        legend('show');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4.normalized step size [row 4](only makes sense for sphere functions)
        if(fname<4)
            subplot(subplot_ROW,subplot_COL,fname+15);
            if lambda==0
                plot(1:1:T_med,sigma_star_matrix_med(1:T_med),'DisplayName',d);hold on;
            else 
                sigma_star_range1 = sigma_star_matrix_med(1:t_start);
                sigma_star_range2 = sigma_star_matrix_med(t_start+1:t_med);
                plot([t_range1 t_range2], [sigma_star_range1 sigma_star_range2],'DisplayName',d);hold on;
            end
            if(fname==1)
                ylabel('normalized step size \sigma*','FontSize',15);%
            end
            xlabel('objective function calls','FontSize',15); 
            set(gca, 'YScale', 'log');
        %     title('normalized step size \sigma*','FontSize',20);
             ylim([0.05 15]);
            legend('-DynamicLegend'); 
            legend('show');
        end
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
            subplot_ROW = 4;
            subplot_COL = 5; 

            xNameSprintf = sprintf('onvergence rate');
            xLimit = [0 1];
            plot_pdf(convergence_array,T_med,201,subplot_ROW,subplot_COL,fname,lambda,xNameSprintf,xLimit);
            xNameSprintf = sprintf('success rate');
            xLimit = [0 1];
            plot_pdf(success_med,T_med,201,subplot_ROW,subplot_COL,fname+subplot_COL*2,lambda,xNameSprintf,xLimit);
%             xNameSprintf = sprintf('fitGain per iteration');
%             xLimit = [0 1];
    %         plot_pdf(delta_matrix_temp,T_med,300,3,5,fname+10,lambda,xNameSprintf,xLimit);
            legend('-DynamicLegend'); 
            legend('show');
            saveas(gcf,'success_convergence_v2.fig');  
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