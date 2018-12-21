% Objective [for normalized step size to update step size]
% 1 Plot
%     Scatter plot: [# of iterations] vs. [success rate]
%     NOTE: each individual run uses one SIGMA_STAR 
% 2. Save
%     Plots & Data     
%
% difficulty: 
%            first few trainijng iterations
%            histogram does not show data properly for sigmaStar and model
%            error
%            save file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = fun_multi_run_SIGMA_STAR(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE,FIGURE_NUM,SIGMA_STAR_array,subplot_ROW,subplot_COL,ROW_num)
%Input:
%   fname:            an index 
%                        1 for linear
%                        2 for quadratic 
%                        3 for cubic 
%                        4 for schwefel
%                        5 for quartic
%    NUM_OF_RUNS      # of runs to take median 
%    lambda           0 for (1+1)-ES
%                     10 for (3/3,10)-ES
%                     20 for (5/5,20)-ES
%                     40 for (10/10,40)-ES
%    TRAINING_SIZE    GP training size
%    LENGTH_SCALE     length scale factor for GP
%    SIGMA_STAR_array array of normalized step size used to update step size
%    subplot_ROW      # of lambda used (each lambda per row) 
%    subplot_COL      # of test problems used (each test problem per col) 
%    ROW_num          specify the row number according to the lambda
% 
%Return:
% 	1.t_array
%   2.sigma_matrix
%   3.T_array
%   4.f_x_matrix
%   5.convergence_rate_array
%   6.GP_error_matrix
%   7.sigma_star_matrix
%   8.success_rate_array
%   9.delta_matrix
   
% iteration number for [mmlWithGP,mmlNoGP,1+1WithGP,1+1NoGP]    
NUM_OF_ITERATIONS = 4000;           % MAX # of iterations       
LEN_SIGMA_STAR = length(SIGMA_STAR_array);

% lambda = 10;
% mu = 3;
n = 10;             % dimension of the problem



% For compact subplots 
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.05], [0.1 0.05], [0.1 0.05]);
if ~make_it_tight,  clear subplot;  end

% GP smooth 
window_length = 40;
kernel = exp(-(-3*window_length:3*window_length).^2/window_length^2/2);
kernel = kernel/sum(kernel);        % Normalized    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Def of variables

% (mu/mu,lambda)-ES with GP
t_array = zeros(LEN_SIGMA_STAR,NUM_OF_RUNS,1);                    % # of iterations for the stop creteria
sigma_matrix = zeros(LEN_SIGMA_STAR,NUM_OF_RUNS,10000);           % store all sigma
T_array = zeros(LEN_SIGMA_STAR,NUM_OF_RUNS,1);                    % # of objective function evaluations for the stop creteria
f_x_matrix = zeros(LEN_SIGMA_STAR,NUM_OF_RUNS,10000);             % store all fx
success_rate_array = zeros(LEN_SIGMA_STAR,NUM_OF_RUNS,1);         % success rate 


med_index = zeros(LEN_SIGMA_STAR,1);
sorted_T = zeros(NUM_OF_RUNS,1);
sigma_matrix_med = zeros(LEN_SIGMA_STAR,10000);
f_x_med = zeros(LEN_SIGMA_STAR,10000);
t_med = zeros(LEN_SIGMA_STAR,1);
success_rate_med = zeros(LEN_SIGMA_STAR,1);

% (mu/mu,lambda)-ES
t_array_noGP = zeros(LEN_SIGMA_STAR,NUM_OF_RUNS,1);                    % # of iterations for the stop creteria
sigma_matrix_noGP = zeros(LEN_SIGMA_STAR,NUM_OF_RUNS,10000);           % store all sigma
T_array_noGP = zeros(LEN_SIGMA_STAR,NUM_OF_RUNS,1);                    % # of objective function evaluations for the stop creteria
f_x_matrix_noGP = zeros(LEN_SIGMA_STAR,NUM_OF_RUNS,10000);             % store all fx
success_rate_array_noGP = zeros(LEN_SIGMA_STAR,NUM_OF_RUNS,1);         % success rate 


med_index_noGP = zeros(LEN_SIGMA_STAR,1);
sorted_T_noGP = zeros(LEN_SIGMA_STAR,1);
sigma_matrix_med_noGP = zeros(LEN_SIGMA_STAR,10000);
f_x_med_noGP = zeros(LEN_SIGMA_STAR,10000);
t_med_noGP = zeros(LEN_SIGMA_STAR,1);
success_rate_med_noGP = zeros(LEN_SIGMA_STAR,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replicates for NUM_OF_RUNS
mu = ceil(lambda/4);
for j = 1:1:LEN_SIGMA_STAR
    SIGMA_STAR = SIGMA_STAR_array(j);

    for i = 1:NUM_OF_RUNS
        x0 = randn(n,mu);
        a = mml_GP_SIGMA_STAR(fname,x0,SIGMA_STAR,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE);
        t_array(j,i) = cell2mat(a(1));
        sigma_matrix(j,i,:) = cell2mat(a(4));
        T_array(j,i) = cell2mat(a(5));
        f_x_matrix(j,i,:) = cell2mat(a(6));
        success_rate_array(j,i) = cell2mat(a(10));
        
        b = mml_noGP_SIGMA_STAR(fname,x0,SIGMA_STAR,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE);
        t_array_noGP(j,i) = cell2mat(b(1));  
        sigma_matrix_noGP(j,i,:) = cell2mat(b(4));
        T_array_noGP(j,i) = cell2mat(b(5));
        f_x_matrix_noGP(j,i,:) = cell2mat(b(6));
        success_rate_array_noGP(j,i) = cell2mat(b(10));

    end 
    
    % Take median run
    sorted_T = sort(T_array(j,:));
    temp_index = find(T_array(j,:) == sorted_T(ceil(length(sorted_T)/2)));
    med_index(j) = temp_index(1);
    sigma_matrix_med(j,:) = sigma_matrix(j,med_index(j),:);
    f_x_med(j,:) = f_x_matrix(j,med_index(j),:);
    t_med(j) = t_array(j,med_index(j));
    success_rate_med(j) = success_rate_array(j,med_index(j));
    
    sorted_T_noGP = sort(T_array_noGP(j,:));
    temp_index = find(T_array_noGP(j,:) == sorted_T_noGP(ceil(length(sorted_T_noGP)/2)));
    med_index_noGP(j) = temp_index(1);
    sigma_matrix_med_noGP(j,:) = sigma_matrix_noGP(j,med_index_noGP(j),:);
    f_x_med_noGP(j,:) = f_x_matrix_noGP(j,med_index_noGP(j),:);    
    t_med_noGP(j) = t_array_noGP(j,med_index_noGP(j));
    success_rate_med_noGP(j) = success_rate_array_noGP(j,med_index_noGP(j));

    % counter
    fprintf('row %d finished \n',j)

end

    





        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graph


% For f(x), sigma, sigmaStar, GP_error over objective function evaluations
figure(FIGURE_NUM);
legend('-DynamicLegend'); 
hold on;
subplot(subplot_ROW,subplot_COL,(ROW_num-1)*subplot_COL+fname);
    mu = ceil(lambda/4);
    d = sprintf('GP-(%d/%d,%d)-ES',mu,mu,lambda);
    d_noGP = sprintf('(%d/%d,%d)-ES',mu,mu,lambda);
    scatter(success_rate_med,t_med,'DisplayName',d);hold on;
    scatter(success_rate_med_noGP,t_med_noGP,'DisplayName',d_noGP);hold on;
    if(fname==1)
        ylabel('number of iterations','FontSize',15);%
    end
    xlabel('success rate','FontSize',15); 

    if lambda==10 
        if(fname == 1)
            dt =sprintf('linear sphere');
            title(dt,'fontsize',15);
        elseif(fname == 2)
            dt =sprintf('quadratic sphere');
            title(dt,'fontsize',15);
        elseif(fname == 3)
            dt =sprintf('cubic sphere');
            title(dt,'fontsize',15);
        end
    end
    legend('-DynamicLegend'); 
    legend('show');
    saveas(gcf,'iteration_success_SIGMA_STAR.fig'); 

% 
% scatter();
% for i=1:1:LEN_SIGMA_STAR
%     % Mamually specify 
%     
%    
%         
%     plot(1:t_med(i), f_x_med(i,1:t_med(i)),'DisplayName',d);hold on; % mml with GP
%     plot(1:t_med_noGP(i), f_x_med_noGP(i,1:t_med_noGP(i)),'DisplayName',d_noGP);hold on; % mml with GP
%     if(fname==1)
%         ylabel('objective function value','FontSize',15);%
%     end
%     xlabel('iteration','FontSize',15); 
%     set(gca, 'YScale', 'log');
%     
%     if lambda==10 && i==1
%         if(fname == 1)
%             dt =sprintf('linear sphere');
%             title(dt,'fontsize',15);
%         elseif(fname == 2)
%             dt =sprintf('quadratic sphere');
%             title(dt,'fontsize',15);
%         elseif(fname == 3)
%             dt =sprintf('cubic sphere');
%             title(dt,'fontsize',15);
% %         elseif(fname == 4)
% %             dt =sprintf('Schwefel function');
% %             title(dt,'fontsize',15);
% %         elseif(fname == 5)
% %             dt =sprintf('quartic function');
% %             title(dt,'fontsize',15); 
%         end
%     end
%     
%     legend('-DynamicLegend'); 
%     legend('show');
% end
% 
% 
% % % success rate
% % figure(FIGURE_NUM+1);
% % legend('-DynamicLegend'); 
% % hold on;
% % d = sprintf('GP-(%d/%d,%d)-ES',mu,mu,lambda);
% % d_noGP = sprintf('(%d/%d,%d)-ES',mu,mu,lambda);
% % for i=1:1:LEN_SIGMA_STAR
% %     
% %     if i==1
% %         if(fname == 1)
% %             dt =sprintf('linear sphere');
% %             title(dt,'fontsize',15);
% %         elseif(fname == 2)
% %             dt =sprintf('quadratic sphere');
% %             title(dt,'fontsize',15);
% %         elseif(fname == 3)
% %             dt =sprintf('cubic sphere');
% %             title(dt,'fontsize',15);
% %         elseif(fname == 4)
% %             dt =sprintf('Schwefel function');
% %             title(dt,'fontsize',15);
% %         elseif(fname == 5)
% %             dt =sprintf('quartic function');
% %             title(dt,'fontsize',15); 
% %         end
% %     end
% %     
% %     subplot(subplot_ROW,subplot_COL,(i-1)*subplot_COL+fname);
% %     bar(success_rate_med,'DisplayName',d);hold on; % mml with GP
% %     bar(success_rate_med_noGP,'DisplayName',d_noGP);hold on; % mml with GP
% %     
% % %     plot(1:t_med(i), f_x_med(i,1:t_med(i))
% % %     plot(1:t_med_noGP(i), f_x_med_noGP(i,1:t_med_noGP(i)),'DisplayName',d_noGP);hold on; % mml with GP
% %     if(fname==1)
% %         ylabel('Success rate','FontSize',15);%
% %     end
% % %     xlabel('iteration','FontSize',15); 
% % %     set(gca, 'YScale', 'log');
% %     legend('-DynamicLegend'); 
% %     legend('show');
% % end
% 
% %     
% % % end
% % % if lambda==0
% % %     d = sprintf('(1+1)-ES');
% % % else
% % %     d = sprintf('(%d/%d,%d)-ES',mu,mu,lambda);
% % % end
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 1.objective function [row 1]
% % subplot(subplot_ROW,subplot_COL,fname);
% % h = histogram(T_array,'Normalization','probability','DisplayName',d);hold on;
% % % h.BinWidth = 20;
% 
% % % if(fname==1)
% % %     ylabel('probability','FontSize',15);%
% % % end
% % xlabel('objective function calls','FontSize',15); 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % % 2.objective function [row 2]
% % subplot(subplot_ROW,subplot_COL,fname+5);
% % if lambda==0
% %     t_start = TRAINING_SIZE+1+2;
% %     plot(1:1:t_med,f_x_med(1:T_med),'DisplayName',d);hold on;
% % else 
% %     t_start = ceil(TRAINING_SIZE/lambda);
% %     fx_range1 = f_x_med(1:t_start);
% %     fx_range2 = f_x_med(t_start+1:t_med);
% %     t_range1 = 1:lambda:lambda*t_start;
% %     t_range2 = lambda*t_start+1:lambda*t_start+length(fx_range2);
% %     plot([t_range1 t_range2], [fx_range1 fx_range2],'DisplayName',d);hold on;% mml with GP
% % end
% % if(fname==1)
% %     ylabel('objective function value','FontSize',15);%
% % end
% % xlabel('objective function calls','FontSize',15); 
% % set(gca, 'YScale', 'log');
% % legend('-DynamicLegend'); 
% % legend('show');
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 3.GP error [row 3]
% % subplot(subplot_ROW,subplot_COL,fname+10);
% % if lambda==0
% %     plot(t_start:1:T_med,GP_error_matrix_med(t_start:T_med),'DisplayName',d);hold on;
% %     d1 = sprintf('(1+1)-ES[S]');
% %     plot(t_start:1:T_med,exp(conv(log(GP_error_matrix_med(t_start:T_med)), kernel, 'same')),'DisplayName',d1,'LineWidth',2);hold on;
% % %     plot(t_start:1:T_med,smoothdata(GP_error_matrix_med(t_start:T_med),'gaussian',40),'DisplayName',d1,'LineWidth',2);hold on;
% % else 
% % %     GP_error_range1 = GP_error_matrix_med(1:t_start);
% %     GP_error_range2 = GP_error_matrix_med(t_start+1:t_med);
% %     plot(t_range2,GP_error_range2,'DisplayName',d);hold on;
% %     % GP smoothed
% % %     smoothed_GP_range1 = smoothdata(GP_error_range1,'gaussian',40);
% %     smoothed_GP_range2 = exp(conv(log(GP_error_range2), kernel, 'same'));
% %     d1 = sprintf('(%d/%d,%d)-ES[S]',mu,mu,lambda);
% %     plot(t_range2, smoothed_GP_range2,'DisplayName',d1,'LineWidth',2);hold on;
% % %     plot([t_range1 t_range2],[smoothed_GP_range1 smoothed_GP_range2],'DisplayName',d1,'LineWidth',2);hold on;
% % end
% % if(fname==1)
% %     ylabel('relative model error','FontSize',15);%
% % end
% % xlabel('objective function calls','FontSize',15); 
% % set(gca, 'YScale', 'log');
% % 
% % legend('-DynamicLegend'); 
% % legend('show');
% % % title('Logarithmic relative model error','FontSize',20);
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % 4.normalized step size [row 4](only makes sense for sphere functions)
% % % if(fname<4)
% % %     subplot(subplot_ROW,subplot_COL,fname+15);
% % %     if lambda==0
% % %         plot(1:1:T_med,sigma_star_matrix_med(1:T_med),'DisplayName',d);hold on;
% % %     else 
% % %         sigma_star_range1 = sigma_star_matrix_med(1:t_start);
% % %         sigma_star_range2 = sigma_star_matrix_med(t_start+1:t_med);
% % %         plot([t_range1 t_range2], [sigma_star_range1 sigma_star_range2],'DisplayName',d);hold on;
% % %     end
% % %     if(fname==1)
% % %         ylabel('normalized step size \sigma*','FontSize',15);%
% % %     end
% % %     xlabel('objective function calls','FontSize',15); 
% % %     set(gca, 'YScale', 'log');
% % % 
% % %     legend('-DynamicLegend'); 
% % %     legend('show');
% % % end
% % saveas(gcf,'merged_plot_NO_emergency_v2.fig'); 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% %     % Plot success rate for a good step 
% %     subplot_ROW = 1;
% %     subplot_COL = 5;    
% %     xNameSprintf = sprintf('success rate');
% %     % Plot success rate for a good step
% %     xLimit = [0 1];
% %     plot_pdf(success_rate_array,T_med,FIGURE_NUM+1,subplot_ROW,subplot_COL,fname,lambda,xNameSprintf,xLimit);
% %     legend('-DynamicLegend'); 
% %     legend('show');

function plot_pdf(data,T_med,figureName,fig_row,fig_col,fig_index,lambda,xNameSprintf,xLimit)
%  Input: 
%       handler for histogram
%  Plot:
    % Histogram and pdf function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For compact subplots 
    make_it_tight = true;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.03], [0.05 0.03], [0.05 0.03]);
    if ~make_it_tight,  clear subplot;  end
    
    
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
    val = {};
% val = {t_array,sigma_matrix,T_array,f_x_matrix,convergence_rate_array,GP_error_matrix,sigma_star_matrix,success_rate_array,delta_matrix};
end