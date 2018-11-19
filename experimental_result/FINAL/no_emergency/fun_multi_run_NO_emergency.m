% Objective [for CSA with emergency store emergency rate]
% 1 Plot
%     1.1 relative model noise to lambda, mu = ceil(lambda/4) 
%     1.2 then effect of changing theta
% 2. Save
%     Data     
% save objective function calls data & plot into files
%
% difficulty: 
%            first few trainijng iterations
%            histogram does not show data properly for sigmaStar and model
%            error
%            save file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = fun_multi_run_NO_emergency(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE,DECREASE_FACTOR,FIGURE_NUM)
%Input:
%   fname:          an index 
%                       1 for linear
%                       2 for quadratic 
%                       3 for cubic 
%                       4 for schwefel
%                       5 for quartic
%    NUM_OF_RUNS    # of runs to average
%    lambda         0 for (1+1)-ES
%                   10 for (3/3,10)-ES
%                   20 for (5/5,20)-ES
%                   40 for (10/10,40)-ES
%    TRAINING_SIZE  GP training size
%    LENGTH_SCALE   length scale factor for GP
%    DECREASE_FACTOR rate of step size to be decreased
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
sigma0 = 1;
NUM_OF_ITERATIONS = 2000;        

% lambda = 10;
% mu = 3;
n = 10;

subplot_ROW = 4;
subplot_COL = 5; 

% For compact subplots 
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.05 0.02], [0.05 0.02]);
if ~make_it_tight,  clear subplot;  end

% GP smooth 
window_length = 40;
kernel = exp(-(-3*window_length:3*window_length).^2/window_length^2/2);
kernel = kernel/sum(kernel);        % Normalized    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Def of variables

% (mu/mu,lambda)-ES with GP
t_array = zeros(NUM_OF_RUNS,1);                    % # of iterations for the stop creteria
sigma_matrix = zeros(NUM_OF_RUNS,10000);           % store all sigma
T_array = zeros(NUM_OF_RUNS,1);                    % # of objective function evaluations for the stop creteria
f_x_matrix = zeros(NUM_OF_RUNS,10000);             % store all fx
convergence_rate_array = zeros(NUM_OF_RUNS,1);     % convergence rate 
GP_error_matrix = zeros(NUM_OF_RUNS,10000);        % store similar to noise-to-signal ratio
sigma_star_matrix = zeros(NUM_OF_RUNS,10000);      % normalized step size 
success_rate_array = zeros(NUM_OF_RUNS,1);         % success rate 
delta_matrix = zeros(NUM_OF_RUNS,10000);           % each [i,j] stores a delta array 
emergency_rate_array = zeros(NUM_OF_RUNS,1);                    % # of objective function evaluations for the stop creteria

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replicates for NUM_OF_RUNS

mu = ceil(lambda/4);
for i = 1:NUM_OF_RUNS

    if(lambda==0)
        x0 = randn(n,1);
        a = withGP(fname,x0,sigma0,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE);
        t_array(i) = cell2mat(a(5));
    else 
        x0 = randn(n,mu);
        a = mml_GP_CSA_Niko(fname,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE);
        t_array(i) = cell2mat(a(1));
%         emergency_rate_array(i) = cell2mat(a(12));

    end
    
    
    sigma_matrix(i,:) = cell2mat(a(4));
    T_array(i) = cell2mat(a(5));
    f_x_matrix(i,:) = cell2mat(a(6));
    convergence_rate_array(i) = cell2mat(a(7));
    GP_error_matrix(i,:) = cell2mat(a(8));
    sigma_star_matrix(i,:) = cell2mat(a(9));
    success_rate_array(i) = cell2mat(a(10));
    delta_matrix(i,:) = cell2mat(a(11));


    % counter
    fprintf('# of runs %d',i)
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Save data
%     if fname == 1
%         f = @(x) (x'*x)^(1/2);
%         save('linear_over_lambda.mat','fname','f','fname','NUM_OF_RUNS','lambda_array',...
%             'LENGTH_SCALE','TRAINING_SIZE',...
%         't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
%             'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
%             'delta_matrix');
%     elseif fname == 2
%         f = @(x) (x'*x);
%         save('quadratic_over_lambda.mat','fname','f','fname','NUM_OF_RUNS','lambda_array',...
%             'LENGTH_SCALE','TRAINING_SIZE',...
%         't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
%             'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
%             'delta_matrix');
%     elseif fname == 3
%         f = @(x) (x'*x)^(3/2);
%         save('cubic_over_lambda.mat','fname','f','fname','NUM_OF_RUNS','lambda_array',...
%             'LENGTH_SCALE','TRAINING_SIZE',...
%         't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
%             'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
%             'delta_matrix');
%     elseif fname == 4
%         f = @f4;
%         save('schwefel_over_lambda.mat','fname','f','fname','NUM_OF_RUNS','lambda_array',...
%             'LENGTH_SCALE','TRAINING_SIZE',...
%         't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
%             'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
%             'delta_matrix');  
%     elseif fname == 5
%         f = @f5;
%         save('quartic_over_lambda.mat','f','fname','f','fname','NUM_OF_RUNS','lambda_array',...
%             'LENGTH_SCALE','TRAINING_SIZE',...
%         't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
%             'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
%             'delta_matrix');       
%     end
% 
% end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take median

t_max = max(t_array);
t_min = min(t_array);

% Find the index of median T 
sorted_T = sort(T_array);
med_index = find(T_array == sorted_T(ceil(length(sorted_T)/2)));
med_index = med_index(1);   % randomly pick one
T_med = T_array(med_index);
t_med = t_array(med_index);


sigma_matrix_med = sigma_matrix(med_index,:);
f_x_med = f_x_matrix(med_index,:);
convergence_med = convergence_rate_array(med_index);
GP_error_matrix_med = GP_error_matrix(med_index,:);
sigma_star_matrix_med = sigma_star_matrix(med_index,:);
success_med = success_rate_array(med_index);
delta_matrix_med = delta_matrix(med_index,:);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graph

% For f(x), sigma, sigmaStar, GP_error over objective function evaluations
figure(FIGURE_NUM);
legend('-DynamicLegend'); 
hold on;
mu = ceil(lambda/4);
if lambda==0
    d = sprintf('(1+1)-ES');
else
    d = sprintf('(%d/%d,%d)-ES',mu,mu,lambda);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.objective function [row 1]
subplot(subplot_ROW,subplot_COL,fname);
h = histogram(T_array,'Normalization','probability','DisplayName',d);hold on;
% h.BinWidth = 20;
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
        d3 =sprintf('Schwefel function');
        title(d3,'fontsize',15);
    elseif(fname == 5)
        d3 =sprintf('quartic function');
        title(d3,'fontsize',15); 
    end
end
% if(fname==1)
%     ylabel('probability','FontSize',15);%
% end
xlabel('objective function calls','FontSize',15); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
legend('-DynamicLegend'); 
legend('show');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.GP error [row 3]
subplot(subplot_ROW,subplot_COL,fname+10);
if lambda==0
    plot(t_start:1:T_med,GP_error_matrix_med(t_start:T_med),'DisplayName',d);hold on;
    d1 = sprintf('(1+1)-ES[S]');
    plot(t_start:1:T_med,exp(conv(log(GP_error_matrix_med(t_start:T_med)), kernel, 'same')),'DisplayName',d1,'LineWidth',2);hold on;
%     plot(t_start:1:T_med,smoothdata(GP_error_matrix_med(t_start:T_med),'gaussian',40),'DisplayName',d1,'LineWidth',2);hold on;
else 
%     GP_error_range1 = GP_error_matrix_med(1:t_start);
    GP_error_range2 = GP_error_matrix_med(t_start+1:t_med);
    plot(t_range2,GP_error_range2,'DisplayName',d);hold on;
    % GP smoothed
%     smoothed_GP_range1 = smoothdata(GP_error_range1,'gaussian',40);
    smoothed_GP_range2 = exp(conv(log(GP_error_range2), kernel, 'same'));
    d1 = sprintf('(%d/%d,%d)-ES[S]',mu,mu,lambda);
    plot(t_range2, smoothed_GP_range2,'DisplayName',d1,'LineWidth',2);hold on;
%     plot([t_range1 t_range2],[smoothed_GP_range1 smoothed_GP_range2],'DisplayName',d1,'LineWidth',2);hold on;
end
if(fname==1)
    ylabel('relative model error','FontSize',15);%
end
xlabel('objective function calls','FontSize',15); 
set(gca, 'YScale', 'log');

legend('-DynamicLegend'); 
legend('show');
% title('Logarithmic relative model error','FontSize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    legend('-DynamicLegend'); 
    legend('show');
end
saveas(gcf,'merged_plot_NO_emergency_v2.fig'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Plot success rate and convergence rate for sphere functions
% if(fname == 1 || fname == 2 || fname == 3||)
%     subplot_ROW = 4;
%     subplot_COL = 3;
%     % Plot convergence rate 
%     xNameSprintf = sprintf('onvergence rate');
%     xLimit = [0 1];
%     plot_pdf(convergence_rate_array,T_med,FIGURE_NUM+1,subplot_ROW,subplot_COL,fname,lambda,xNameSprintf,xLimit);
%     xNameSprintf = sprintf('success rate');
%     % Plot success rate for a good step
%     xLimit = [0 1];
%     plot_pdf(success_rate_array,T_med,FIGURE_NUM+1,subplot_ROW,subplot_COL,fname+subplot_COL*2,lambda,xNameSprintf,xLimit);
%     legend('-DynamicLegend'); 
%     legend('show');
%     saveas(gcf,'success_convergence_emergency_v1.fig'); 
% end

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

val = {t_array,sigma_matrix,T_array,f_x_matrix,convergence_rate_array,GP_error_matrix,sigma_star_matrix,success_rate_array,delta_matrix};
end