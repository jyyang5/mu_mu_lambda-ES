% Objective 
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
function val = fun_multi_over_decreaseFactor(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE,DECREASE_FACTOR_array,strategyName)
%Input:
%       f:          an index 1 for linear, 2 for quadratic and 3 for cubic 
%    name:          a number specify f  
%    NUM_OF_RUNS    # of runs to average
%    DECREASE_FACTOR_array   an array of lambda
%    TRAINING_SIZE  GP training size
%    LENGTH_SCALE   length scale factor for GP
%Return:
%    1. T_med         
%    2. convergence_med
%    3. success_med
%    4. GP_error_matrix_med
%    5. delta_matrix_med  
%    6. emergency_rate_matrix_med 
    
sigma0 = 1;
NUM_OF_ITERATIONS = 2000;
DECREASE_FATCOR_LENGTH = length(DECREASE_FACTOR_array);
% lambda = 10;
% mu = 3;
n = 10;

% close all graph
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Def of variables

% (mu/mu,lambda)-ES with GP
t_array = zeros(DECREASE_FATCOR_LENGTH,NUM_OF_RUNS,1);                    % # of iterations for the stop creteria
sigma_matrix = zeros(DECREASE_FATCOR_LENGTH,NUM_OF_RUNS,10000);           % store all sigma
T_array = zeros(DECREASE_FATCOR_LENGTH,NUM_OF_RUNS,1);                    % # of objective function evaluations for the stop creteria
f_x_matrix = zeros(DECREASE_FATCOR_LENGTH,NUM_OF_RUNS,10000);             % store all fx
convergence_rate_array = zeros(DECREASE_FATCOR_LENGTH,NUM_OF_RUNS,1);     % convergence rate 
GP_error_matrix = zeros(DECREASE_FATCOR_LENGTH,NUM_OF_RUNS,10000);        % store similar to noise-to-signal ratio
sigma_star_matrix = zeros(DECREASE_FATCOR_LENGTH,NUM_OF_RUNS,10000);      % normalized step size 
success_rate_array = zeros(DECREASE_FATCOR_LENGTH,NUM_OF_RUNS,1);         % success rate 
delta_matrix = zeros(DECREASE_FATCOR_LENGTH,NUM_OF_RUNS,10000);           % each [i,j] stores a delta array 
emergency_rate_matrix = zeros(DECREASE_FATCOR_LENGTH,NUM_OF_RUNS,1);      % count the occurance of emergencies 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterate over lambda and replicate multiple runs

for i=1:1:DECREASE_FATCOR_LENGTH
    DECREASE_FACTOR_temp = DECREASE_FACTOR_array(i);
    mu = ceil(lambda/4);
    for j = 1:NUM_OF_RUNS
        x0 = randn(n,mu);
    
        % (mu/mu,lambda)-ES with GP
        a = strategyName(fname,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,DECREASE_FACTOR_temp);
        t_array(i,j) = cell2mat(a(1));
        sigma_matrix(i,j,:) = cell2mat(a(4));
        T_array(i,j) = cell2mat(a(5));
        f_x_matrix(i,j,:) = cell2mat(a(6));
        convergence_rate_array(i,j) = cell2mat(a(7));
        GP_error_matrix(i,j,:) = cell2mat(a(8));
        sigma_star_matrix(i,j,:) = cell2mat(a(9));
        success_rate_array(i,j) = cell2mat(a(10));
        delta_matrix(i,j,:) = cell2mat(a(11));
        emergency_rate_matrix(i,j) = cell2mat(a(12));
        
    end 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save data
    if fname == 1
        f = @(x) (x'*x)^(1/2);
        save('linear_over_lambda.mat','fname','f','fname','NUM_OF_RUNS','DECREASE_FACTOR_array',...
            'DECREASE_FATCOR_LENGTH','LENGTH_SCALE','TRAINING_SIZE',...
            't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
            'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
            'delta_matrix','emergency_rate_matrix');
    elseif fname == 2
        f = @(x) (x'*x);
        save('quadratic_over_lambda.mat','fname','f','fname','NUM_OF_RUNS','DECREASE_FACTOR_array',...
            'DECREASE_FATCOR_LENGTH','LENGTH_SCALE','TRAINING_SIZE',...
            't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
            'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
            'delta_matrix','emergency_rate_matrix');
    elseif fname == 3
        f = @(x) (x'*x)^(3/2);
        save('cubic_over_lambda.mat','fname','f','fname','NUM_OF_RUNS','DECREASE_FACTOR_array',...
            'DECREASE_FATCOR_LENGTH','LENGTH_SCALE','TRAINING_SIZE',...
            't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
            'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
            'delta_matrix','emergency_rate_matrix');
    elseif fname == 4
        f = @f4;
        save('schwefel_over_lambda.mat','fname','f','fname','NUM_OF_RUNS','DECREASE_FACTOR_array',...
            'DECREASE_FATCOR_LENGTH','LENGTH_SCALE','TRAINING_SIZE',...
            't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
            'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
            'delta_matrix','emergency_rate_matrix');  
    elseif fname == 5
        f = @f5;
        save('quartic_over_lambda.mat','f','fname','f','fname','NUM_OF_RUNS','DECREASE_FACTOR_array',...
            'DECREASE_FATCOR_LENGTH','LENGTH_SCALE','TRAINING_SIZE',...
            't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
            'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
            'delta_matrix','emergency_rate_matrix');       
    end

end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take median

    t_med = median(t_array,2);
    t_max = max(t_array,2);
    t_min = min(t_array,2);
    sigma_matrix_med = squeeze(median(sigma_matrix,2));
    T_med = median(T_array,2);
    f_x_med = squeeze(median(f_x_matrix,2));
    convergence_med = median(convergence_rate_array,2);
    GP_error_matrix_med = squeeze(median(GP_error_matrix,2));
    sigma_star_matrix_med = squeeze(median(sigma_star_matrix,2));
    success_med = median(success_rate_array,2);
    delta_matrix_med = squeeze(median(delta_matrix,2));
    emergency_rate_matrix_med = median(emergency_rate_matrix,2);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graph

% For f(x), sigma, sigmaStar, GP_error over objective function evaluations
figure(fname);
legend('-DynamicLegend'); 
hold on;
for i = 1:1:DECREASE_FATCOR_LENGTH
    figure(fname);
    hold on;
    DECREASE_FACTOR_temp = DECREASE_FACTOR_array(i);
    mu = ceil(lambda/4);
    d =sprintf('DECREASE_FACTOR=%',DECREASE_FACTOR_temp);
    f_x_med_temp = f_x_med(i,:);
    sigma_med = sigma_matrix_med(i,:);
    sigma_star_med = sigma_star_matrix_med(i,:);
    GP_error_med = GP_error_matrix_med(i,:);
    
    
    % objective function
    subplot(1,4,1);
    t_start = ceil(TRAINING_SIZE/lambda);
    fx_range1 = f_x_med_temp(1:t_start);
    fx_range2 = f_x_med_temp(t_start+1:t_med(i));
    t_range1 = 1:lambda:lambda*t_start;
    t_range2 = lambda*t_start+1:lambda*t_start+length(fx_range2);
    semilogy([t_range1 t_range2], [fx_range1 fx_range2],'DisplayName',d);hold on;% mml with GP
    ylabel('objective function value f(y)','FontSize',15);%
    xlabel('objective function evaluations','FontSize',15); 
    set(gca, 'YScale', 'log');
    title('objective function value','FontSize',20);
    if(i==DECREASE_FATCOR_LENGTH)
        legend('-DynamicLegend'); 
        legend('show');
    end

    % step size 
    subplot(1,4,2);
    sigma_range1 = sigma_med(1:t_start);
    sigma_range2 = sigma_med(t_start+1:t_med(i));
    plot([t_range1 t_range2], [sigma_range1 sigma_range2],'DisplayName',d);hold on;
    ylabel('step size \sigma','FontSize',15);%
    xlabel('objective function evaluations','FontSize',15); 
    set(gca, 'YScale', 'log');
    title('step size \sigma','FontSize',20);
    if(i==DECREASE_FATCOR_LENGTH)
        legend('-DynamicLegend'); 
        legend('show');
    end
    
    % normalized step size
    subplot(1,4,3);
    sigma_star_range1 = sigma_star_med(1:t_start);
    sigma_star_range2 = sigma_star_med(t_start+1:t_med(i));
    plot([t_range1 t_range2], [sigma_star_range1 sigma_star_range2],'DisplayName',d);hold on;
    ylabel('normalized step size \sigma*','FontSize',15);%
    xlabel('objective function evaluations','FontSize',15); 
    set(gca, 'YScale', 'log');
    title('normalized step size \sigma*','FontSize',20);
    if(i==DECREASE_FATCOR_LENGTH)
        legend('-DynamicLegend'); 
        legend('show');
    end
    
    % GP error
    subplot(1,4,4);
    GP_error_range1 = GP_error_med(1:t_start);
    GP_error_range2 = GP_error_med(t_start+1:t_med(i));
    plot([t_range1 t_range2], [GP_error_range1 GP_error_range2],'DisplayName',d);hold on;
    ylabel('relative model error','FontSize',15);%
    xlabel('objective function evaluations','FontSize',15); 
    set(gca, 'YScale', 'log');
    title('Logarithmic relative model error','FontSize',20);
    if(i==DECREASE_FATCOR_LENGTH)
        legend('-DynamicLegend'); 
        legend('show');
    end
    
    
    
end
    hold off;
    % svae plot    
    if fname == 1
        saveas(gcf,'linear_sphere_funCall.fig');
    elseif fname == 2
        saveas(gcf,'quadratic_sphere_funCall.fig');
    elseif fname == 3
        saveas(gcf,'cubic_sphere_funCall.fig');
    elseif fname == 4 
        saveas(gcf,'schwefel_function_funCall.fig');
    elseif fname == 5
        saveas(gcf,'quartic_function_funCall.fig');
    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each lambda plot pmf fitness gain/iteration hitogram
% assume 3 rows 
    numOfCol = ceil(DECREASE_FATCOR_LENGTH/3);
    figure(10+i+fname);
    for i = 1:1:DECREASE_FATCOR_LENGTH     
        subplot(3,numOfCol,i);
        DECREASE_FACTOR_temp = DECREASE_FACTOR_array(i);
        mu = ceil(lambda/4);
    
        delta_array = delta_matrix_med(i,1:t_med(i));
        histogram(delta_array(1:t_med(i)),'Normalization','probability');
        xlabel('number of iterations','fontsize',15);
        ylabel('prob','fontsize',15);
        d =sprintf('Normalized fitGain pdf DECREASE_FACTOR=%d',DECREASE_FACTOR_temp);
        title(d,'fontsize',20);
   
    end    
 
    if fname == 1
        d =sprintf('hist_linear_%d.fig',lambda);
    elseif fname == 2
        d =sprintf('hist_quadratic_%d.fig',lambda);
    elseif fname == 3
        d =sprintf('hist_cubic_%d.fig',lambda);
    elseif fname == 4 
        d =sprintf('hist_schwefel_%d.fig',lambda);
    elseif fname == 5
        d =sprintf('hist_quartic_%d.fig',lambda);
        saveas(gcf,'hist_quartic_function.fig');
    end
    saveas(gcf,d);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot by taking the mean of relative model error in median run
    mean_GP_error_array = zeros(DECREASE_FATCOR_LENGTH,1);
    for i = 1:1:DECREASE_FATCOR_LENGTH   
        mean_GP_error_array(i) = mean(GP_error_matrix_med(i,1:t_med(i)));
    end
    figure(50+fname)
    plot(DECREASE_FACTOR_array,mean_GP_error_array);
    ylabel('mean of relative model error in median run','FontSize',15);%
    xlabel('DECREASE_FACTOR','FontSize',15); 
    set(gca, 'YScale', 'log');
    title('Relative model error (mean)','FontSize',20);
    
    if fname == 1
        d3 = sprintf('linear_modelError_lambda.fig');
    elseif fname == 2
        d3 = sprintf('quadratic_modelError_lambda.fig');
    elseif fname == 3
        d3 = sprintf('cubic_modelError_lambda.fig');
    elseif fname == 4 
        d3 = sprintf('schwefel_modelError_lambda.fig');
    elseif fname == 5
        d3 = sprintf('quartic_modelError_lambda.fig');
    end
    saveas(gcf,d3);
    
    
    val = {T_med,convergence_med,success_med,GP_error_matrix_med,delta_matrix_med,emergency_rate_matrix_med};

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val] = f4(x)
    val = 0;
    for i = 1:length(x)
        temp = 0;
        for j = 1:i
            temp = temp + x(j);
        end
        val = val + temp^2;          
    end
end

function [val] = f5(x)
    val = 0;
    beta = 1;
    for i = 1:length(x)-1
    val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
    end
end