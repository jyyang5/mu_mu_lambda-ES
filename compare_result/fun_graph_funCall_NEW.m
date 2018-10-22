% plot 
%     1. step size, objective function value, normalized step size over
%       funCalls
%     2. histgogram of objective function calls, convergence, success rate
%
% compare mml-ES with (1+1)-ES with model assistance
% save objective function calls data & plot into files
%
% difficulty: 
%            first few trainijng iterations
%            histogram does not show data properly for sigmaStar and model
%            error
%            save file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = fun_graph_funCall_NEW(f,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE)
%Input:
%       f:          objective function 
%    name:          a number specify f  
%    NUM_OF_RUNS    # of runs to average
%Return:
%    iteration number for [mmlWithGP,mmlNoGP,1+1WithGP,1+1NoGP]    
sigma0 = 1;
NUM_OF_ITERATIONS = 10000;
% lambda = 10;
% mu = 3;
n = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Def of variables

% (mu/mu,lambda)-ES with GP
t_array = zeros(NUM_OF_RUNS,1);                    % # of iterations for the stop creteria
sigma_matrix = zeros(NUM_OF_RUNS,10000);           % store all sigma
T_array = zeros(NUM_OF_RUNS,1);                    % # of objective function evaluations for the stop creteria
f_x_matrix = zeros(NUM_OF_RUNS,10000);             % store all fx
convergence_rate_array = zeros(NUM_OF_RUNS,1);     % convergence rate 
GP_error_matrix = zeros(NUM_OF_RUNS,10000);           % store similar to noise-to-signal ratio
sigma_star_matrix = zeros(NUM_OF_RUNS,10000);      % normalized step size 
success_rate_array = zeros(NUM_OF_RUNS,1);         % success rate 

% (1+1)-ES with GP
t_array1 = zeros(NUM_OF_RUNS,1);                    % # of iterations for the stop creteria
sigma_matrix1 = zeros(NUM_OF_RUNS,10000);           % store all sigma
T_array1 = zeros(NUM_OF_RUNS,1);                    % # of objective function evaluations for the stop creteria
f_x_matrix1 = zeros(NUM_OF_RUNS,10000);             % store all fx
convergence_rate_array1 = zeros(NUM_OF_RUNS,1);     % convergence rate 
GP_error_matrix1 = zeros(NUM_OF_RUNS,10000);        % store similar to noise-to-signal ratio
sigma_star_matrix1 = zeros(NUM_OF_RUNS,10000);      % normalized step size 
success_rate_array1 = zeros(NUM_OF_RUNS,1);         % success rate 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replicates and take the median 

for i = 1:NUM_OF_RUNS
    
    x0 = randn(n,mu);
    
    % (mu/mu,lambda)-ES with GP
    a = mml_GP_CSA(f,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE);
    t_array(i) = cell2mat(a(1));
    sigma_matrix(i,:) = cell2mat(a(4));
    T_array(i) = cell2mat(a(5));
    f_x_matrix(i,:) = cell2mat(a(6));
    convergence_rate_array(i) = cell2mat(a(7));
    GP_error_matrix(i,:) = cell2mat(a(8));
    sigma_star_matrix(i,:) = cell2mat(a(9));
    success_rate_array(i) = cell2mat(a(10));
    

    x0 = randn(n,1);
    
    % (1+1)-ES with GP
    TRAINING_SIZE = 40;
    b = withGP(f,x0,sigma0,NUM_OF_ITERATIONS,TRAINING_SIZE);
    t_array1(i) = cell2mat(b(1));
    sigma_matrix1(i,:) = cell2mat(b(4));
    T_array1(i) = cell2mat(b(5));
    f_x_matrix1(i,:) = cell2mat(b(6));
    convergence_rate_array1(i) = cell2mat(b(7));
    GP_error_matrix1(i,:) = cell2mat(b(8));
    sigma_star_matrix1(i,:) = cell2mat(b(9));
    success_rate_array1(i) = cell2mat(b(10));
    
    % counter
    di = sprintf('run %d',i);
    disp(di);
 
end
    t_med = median(t_array);
    t_max = max(t_array);
    t_min = min(t_array);
    sigma_med = median(sigma_matrix);
    T_med = median(T_array);
    f_x_med = median(f_x_matrix);
    convergence_med = median(convergence_rate_array);
    GP_error_med = median(GP_error_matrix);
    sigma_star_med = median(sigma_star_matrix);
    success_med = median(success_rate_array);
    
    
    t_med1 = median(t_array1);
    t_max1 = max(t_array1);
    t_min1 = min(t_array1);
    sigma_med1 = median(sigma_matrix1);
    T_med1 = median(T_array1);
    f_x_med1 = median(f_x_matrix1);
    convergence_med1 = median(convergence_rate_array1);
    GP_error_med1 = median(GP_error_matrix1);
    sigma_star_med1 = median(sigma_star_matrix1);
    success_med1 = median(success_rate_array1);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graph

% For f(x), sigma, sigmaStar, GP_error over objective function evaluations
figure(name);
      
    subplot(1,4,1);
    t_start = ceil(TRAINING_SIZE/lambda);
    fx_range1 = f_x_med(1:t_start);
    fx_range2 = f_x_med(t_start+1:t_med);
    t_range1 = 1:lambda:lambda*t_start;
    t_range2 = lambda*t_start+1:lambda*t_start+length(fx_range2);
    semilogy([t_range1 t_range2], [fx_range1 fx_range2]);hold on;% mml with GP
    % plot([lambda lambda+1:T-lambda-1],[fcentroid_array(1) fcentroid_array(2:t)]);hold on;
    plot(1:1:t_med1,f_x_med1(1:t_med1));
    ylabel('objective function value f(y)','FontSize',15);%
    xlabel('objective function evaluations','FontSize',15); 
    set(gca, 'YScale', 'log');
    legend({'mml-ES','(1+1)-ES'},'FontSize',10); %
    title('objective function value','FontSize',20);
    hold off;


    subplot(1,4,2);
    sigma_range1 = sigma_med(1:t_start);
    sigma_range2 = sigma_med(t_start+1:t_med);
    plot([t_range1 t_range2], [sigma_range1 sigma_range2]);hold on;
    % plot(1:1:t,sigma_star_array(1:t));hold on;
    plot(1:1:t_med1,sigma_med1(1:t_med1));
    ylabel('step size \sigma','FontSize',15);%
    xlabel('objective function evaluations','FontSize',15); 
    set(gca, 'YScale', 'log');
    legend({'mml-ES','(1+1)-ES'},'FontSize',10); %
    title('step size \sigma','FontSize',20);
    hold off;

    subplot(1,4,3);
    sigma_star_range1 = sigma_star_med(1:t_start);
    sigma_star_range2 = sigma_star_med(t_start+1:t_med);
    plot([t_range1 t_range2], [sigma_star_range1 sigma_star_range2]);hold on;
    % plot(1:1:t,sigma_star_array(1:t));hold on;
    plot(1:1:t_med1,sigma_star_med1(1:t_med1));
    ylabel('normalized step size \sigma*','FontSize',15);%
    xlabel('objective function evaluations','FontSize',15); 
    set(gca, 'YScale', 'log');
    legend({'mml-ES','(1+1)-ES'},'FontSize',10); %
    title('normalized step size \sigma*','FontSize',20);
    hold off;

    subplot(1,4,4);
    GP_error_range1 = GP_error_med(1:t_start);
    GP_error_range2 = GP_error_med(t_start+1:t_med);
    plot([t_range1 t_range2], [GP_error_range1 GP_error_range2]);hold on;
    % plot(1:t,GP_error(1:t));hold on;
    plot(1:t_med1,GP_error_med1(1:t_med1));
    ylabel('relative model error','FontSize',15);%
    xlabel('objective function evaluations','FontSize',15); 
    set(gca, 'YScale', 'log');
    legend({'mml-ES','(1+1)-ES'},'FontSize',10); %
    title('relative model error','FontSize',20);
    hold off;
    if name == 6
        saveas(gcf,'linear_sphere_funCall.fig');
    elseif name == 7
        saveas(gcf,'quadratic_sphere_funCall.fig');
    elseif name == 8
        saveas(gcf,'cubic_sphere_funCall.fig');
    elseif name == 9 
        saveas(gcf,'schwefel_function_funCall.fig');
    elseif name == 10
        saveas(gcf,'quartic_function_funCall.fig');
    end
    

% For percentiles convergence rate, success rate, # objective function evaluation 

figure(10+name)
    
    % objective function evaluations
    subplot(1,3,1);      
    histogram(T_array,'Normalization','probability');hold on;
    histogram(T_array1,'Normalization','probability');
    legend({'mml-ES','(1+1)-ES'},'FontSize',10); 
    title('Objective function evaluations','FontSize',20);

    % convergence rate
    subplot(1,3,2);
    histogram(convergence_rate_array,'Normalization','probability');hold on;
    histogram(convergence_rate_array1,'Normalization','probability');
    legend({'mml-ES','(1+1)-ES'},'FontSize',10); 
    title('Convergence rate','FontSize',20);

    % success rate
    subplot(1,3,3);
    histogram(success_rate_array,'Normalization','probability');hold on;
    histogram(success_rate_array1,'Normalization','probability');
    legend({'mml-ES','(1+1)-ES'},'FontSize',10); 
    title('Success rate','FontSize',20);
    hold off;
    
    if name == 6
        saveas(gcf,'hist_linear_sphere.fig');
    elseif name == 7
        saveas(gcf,'hist_quadratic_sphere.fig');
    elseif name == 8
        saveas(gcf,'hist_cubic_sphere.fig');
    elseif name == 9 
        saveas(gcf,'hist_schwefel_function.fig');
    elseif name == 10
        saveas(gcf,'hist_quartic_function.fig');
    end
    
    
%     % normalized step size 
%     subplot(2,3,4);
%     histogram(sigma_star_matrix,'Normalization','probability');hold on;
%     histogram(sigma_star_matrix1,'Normalization','probability');
%     legend({'mml-ES','(1+1)-ES'},'FontSize',10); 
%     title('Normalized step size \sigma*','FontSize',20);
%     set(gca, 'XScale', 'log');
%     hold off;
%     
%     % relative GP error
%     subplot(2,3,5);
%     histogram(GP_error_matrix,'Normalization','probability');hold on;
%     histogram(GP_error_matrix1,'Normalization','probability');
%     legend({'mml-ES','(1+1)-ES'},'FontSize',10);
%     title('Relative model error','FontSize',20);
%     set(gca, 'XScale', 'log');
%     hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save data
    
    if name == 6
        save('linear_sphere.mat','f','NUM_OF_RUNS','mu','lambda','TRAINING_SIZE',...
        't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
        'GP_error_matrix','sigma_star_matrix','success_rate_array','t_array1',...
        'sigma_matrix1','T_array1','f_x_matrix1','convergence_rate_array1',...
        'GP_error_matrix1','sigma_star_matrix1','success_rate_array1'); 
    elseif name == 7
        save('quadratic_sphere.mat','f','NUM_OF_RUNS','mu','lambda','TRAINING_SIZE',...
        't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
        'GP_error_matrix','sigma_star_matrix','success_rate_array','t_array1',...
        'sigma_matrix1','T_array1','f_x_matrix1','convergence_rate_array1',...
        'GP_error_matrix1','sigma_star_matrix1','success_rate_array1'); 
    elseif name == 8
        save('cubic_sphere.mat','f','NUM_OF_RUNS','mu','lambda','TRAINING_SIZE',...
        't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
        'GP_error_matrix','sigma_star_matrix','success_rate_array','t_array1',...
        'sigma_matrix1','T_array1','f_x_matrix1','convergence_rate_array1',...
        'GP_error_matrix1','sigma_star_matrix1','success_rate_array1');        
    elseif name == 9
        save('schwefel.mat','f','NUM_OF_RUNS','mu','lambda','TRAINING_SIZE',...
        't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
        'GP_error_matrix','sigma_star_matrix','success_rate_array','t_array1',...
        'sigma_matrix1','T_array1','f_x_matrix1','convergence_rate_array1',...
        'GP_error_matrix1','sigma_star_matrix1','success_rate_array1');        
    elseif name == 10
        save('quartic.mat','f','NUM_OF_RUNS','mu','lambda','TRAINING_SIZE',...
        't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
        'GP_error_matrix','sigma_star_matrix','success_rate_array','t_array1',...
        'sigma_matrix1','T_array1','f_x_matrix1','convergence_rate_array1',...
        'GP_error_matrix1','sigma_star_matrix1','success_rate_array1');        
    end
    
    
    val = {T_med,T_med1,convergence_med,convergence_med1,success_med,success_med1};

    
end