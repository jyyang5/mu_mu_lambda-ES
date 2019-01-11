% Objective [update step size using different success rate]
% 1 Plot
%     1. convergence plot
%     2. step size 
%     3. normalized step size
%     4. success rate [hist]
%     5. success rate [mid]
%        
% 2. Save
%     Plots
%     Data     
%
% difficulty: 
%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = fun_multi_run_change_success_rate_FUNCALLS(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE,PROB_RATE_array,FIGURE_NUM,subplot_ROW,subplot_COL,str_cell_SIGMA_STAR)
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
%    SUCCESS_RATE   success rate used to update step size
%    FIGURE_NUM     the first fig. number
%    subplot_ROW    # of rows in each plot [here # of lambda used]
%    subplot_COL    # of cols in each plot [here # of test functions]
%    str_cell_SIGMA_STAR   for individual name in x-axes
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
NUM_OF_ITERATIONS = 2500;        
% SIGMA_STAR_array = SIGMA_STAR_array;
% LEN_SIGMA_STAR = length(SIGMA_STAR_array);

% lambda = 10;
% mu = 3;
n = 10;
mu = ceil(lambda/4);
% For compact subplots 
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.03], [0.05 0.03], [0.05 0.03]);
if ~make_it_tight,  clear subplot;  end

% GP smooth 
window_length = 40;
kernel = exp(-(-3*window_length:3*window_length).^2/window_length^2/2);
kernel = kernel/sum(kernel);        % Normalized    


% Four probs for each SIGMA_STAR in SIGMA_STAR_array
four_prob_final_med = zeros(length(PROB_RATE_array),4);
success_rate_final_med = zeros(length(PROB_RATE_array),1);
T_final_med = zeros(length(PROB_RATE_array),1);
T_array_all = zeros(length(PROB_RATE_array),NUM_OF_RUNS);

figure(FIGURE_NUM);
for q = 1:1:length(PROB_RATE_array)
    PROB_RATE = PROB_RATE_array(q);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Def of variables

    % (mu/mu,lambda)-ES with GP
    t_array = zeros(NUM_OF_RUNS,1);                    % # of iterations for the stop creteria
    sigma_matrix = zeros(NUM_OF_RUNS,10000);           % store all sigma
    T_array = zeros(NUM_OF_RUNS,1);                    % # of objective function evaluations for the stop creteria
    f_x_matrix = zeros(NUM_OF_RUNS,10000);             % store all fx
    success_rate_array = zeros(NUM_OF_RUNS,1);         % success rate 
    % convergence_rate_array = zeros(NUM_OF_RUNS,1);     % convergence rate 
    % GP_error_matrix = zeros(NUM_OF_RUNS,10000);        % store similar to noise-to-signal ratio
    sigma_star_matrix = zeros(NUM_OF_RUNS,10000);      % normalized step size 
    % success_rate_array = zeros(NUM_OF_RUNS,1);         % success rate 
    % delta_matrix = zeros(NUM_OF_RUNS,10000);           % each [i,j] stores a delta array 
    % emergency_rate_array = zeros(NUM_OF_RUNS,1);                    % # of objective function evaluations for the stop creteria

    sigma_matrix_med = zeros(1,10000);
    f_x_med = zeros(1,10000);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Replicates for NUM_OF_RUNS

    for i = 1:NUM_OF_RUNS
        x0 = randn(n,1);
        sigma0 = 1;
        a = bestOfTwo_GP_change_success_rate(fname,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,PROB_RATE);
        t_array(i) = cell2mat(a(1));
        sigma_matrix(i,:) = cell2mat(a(4));
        T_array(i) = cell2mat(a(5));
        T_array_all(q,i) = cell2mat(a(5));
        f_x_matrix(i,:) = cell2mat(a(6));
        success_rate_array(i) = cell2mat(a(10));
        sigma_star_matrix(i,:) = cell2mat(a(9)); 
%         four_prob_matrix(i,:) = cell2mat(a(12));
    %     b = mml_noGP_SIGMA_STAR(fname,x0,SIGMA_STAR,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE);
    %     t_array_noGP(i) = cell2mat(b(1));  
    %     sigma_matrix_noGP(i,:) = cell2mat(b(4));
    %     T_array_noGP(i) = cell2mat(b(5));
    %     f_x_matrix_noGP(i,:) = cell2mat(b(6));
    %     success_rate_array_noGP(i) = cell2mat(b(10));

    end 

    % Take median run
    sorted_T = sort(T_array);
    temp_index = find(T_array == sorted_T(ceil(length(sorted_T)/2)));
    med_index = temp_index(1);
    sigma_matrix_med(:) = sigma_matrix(med_index,:);
    f_x_med(:) = f_x_matrix(med_index,:);
    t_med = t_array(med_index);
    success_rate_med = success_rate_array(med_index);
    sigma_star_med = sigma_star_matrix(med_index,:);
	success_rate_final_med(q) = success_rate_med;
%     four_prob_final_med(q,:) = four_prob_matrix(med_index,:);
    T_final_med(q) = T_array(med_index); 

    % sorted_T_noGP = sort(T_array_noGP);
    % temp_index_noGP = find(T_array_noGP == sorted_T_noGP(ceil(length(sorted_T_noGP)/2)));
    % med_index_noGP = temp_index_noGP(1);
    % sigma_matrix_med_noGP(:) = sigma_matrix_noGP(med_index_noGP,:);
    % f_x_med_noGP(:) = f_x_matrix_noGP(med_index_noGP,:);
    % t_med_noGP = t_array_noGP(med_index_noGP);
    % success_rate_med_noGP = success_rate_array_noGP(med_index_noGP);
    
    % Range of T for ploting
    t_start = TRAINING_SIZE/lambda;
    T_range_1 = (0:1:t_start).*(lambda+1)+1;
    T_range_2 = (t_start+2:t_med)+lambda*t_start;
    T_range = [T_range_1 T_range_2];

    % counter
    fprintf('FIG %d finished \n',fname);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % graph

    % For f(x), sigma, sigmaStar, GP_error over objective function evaluations
    
    legend('-DynamicLegend'); 
    hold on;
    
    d = sprintf('(%d/%d,%d)-ES,S=%.3f',mu,mu,lambda,PROB_RATE);
    % d_noGP = sprintf('(%d/%d,%d)-ES',mu,mu,lambda);


    % Fig 0:objFun calls hist.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(subplot_ROW,subplot_COL,(0)*subplot_COL+fname);
    histogram(T_array,'DisplayName',d);hold on; % mml with GP
    if(fname == 1)
        dt =sprintf('linear sphere');
        title(dt,'fontsize',15);
    elseif(fname == 2)
        dt =sprintf('quadratic sphere');
        title(dt,'fontsize',15);
    elseif(fname == 3)
        dt =sprintf('cubic sphere');
        title(dt,'fontsize',15);
    elseif(fname == 4)
        dt =sprintf('Schwefel function');
        title(dt,'fontsize',15);
    elseif(fname == 5)
        dt =sprintf('quartic function');
        title(dt,'fontsize',15); 
    end
    xlabel('Objective function calls','FontSize',15); 
    % plot(1:t_med_noGP, f_x_med_noGP(i,1:t_med_noGP),'DisplayName',d_noGP);hold on; % mml with GP
        
    

    % Fig 1: convergence plots [FIGURE_NUM.fig]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(subplot_ROW,subplot_COL,(1)*subplot_COL+fname);
    plot(T_range, f_x_med(1:t_med),'DisplayName',d);hold on; % mml with GP
    % plot(1:t_med_noGP, f_x_med_noGP(i,1:t_med_noGP),'DisplayName',d_noGP);hold on; % mml with GP
    if(fname==1)
        ylabel(sprintf('Objective function value'),'FontSize',15);
    end
        
    xlabel('function calls','FontSize',15); 
    set(gca, 'YScale', 'log');

    legend('-DynamicLegend'); 
    legend('show');


    % Fig 2: step size [FIGURE_NUM+1.fig]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subplot(subplot_ROW,subplot_COL,(2)*subplot_COL+fname);
    plot(T_range, sigma_matrix_med(1:t_med),'DisplayName',d);hold on; % mml with GP
    % plot(1:t_med_noGP, f_x_med_noGP(i,1:t_med_noGP),'DisplayName',d_noGP);hold on; % mml with GP
    if(fname==1)
        ylabel(sprintf('Step size'),'FontSize',15);
    end
    xlabel('function calls','FontSize',15);
    set(gca, 'YScale', 'log');

    legend('-DynamicLegend'); 
    legend('show');


    % Fig 3: Normalized step size [FIGURE_NUM+2.fig]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if fname<=3
        subplot(subplot_ROW,subplot_COL,(3)*subplot_COL+fname);
        plot(T_range, sigma_star_med(1:t_med),'DisplayName',d);hold on; % mml with GP
        % plot(1:t_med_noGP, f_x_med_noGP(i,1:t_med_noGP),'DisplayName',d_noGP);hold on; % mml with GP
        if(fname==1)
            ylabel(sprintf('Normalized step size'),'FontSize',15);
        end
        xlabel('function calls','FontSize',15);
        set(gca, 'YScale', 'log');
        legend('-DynamicLegend'); 
        legend('show');
%     end


    % % Fig 4:Success rate
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % subplot(subplot_ROW,subplot_COL,(3-1)*subplot_COL+fname);
    % histogram(success_rate_array,'DisplayName',d);hold on;
    % % plot(1:t_med_noGP, f_x_med_noGP(i,1:t_med_noGP),'DisplayName',d_noGP);hold on; % mml with GP
    % if(fname==1)
    %     ylabel(sprintf('Success rate'),'FontSize',15);
    % end
    % xlabel('number of iterations','FontSize',15); 
    % % set(gca, 'YScale', 'log');
    % legend('-DynamicLegend'); 
    % legend('show');
    % 
    % saveas(gcf,'hist_success_rate.fig'); 

    
    
    

% legend('-DynamicLegend'); 
% legend('show');
end
    legend('-DynamicLegend'); 
    legend('show');
    
    figure(FIGURE_NUM);
    subplot(subplot_ROW,subplot_COL,(0)*subplot_COL+fname);   
    legend({'S=0.7','S=0.8','S=0.9'});
    % Fig 4:med success rate[bar]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(subplot_ROW,subplot_COL,(4)*subplot_COL+fname);
    bar(success_rate_final_med);hold on;
    set(gca,'xticklabel',str_cell_SIGMA_STAR);
    if(fname==1)
        ylabel(sprintf('Success rate'),'FontSize',15);
    end
    xlabel('S','FontSize',15); 

   
%     % Fig 5:med success rate[bar]
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     subplot(subplot_ROW,subplot_COL,(5-1)*subplot_COL+fname);
%     bar(four_prob_final_med./sum(four_prob_final_med,2));hold on;
%     legend({'TN','FP','FN','TP'})
%     set(gca,'xticklabel',str_cell_SIGMA_STAR);
%     if(fname==1)
%         ylabel(sprintf('Probabilities'),'FontSize',15);
%     end
%     xlabel('S','FontSize',15); 
    fig_name = sprintf('bestOfTwo_merged_%d.fig',lambda);
    saveas(gcf,fig_name); 

    val = {T_final_med};
% val = {t_array,sigma_matrix,T_array,f_x_matrix,convergence_rate_array,GP_error_matrix,sigma_star_matrix,success_rate_array,delta_matrix};
end