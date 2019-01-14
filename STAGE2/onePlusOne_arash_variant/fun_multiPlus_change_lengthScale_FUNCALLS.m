% Objective [update step size C1,C2,C3] add hist of objFunCalls [over DF]
% 1 Plot
%     1. histgoram of objFunCalls
%     2. convergence plot
%     3. step size 
%     4. normalized step size
%     5. success rate and evaluation rate [bar]
%     6. 4 probs [bar]
%        
% 2. Save
%     Plots
%     Data     
%
% difficulty: 
%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = fun_multiPlus_change_lengthScale_FUNCALLS(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE_array,SUCCESS_RATE,DF,C3,FIGURE_NUM,subplot_ROW,subplot_COL,str_cell_SIGMA_STAR)
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
%    SUCCESS_RATE_array  C1/(C1+C2) \approx success rate
%    DF(dampin)     C1 = 1*SUCCESS_RATE*DAMPING_FACTOR, 
%                   C2 = (1-SUCCESS_RATE)*DAMPING_FACTOR
%    C3             exp(-C3) when bad step estimated by model
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

PROB_RATE_array = LENGTH_SCALE_array;

% Four probs for each SIGMA_STAR in SIGMA_STAR_array
four_prob_final_med = zeros(length(PROB_RATE_array),4);
success_rate_final_med = zeros(length(PROB_RATE_array),1);
eval_rate_final_med = zeros(length(PROB_RATE_array),1);
T_final_med = zeros(length(PROB_RATE_array),1);

figure(FIGURE_NUM);
for q = 1:1:length(PROB_RATE_array)
    LENGTH_SCALE = PROB_RATE_array(q);
    C1 = SUCCESS_RATE*DF;
    C2 = (1-SUCCESS_RATE)*DF;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Def of variables

    % (mu/mu,lambda)-ES with GP
    t_array = zeros(NUM_OF_RUNS,1);                    % # of iterations for the stop creteria
    sigma_matrix = zeros(NUM_OF_RUNS,10000);           % store all sigma
    T_array = zeros(NUM_OF_RUNS,1);                    % # of objective function evaluations for the stop creteria
    f_x_matrix = zeros(NUM_OF_RUNS,10000);             % store all fx
    success_rate_array = zeros(NUM_OF_RUNS,1);         % success rate 
    eval_rate_array = zeros(NUM_OF_RUNS,1);            % evaluation rate 
    sigma_star_matrix = zeros(NUM_OF_RUNS,10000);      % normalized step size 
    sigma_matrix_med = zeros(1,10000);
    f_x_med = zeros(1,10000);
    four_prob_matrix = zeros(NUM_OF_RUNS,4);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Replicates for NUM_OF_RUNS

    for i = 1:NUM_OF_RUNS
        x0 = randn(n,1);
        sigma0 = 1;
        a = bestSoFar_fourProb_GP_arashVariant(fname,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,C1,C2,C3);
        t_array(i) = cell2mat(a(1));
        sigma_matrix(i,:) = cell2mat(a(4));
        T_array(i) = cell2mat(a(5));
        f_x_matrix(i,:) = cell2mat(a(6));
        success_rate_array(i) = cell2mat(a(10));
        eval_rate_array(i) = cell2mat(a(13));
        sigma_star_matrix(i,:) = cell2mat(a(9)); 
        four_prob_matrix(i,:) = cell2mat(a(12));

    end 

    % Take median run
    sorted_T = sort(T_array);
    temp_index = find(T_array == sorted_T(ceil(length(sorted_T)/2)));
    med_index = temp_index(1);
    sigma_matrix_med(:) = sigma_matrix(med_index,:);
    f_x_med(:) = f_x_matrix(med_index,:);
    t_med = t_array(med_index);
    success_rate_med = success_rate_array(med_index);
    eval_rate_med = eval_rate_array(med_index);
    sigma_star_med = sigma_star_matrix(med_index,:);
	success_rate_final_med(q) = success_rate_med;
    eval_rate_final_med(q) = eval_rate_med;
    four_prob_final_med(q,:) = four_prob_matrix(med_index,:);
    T_final_med(q) = T_array(med_index); 

    
    % Range of T for ploting
%     t_start = TRAINING_SIZE/lambda;
%     T_range_1 = (0:1:t_start).*(lambda+1)+1;
%     T_range_2 = (t_start+2:t_med)+lambda*t_start;
%     T_range = [T_range_1 T_range_2];
    T_range = 1:1:t_med;
    % counter
    fprintf('FIG %d finished \n',fname);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRPAH
    % For f(x), sigma, sigmaStar, GP_error over objective function evaluations
    
    legend('-DynamicLegend'); 
    hold on;
    
    d = sprintf('(%d/%d,%d)-ES,S=%.1f,DF=%.1f,LS=%d',mu,mu,lambda,SUCCESS_RATE,DF,LENGTH_SCALE);
    
    % Fig 1: histgragm of objFunCalls
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(subplot_ROW,subplot_COL,(0)*subplot_COL+fname);
    h = histogram(T_array);hold on; % mml with GP
    
    if(fname == 1)
        dt =sprintf('linear sphere');
        title(dt,'fontsize',15);
        h.BinWidth = 5;
    elseif(fname == 2)
        dt =sprintf('quadratic sphere');
        title(dt,'fontsize',15);
        h.BinWidth = 2;
    elseif(fname == 3)
        dt =sprintf('cubic sphere');
        title(dt,'fontsize',15);
        h.BinWidth = 2;
    elseif(fname == 4)
        dt =sprintf('Schwefel function');
        title(dt,'fontsize',15);
        h.BinWidth = 10;
    elseif(fname == 5)
        dt =sprintf('quartic function');
        title(dt,'fontsize',15); 
        h.BinWidth = 10;
    end
    xlabel('Objective function calls','FontSize',15); 
    
    
    
    % Fig 2: convergence plots [FIGURE_NUM.fig]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(subplot_ROW,subplot_COL,(1)*subplot_COL+fname);
    plot(T_range, f_x_med(1:t_med),'DisplayName',d);hold on; % mml with GP
    if(fname==1)
        ylabel(sprintf('Objective function value'),'FontSize',15);
    end
    xlabel('function calls','FontSize',15); 
    set(gca, 'YScale', 'log');
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
    legend('-DynamicLegend'); 
    legend('show');


    % Fig 3: step size [FIGURE_NUM+1.fig]
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

    % Fig 4: Normalized step size [FIGURE_NUM+2.fig]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

end
    figure(FIGURE_NUM);
    subplot(subplot_ROW,subplot_COL,(0)*subplot_COL+fname);
    legendCell = {};
    for i = 1:1:length(PROB_RATE_array)
        legendCell{i} = sprintf('S=%.1f,LS=%.1f',SUCCESS_RATE,PROB_RATE_array(i));
    end
    legend(legendCell);
    
    % Fig 5:med success rate & evaluation rate [bar][bar]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(subplot_ROW,subplot_COL,(4)*subplot_COL+fname);
    bar([success_rate_final_med, eval_rate_final_med]);hold on;
    set(gca,'xticklabel',str_cell_SIGMA_STAR);
    if(fname==1)
        ylabel(sprintf('Probabilities'),'FontSize',15);
    end
    legend({'success rate','evaluation rate'});
    xlabel('length scale (LS)','FontSize',15);  

    % Fig 6: four probs TP, TN, FP, FN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(subplot_ROW,subplot_COL,(5)*subplot_COL+fname);
    bar(four_prob_final_med./sum(four_prob_final_med,2));hold on;
    legend({'TN','FP','FN','TP'})
    set(gca,'xticklabel',str_cell_SIGMA_STAR);
    if(fname==1)
        ylabel(sprintf('Probabilities'),'FontSize',15);
    end
    xlabel('length scale (LS)','FontSize',15); 
    fig_name = sprintf('merged%d_LS.fig',lambda);
    saveas(gcf,fig_name); 

    val = {T_final_med};
% val = {t_array,sigma_matrix,T_array,f_x_matrix,convergence_rate_array,GP_error_matrix,sigma_star_matrix,success_rate_array,delta_matrix};
end