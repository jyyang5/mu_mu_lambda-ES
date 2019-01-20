% Objective [update step size C1,C2,C3] add hist of objFunCalls [over lambda]
% 1 Plot speed-ups of over (1+1)-ES
%      GP-(1+1)-ES
%      GP-(3/3,10)-ES
%      GP-(5/5,20)-ES
%      GP-(10/10,40)-ES
%    over 1.1 test functions: sphere, quartic, ellipsoid, and schewefel
%    over 1.2 dimension of data: [4,8,10,16]
% 2. Save
%      Data 
%      plots
%
% difficulty: 
%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_multi_over_dim(n,NUM_OF_RUNS,f6_range, f7_range, f8_range,TRAINING_SIZE,LS_onePlusOne,LS_mml,NUM_OF_ITERATIONS,FIGURE_NUM,subplot_ROW,subplot_COL,fig_row_index,C1,C2,C3)
%Input:
%    n:                     dim of data
%    NUM_OF_RUNS:           # of replicates 
%    f6_range:              range of exponents used in sphere functions 
%    f7_range:              range of \alpha used in quartic functions 
%    f8_range:              range of \beta used in ellipsoids functions 

%    TRAINING_SIZE:         GP training size
%    LS_onePlusOne:         length scale factor for GP [GP-(1+1)-ES]
%    LS_mml:                length scale factor for GP [mml-(1+1)-ES]
%    NUM_OF_ITERATIONS:     max number of objective function calls   
%    FIGURE_NUM:            the first fig. number
%    subplot_ROW:           # of rows in each plot [here # of n]
%    subplot_COL:           # of cols in each plot [here # of test functions]
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


% For compact subplots 
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.03], [0.05 0.03], [0.05 0.03]);
if ~make_it_tight,  clear subplot;  end

% GP smooth 
window_length = 40;
kernel = exp(-(-3*window_length:3*window_length).^2/window_length^2/2);
kernel = kernel/sum(kernel);        % Normalized    

NUM_OF_STRATEGIES = 5;

T_med_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range));
f_x_med_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range),50000);
sigma_med_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range),50000);
sigma_star__med_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range),50000);
success_med_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range),50000);
four_prob_med_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range),4);
eval_rate_med_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range));

T_med_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range));
f_x_med_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range),50000);
sigma_med_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range),50000);
sigma_star__med_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range),50000);
success_med_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range),50000);
four_prob_med_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range),4);
eval_rate_med_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range));

T_med_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range));
f_x_med_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range),50000);
sigma_med_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range),50000);
sigma_star__med_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range),50000);
success_med_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range),50000);
four_prob_med_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range),4);
eval_rate_med_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range));

f9_range = [1];
T_med_f9 = zeros(NUM_OF_STRATEGIES,length(f9_range));
f_x_med_f9 = zeros(NUM_OF_STRATEGIES,length(f9_range),50000);
sigma_med_f9 = zeros(NUM_OF_STRATEGIES,length(f9_range),50000);
sigma_star__med_f9 = zeros(NUM_OF_STRATEGIES,length(f9_range),50000);
success_med_f9 = zeros(NUM_OF_STRATEGIES,length(f9_range),50000);
four_prob_med_f9 = zeros(NUM_OF_STRATEGIES,length(f9_range),4);
eval_rate_med_f9 = zeros(NUM_OF_STRATEGIES,length(f9_range));

sigma0 = 1;


for fname = 6:9
    
    if fname == 6
        para_arrya = f6_range;
    elseif fname == 7
        para_arrya = f7_range;
    elseif fname ==8
        para_arrya = f8_range;
    elseif fname ==9
        para_arrya = f9_range;
    end
    for para_i = 1:length(para_arrya)
        para = para_arrya(para_i);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Replicates for NUM_OF_RUNS
        sigma_matrix = zeros(NUM_OF_STRATEGIES,NUM_OF_RUNS,50000);           % store all sigma
        T_array = zeros(NUM_OF_STRATEGIES,NUM_OF_RUNS,1);                    % # of objective function evaluations for the stop creteria
        f_x_matrix = zeros(NUM_OF_STRATEGIES,NUM_OF_RUNS,50000);             % store all fx
        success_rate_array = zeros(NUM_OF_STRATEGIES,NUM_OF_RUNS,1);         % success rate 
        eval_rate_array = zeros(NUM_OF_STRATEGIES,NUM_OF_RUNS,1);            % evaluation rate 
        sigma_star_matrix = zeros(NUM_OF_STRATEGIES,NUM_OF_RUNS,50000);      % normalized step size 
        four_prob_matrix = zeros(NUM_OF_STRATEGIES,NUM_OF_RUNS,4);           % four probs.
        
        for i = 1:NUM_OF_RUNS
            x0 = randn(n,1);
            a = onePlusOne(fname,para,x0,sigma0,NUM_OF_ITERATIONS);
            b = withGP(fname,para,x0,sigma0,NUM_OF_ITERATIONS,TRAINING_SIZE,LS_onePlusOne);
            x0 = randn(n,1);
            lambda = 10;
            c1 = bestSoFar_arashVariant(fname,para,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LS_mml,C1,C2,C3);
            lambda = 20;
            c2 = bestSoFar_arashVariant(fname,para,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LS_mml,C1,C2,C3);
            lambda = 40;
            c3 = bestSoFar_arashVariant(fname,para,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LS_mml,C1,C2,C3);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Assignment
            sigma_matrix(1,i,:) = cell2mat(a(4));
            T_array(1,i) = cell2mat(a(5));
            f_x_matrix(1,i,:) = cell2mat(a(6));
            success_rate_array(1,i) = cell2mat(a(10));
            sigma_star_matrix(1,i,:) = cell2mat(a(9)); 
            
            sigma_matrix(2,i,:) = cell2mat(b(4));
            T_array(2,i) = cell2mat(b(5));
            f_x_matrix(2,i,:) = cell2mat(b(6));
            success_rate_array(2,i) = cell2mat(b(10));
            sigma_star_matrix(2,i,:) = cell2mat(b(9)); 
            four_prob_matrix(2,i,:) = cell2mat(b(12));
            eval_rate_array(2,i) = cell2mat(b(13));
            
            
            sigma_matrix(3,i,:) = cell2mat(c1(4));
            T_array(3,i) = cell2mat(c1(5));
            f_x_matrix(3,i,:) = cell2mat(c1(6));
            success_rate_array(3,i) = cell2mat(c1(10));
            sigma_star_matrix(3,i,:) = cell2mat(c1(9));
            four_prob_matrix(3,i,:) = cell2mat(c1(12));
            eval_rate_array(3,i) = cell2mat(c1(13));
            
            
            sigma_matrix(4,i,:) = cell2mat(c2(4));
            T_array(4,i) = cell2mat(c2(5));
            f_x_matrix(4,i,:) = cell2mat(c2(6));
            success_rate_array(4,i) = cell2mat(c2(10));
            sigma_star_matrix(4,i,:) = cell2mat(c2(9));
            four_prob_matrix(4,i,:) = cell2mat(c2(12));
            eval_rate_array(4,i) = cell2mat(c2(13));
            
            sigma_matrix(5,i,:) = cell2mat(c3(4));
            T_array(5,i) = cell2mat(c3(5));
            f_x_matrix(5,i,:) = cell2mat(c3(6));
            success_rate_array(5,i) = cell2mat(c3(10));
            sigma_star_matrix(5,i,:) = cell2mat(c3(9));
            four_prob_matrix(5,i,:) = cell2mat(c3(12));
            eval_rate_array(5,i) = cell2mat(c3(13));
        end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Take median run
        for strategy_i = 1:1:5
            temp_T_array = squeeze(T_array(strategy_i,:));
            sorted_T = sort(temp_T_array);
            temp_index = find(temp_T_array == sorted_T(ceil(length(sorted_T)/2)));
            med_index = temp_index(1);
            if fname == 6
                T_med_f6(strategy_i,para_i) = T_array(strategy_i,med_index);
                f_x_med_f6(strategy_i,para_i,:) = f_x_matrix(strategy_i,med_index,:);
                sigma_med_f6(strategy_i,para_i,:) = sigma_matrix(strategy_i,med_index,:);
                sigma_star__med_f6(strategy_i,para_i,:) = sigma_star_matrix(strategy_i,med_index,:);
                success_med_f6(strategy_i,para_i) = success_rate_array(strategy_i,med_index);
                if strategy_i ~= 1
                    four_prob_med_f6(strategy_i,para_i,:) = four_prob_matrix(strategy_i,med_index,:);
                    eval_rate_med_f6(strategy_i,para_i) = eval_rate_array(strategy_i,med_index);
                end
            elseif fname == 7
                T_med_f7(strategy_i,para_i) = T_array(strategy_i,med_index);
                f_x_med_f7(strategy_i,para_i,:) = f_x_matrix(strategy_i,med_index,:);
                sigma_med_f7(strategy_i,para_i,:) = sigma_matrix(strategy_i,med_index,:);
                sigma_star__med_f7(strategy_i,para_i,:) = sigma_star_matrix(strategy_i,med_index,:);
                success_med_f7(strategy_i,para_i) = success_rate_array(strategy_i,med_index);
                if strategy_i ~= 1
                    four_prob_med_f7(strategy_i,para_i,:) = four_prob_matrix(strategy_i,med_index,:);
                    eval_rate_med_f7(strategy_i,para_i) = eval_rate_array(strategy_i,med_index);
                end
            elseif fname == 8
                T_med_f8(strategy_i,para_i) = T_array(strategy_i,med_index);
                f_x_med_f8(strategy_i,para_i,:) = f_x_matrix(strategy_i,med_index,:);
                sigma_med_f8(strategy_i,para_i,:) = sigma_matrix(strategy_i,med_index,:);
                sigma_star__med_f8(strategy_i,para_i,:) = sigma_star_matrix(strategy_i,med_index,:);
                success_med_f8(strategy_i,para_i) = success_rate_array(strategy_i,med_index);
                if strategy_i ~= 1
                    four_prob_med_f8(strategy_i,para_i,:) = four_prob_matrix(strategy_i,med_index,:);
                    eval_rate_med_f8(strategy_i,para_i) = eval_rate_array(strategy_i,med_index);
                end
            elseif fname == 9
                T_med_f9(strategy_i,para_i) = T_array(strategy_i,med_index);
                f_x_med_f9(strategy_i,para_i,:) = f_x_matrix(strategy_i,med_index,:);
                sigma_med_f9(strategy_i,para_i,:) = sigma_matrix(strategy_i,med_index,:);
                sigma_star__med_f9(strategy_i,para_i,:) = sigma_star_matrix(strategy_i,med_index,:);
                success_med_f9(strategy_i,para_i) = success_rate_array(strategy_i,med_index);
                if strategy_i ~= 1
                    four_prob_med_f9(strategy_i,para_i,:) = four_prob_matrix(strategy_i,med_index,:);
                    eval_rate_med_f9(strategy_i,para_i) = eval_rate_array(strategy_i,med_index);
                end
            end            
        end
    end
    fprintf('fname = %d done\n',fname);
end

save_name_sprint = sprintf('med_dim=%d.mat',n);
save(save_name_sprint,'NUM_OF_STRATEGIES','n','NUM_OF_RUNS','f6_range','f7_range', 'f8_range',...
    'TRAINING_SIZE','LS_onePlusOne','LS_mml', 'NUM_OF_ITERATIONS','C1','C2','C3',...
    'T_med_f6','f_x_med_f6','sigma_med_f6','sigma_star__med_f6','success_med_f6','four_prob_med_f6','eval_rate_med_f6',...
    'T_med_f7','f_x_med_f7','sigma_med_f7','sigma_star__med_f7','success_med_f7','four_prob_med_f7','eval_rate_med_f7',...
    'T_med_f8','f_x_med_f8','sigma_med_f8','sigma_star__med_f8','success_med_f8','four_prob_med_f8','eval_rate_med_f8',...
    'T_med_f9','f_x_med_f9','sigma_med_f9','sigma_star__med_f9','success_med_f9','four_prob_med_f9','eval_rate_med_f9');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(FIGURE_NUM);



%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig. 1
subplot(subplot_ROW,subplot_COL,(fig_row_index-1)*subplot_COL+1)
plot(f6_range, T_med_f6(1,:)./T_med_f6(2,:));hold on;
plot(f6_range, T_med_f6(1,:)./T_med_f6(3,:));hold on;
plot(f6_range, T_med_f6(1,:)./T_med_f6(4,:));hold on;
plot(f6_range, T_med_f6(1,:)./T_med_f6(5,:));hold on;
legend({'GP-(1+1)-ES','GP-(3/3,10)-ES','GP-(5/5,20)-ES','GP-(10/10,40)-ES'});
if fig_row_index==1
    title('spheres','fontsize',15);
end
xlabel('exponent/2','FontSize',15); 
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
ylim([1 40]);
ylabel(sprintf('speed-up over (1+1)-ES N=%d',n),'FontSize',13);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig. 2
subplot(subplot_ROW,subplot_COL,(fig_row_index-1)*subplot_COL+2)
plot(f7_range, T_med_f7(1,:)./T_med_f7(2,:));hold on;
plot(f7_range, T_med_f7(1,:)./T_med_f7(3,:));hold on;
plot(f7_range, T_med_f7(1,:)./T_med_f7(4,:));hold on;
plot(f7_range, T_med_f7(1,:)./T_med_f7(5,:));hold on;
legend({'GP-(1+1)-ES','GP-(3/3,10)-ES','GP-(5/5,20)-ES','GP-(10/10,40)-ES'});
if fig_row_index==1
    title('quartic','fontsize',15);
end
set(gca, 'YScale', 'log');
ylim([1 40]);
xlabel('\alpha','FontSize',15); 
% ylabel(sprintf('speed-up over (1+1)-ES N=%d',n),'FontSize',13);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig. 3
subplot(subplot_ROW,subplot_COL,(fig_row_index-1)*subplot_COL+3)
plot(f8_range, T_med_f8(1,:)./T_med_f8(2,:));hold on;
plot(f8_range, T_med_f8(1,:)./T_med_f8(3,:));hold on;
plot(f8_range, T_med_f8(1,:)./T_med_f8(4,:));hold on;
plot(f8_range, T_med_f8(1,:)./T_med_f8(5,:));hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
ylim([1 40]);
legend({'GP-(1+1)-ES','GP-(3/3,10)-ES','GP-(5/5,20)-ES','GP-(10/10,40)-ES'});
if fig_row_index==1
    title('ellipsoids','fontsize',15);
end
xlabel('\beta','FontSize',15); 
% ylabel(sprintf('speed-up over (1+1)-ES N=%d',n),'FontSize',13);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig. 4
subplot(subplot_ROW,subplot_COL,(fig_row_index-1)*subplot_COL+4)
bar([T_med_f9(1,:)./T_med_f9(2,:); T_med_f9(1,:)./T_med_f9(3,:);...
    T_med_f9(1,:)./T_med_f9(4,:); T_med_f9(1,:)./T_med_f9(5,:); ]);hold on;
% plot(f6_range, T_med_f7(1,:)./T_med_f9(2,:));hold on;
% plot(f6_range, T_med_f7(1,:)./T_med_f9(3,:));hold on;
% plot(f6_range, T_med_f7(1,:)./T_med_f9(4,:));hold on;
% plot(f6_range, T_med_f7(1,:)./T_med_f9(5,:));
str_cell_SIGMA_STAR = {'GP-(1+1)','GP-(3/3,10)','GP-(5/5,20)','GP-(10/10,40)'};
set(gca,'xticklabel',str_cell_SIGMA_STAR);
set(gca, 'YScale', 'log');
ylim([1 40]);
if fig_row_index==1
    title('Schwefel','fontsize',15);
end
% ylabel(sprintf('speed-up over (1+1)-ES N=%d',n),'FontSize',13);

fig_name = sprintf('merged_speed-up_dim%d.fig',n);
saveas(gcf,fig_name); 


% 
% for fname = 6:1:9
%     if(fname == 6)
%         dt =sprintf('spheres');
%         title(dt,'fontsize',15);
%     elseif(fname == 7)
%         dt =sprintf('');
%         title(dt,'fontsize',15);
%     elseif(fname == 8)
%         dt =sprintf('');
%         title(dt,'fontsize',15);
%     elseif(fname == 9)
%         dt =sprintf('Schwefel function');
%         title(dt,'fontsize',15);
%     end
%     
% end

% 
% for q = 1:1:length(PROB_RATE_array)
%     lambda = PROB_RATE_array(q);
%     mu = ceil(lambda/4);
%     C1 = SUCCESS_RATE*DF;
%     C2 = (1-SUCCESS_RATE)*DF;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Def of variables
% 
%     
% 
% 
%     
% 
%     
%     % Range of T for ploting
% %     t_start = TRAINING_SIZE/lambda;
% %     T_range_1 = (0:1:t_start).*(lambda+1)+1;
% %     T_range_2 = (t_start+2:t_med)+lambda*t_start;
% %     T_range = [T_range_1 T_range_2];
%     T_range = 1:1:t_med;
%     % counter
%     fprintf('FIG %d finished \n',fname);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % GRPAH
%     % For f(x), sigma, sigmaStar, GP_error over objective function evaluations
%     
%     legend('-DynamicLegend'); 
%     hold on;
%     
%     d = sprintf('(%d/%d,%d)-ES,C=%.1f,%.1f,%.1f,LS=%d',mu,mu,lambda,C1,C2,C3,LENGTH_SCALE);
%     
%     % Fig 1: histgragm of objFunCalls
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     subplot(subplot_ROW,subplot_COL,(0)*subplot_COL+fname);
%     h = histogram(T_array);hold on; % mml with GP
%     
%     if(fname == 1)
%         dt =sprintf('linear sphere');
%         title(dt,'fontsize',15);
%         h.BinWidth = 5;
%     elseif(fname == 2)
%         dt =sprintf('quadratic sphere');
%         title(dt,'fontsize',15);
%         h.BinWidth = 2;
%     elseif(fname == 3)
%         dt =sprintf('cubic sphere');
%         title(dt,'fontsize',15);
%         h.BinWidth = 2;
%     elseif(fname == 4)
%         dt =sprintf('Schwefel function');
%         title(dt,'fontsize',15);
%         h.BinWidth = 10;
%     elseif(fname == 5)
%         dt =sprintf('quartic function');
%         title(dt,'fontsize',15); 
%         h.BinWidth = 10;
%     end
%     xlabel('Objective function calls','FontSize',15); 
%     
%     
%     
%     % Fig 2: convergence plots [FIGURE_NUM.fig]
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     subplot(subplot_ROW,subplot_COL,(1)*subplot_COL+fname);
%     plot(T_range, f_x_med(1:t_med),'DisplayName',d);hold on; % mml with GP
%     if(fname==1)
%         ylabel(sprintf('Objective function value'),'FontSize',15);
%     end
%     xlabel('function calls','FontSize',15); 
%     set(gca, 'YScale', 'log');
%     legend('-DynamicLegend'); 
%     legend('show');
% 
% 
%     % Fig 3: step size [FIGURE_NUM+1.fig]
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     subplot(subplot_ROW,subplot_COL,(2)*subplot_COL+fname);
%     plot(T_range, sigma_matrix_med(1:t_med),'DisplayName',d);hold on; % mml with GP
%     % plot(1:t_med_noGP, f_x_med_noGP(i,1:t_med_noGP),'DisplayName',d_noGP);hold on; % mml with GP
%     if(fname==1)
%         ylabel(sprintf('Step size'),'FontSize',15);
%     end
%     xlabel('function calls','FontSize',15);
%     set(gca, 'YScale', 'log');
%     legend('-DynamicLegend'); 
%     legend('show');
% 
%     % Fig 4: Normalized step size [FIGURE_NUM+2.fig]
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     subplot(subplot_ROW,subplot_COL,(3)*subplot_COL+fname);
%     plot(T_range, sigma_star_med(1:t_med),'DisplayName',d);hold on; % mml with GP
%     % plot(1:t_med_noGP, f_x_med_noGP(i,1:t_med_noGP),'DisplayName',d_noGP);hold on; % mml with GP
%     if(fname==1)
%         ylabel(sprintf('Normalized step size'),'FontSize',15);
%     end
%     xlabel('function calls','FontSize',15);
%     set(gca, 'YScale', 'log');
%     legend('-DynamicLegend'); 
%     legend('show');
% 
% end
%     figure(FIGURE_NUM);
%     subplot(subplot_ROW,subplot_COL,(0)*subplot_COL+fname);
%     legendCell = {};
%     for i = 1:1:length(PROB_RATE_array)
%         legendCell{i} = sprintf('(%d/%d,%d)-ES',ceil(PROB_RATE_array(i)/4),ceil(PROB_RATE_array(i)/4),PROB_RATE_array(i));
%     end
%     legend(legendCell);
%     
%     % Fig 5:med success rate & evaluation rate [bar][bar]
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     subplot(subplot_ROW,subplot_COL,(4)*subplot_COL+fname);
%     bar([success_rate_final_med, eval_rate_final_med]);hold on;
%     set(gca,'xticklabel',str_cell_SIGMA_STAR);
%     if(fname==1)
%         ylabel(sprintf('Probabilities'),'FontSize',15);
%     end
%     legend({'success rate','evaluation rate'});
%     xlabel('\lambda','FontSize',15);  
% 
%     % Fig 6: four probs TP, TN, FP, FN
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     subplot(subplot_ROW,subplot_COL,(5)*subplot_COL+fname);
%     bar(four_prob_final_med./sum(four_prob_final_med,2));hold on;
%     legend({'TN','FP','FN','TP'})
%     set(gca,'xticklabel',str_cell_SIGMA_STAR);
%     if(fname==1)
%         ylabel(sprintf('Probabilities'),'FontSize',15);
%     end
%     xlabel('\lambda','FontSize',15); 
%     fig_name = sprintf('merged%d_LS.fig',lambda);
%     saveas(gcf,fig_name); 
% 
%     val = {T_final_med};
% % val = {t_array,sigma_matrix,T_array,f_x_matrix,convergence_rate_array,GP_error_matrix,sigma_star_matrix,success_rate_array,delta_matrix};
end