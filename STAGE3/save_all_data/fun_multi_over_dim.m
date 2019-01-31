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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save all data
T_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range),NUM_OF_RUNS);
f_x_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range),NUM_OF_RUNS,50000);
sigma_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range),NUM_OF_RUNS,50000);
sigma_star_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range),NUM_OF_RUNS,50000);
success_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range),NUM_OF_RUNS);
four_prob_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range),NUM_OF_RUNS,4);
eval_rate_f6 = zeros(NUM_OF_STRATEGIES,length(f6_range),NUM_OF_RUNS);

T_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range),NUM_OF_RUNS);
f_x_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range),NUM_OF_RUNS,50000);
sigma_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range),NUM_OF_RUNS,50000);
sigma_star_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range),NUM_OF_RUNS,50000);
success_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range),NUM_OF_RUNS);
four_prob_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range),NUM_OF_RUNS,4);
eval_rate_f7 = zeros(NUM_OF_STRATEGIES,length(f7_range),NUM_OF_RUNS);

T_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range),NUM_OF_RUNS);
f_x_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range),NUM_OF_RUNS,50000);
sigma_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range),NUM_OF_RUNS,50000);
sigma_star_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range),NUM_OF_RUNS,50000);
success_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range),NUM_OF_RUNS);
four_prob_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range),NUM_OF_RUNS,4);
eval_rate_f8 = zeros(NUM_OF_STRATEGIES,length(f8_range),NUM_OF_RUNS);

f9_range = [1];
T_f9 = zeros(NUM_OF_STRATEGIES,length(f9_range),NUM_OF_RUNS);
f_x_f9 = zeros(NUM_OF_STRATEGIES,length(f9_range),NUM_OF_RUNS,50000);
sigma_f9 = zeros(NUM_OF_STRATEGIES,length(f9_range),NUM_OF_RUNS,50000);
sigma_star_f9 = zeros(NUM_OF_STRATEGIES,length(f9_range),NUM_OF_RUNS,50000);
success_f9 = zeros(NUM_OF_STRATEGIES,length(f9_range),NUM_OF_RUNS);
four_prob_f9 = zeros(NUM_OF_STRATEGIES,length(f9_range),NUM_OF_RUNS,4);
eval_rate_f9 = zeros(NUM_OF_STRATEGIES,length(f8_range),NUM_OF_RUNS);

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
        % Save all datas
        if fname == 6
            T_f6(:,para_i,:) = T_array;
            f_x_f6(:,para_i,:,:) = f_x_matrix;
            sigma_f6(:,para_i,:,:) = sigma_matrix;
            sigma_star_f6(:,para_i,:,:) = sigma_star_matrix;
            success_f6(:,para_i,:) = success_rate_array;
            four_prob_f6(:,para_i,:,:) = four_prob_matrix;
            eval_rate_f6(:,para_i,:) = eval_rate_array;
        elseif fname == 7 
            T_f7(:,para_i,:) = T_array;
            f_x_f7(:,para_i,:,:) = f_x_matrix;
            sigma_f7(:,para_i,:,:) = sigma_matrix;
            sigma_star_f7(:,para_i,:,:) = sigma_star_matrix;
            success_f7(:,para_i,:) = success_rate_array;
            four_prob_f7(:,para_i,:,:) = four_prob_matrix;
            eval_rate_f7(:,para_i,:) = eval_rate_array;
        elseif fname == 8
            T_f8(:,para_i,:) = T_array;
            f_x_f8(:,para_i,:,:) = f_x_matrix;
            sigma_f8(:,para_i,:,:) = sigma_matrix;
            sigma_star_f8(:,para_i,:,:) = sigma_star_matrix;
            success_f8(:,para_i,:) = success_rate_array;
            four_prob_f8(:,para_i,:,:) = four_prob_matrix;
            eval_rate_f8(:,para_i,:) = eval_rate_array;
        elseif fname == 9
            T_f9(:,para_i,:) = T_array;
            f_x_f9(:,para_i,:,:) = f_x_matrix;
            sigma_f9(:,para_i,:,:) = sigma_matrix;
            sigma_star_f9(:,para_i,:,:) = sigma_star_matrix;
            success_f9(:,para_i,:) = success_rate_array;
            four_prob_f9(:,para_i,:,:) = four_prob_matrix;
            eval_rate_f9(:,para_i,:) = eval_rate_array;
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

save_name_sprint = sprintf('allData_dim=%d.mat',n);
save(save_name_sprint,'NUM_OF_STRATEGIES','n','NUM_OF_RUNS','f6_range','f7_range', 'f8_range',...
    'TRAINING_SIZE','LS_onePlusOne','LS_mml', 'NUM_OF_ITERATIONS','C1','C2','C3',...
    'T_med_f6','f_x_med_f6','sigma_med_f6','sigma_star__med_f6','success_med_f6','four_prob_med_f6','eval_rate_med_f6',...
    'T_med_f7','f_x_med_f7','sigma_med_f7','sigma_star__med_f7','success_med_f7','four_prob_med_f7','eval_rate_med_f7',...
    'T_med_f8','f_x_med_f8','sigma_med_f8','sigma_star__med_f8','success_med_f8','four_prob_med_f8','eval_rate_med_f8',...
    'T_med_f9','f_x_med_f9','sigma_med_f9','sigma_star__med_f9','success_med_f9','four_prob_med_f9','eval_rate_med_f9',....
    'T_f6','f_x_f6','sigma_f6','sigma_star_f6','success_f6','four_prob_f6','eval_rate_f6',...
    'T_f7','f_x_f7','sigma_f7','sigma_star_f7','success_f7','four_prob_f7','eval_rate_f7',...
    'T_f8','f_x_f8','sigma_f8','sigma_star_f8','success_f8','four_prob_f8','eval_rate_f8',...
    'T_f9','f_x_f9','sigma_f9','sigma_star_f9','success_f9','four_prob_f9','eval_rate_f9');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(FIGURE_NUM);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig. 1
subplot(subplot_ROW,subplot_COL,(fig_row_index-1)*subplot_COL+1)
for i = 2:1:NUM_OF_STRATEGIES
    plot(f6_range, T_med_f6(1,:)./T_med_f6(i,:));hold on;
end
legend({'GP-(1+1)-ES','GP-(3/3,10)-ES','GP-(5/5,20)-ES','GP-(10/10,40)-ES'});
if fig_row_index==1
    title('spheres','fontsize',15);
end
xlabel('exponent/2','FontSize',15); 
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
ylim([1 40]);
ylabel(sprintf('speed-up over (1+1)-ES N=%d',n),'FontSize',13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig. 2
subplot(subplot_ROW,subplot_COL,(fig_row_index-1)*subplot_COL+2)
for i = 2:1:NUM_OF_STRATEGIES
    plot(f7_range, T_med_f7(1,:)./T_med_f7(i,:));hold on;
end
legend({'GP-(1+1)-ES','GP-(3/3,10)-ES','GP-(5/5,20)-ES','GP-(10/10,40)-ES'});
if fig_row_index==1
    title('quartic','fontsize',15);
end
set(gca, 'YScale', 'log');
ylim([1 40]);
xlabel('\alpha','FontSize',15); 
% ylabel(sprintf('speed-up over (1+1)-ES N=%d',n),'FontSize',13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig. 3
subplot(subplot_ROW,subplot_COL,(fig_row_index-1)*subplot_COL+3)
for i = 2:1:NUM_OF_STRATEGIES
    plot(f8_range, T_med_f8(1,:)./T_med_f8(i,:));hold on;
end
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
ylim([1 40]);
legend({'GP-(1+1)-ES','GP-(3/3,10)-ES','GP-(5/5,20)-ES','GP-(10/10,40)-ES'});
if fig_row_index==1
    title('ellipsoids','fontsize',15);
end
xlabel('\beta','FontSize',15); 
% ylabel(sprintf('speed-up over (1+1)-ES N=%d',n),'FontSize',13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig. 4
subplot(subplot_ROW,subplot_COL,(fig_row_index-1)*subplot_COL+4)
% Refactored code is below
% bar([T_med_f9(1,:)./T_med_f9(2,:); T_med_f9(1,:)./T_med_f9(3,:);...
%     T_med_f9(1,:)./T_med_f9(4,:); T_med_f9(1,:)./T_med_f9(5,:); ]);hold on;
bar(repmat(T_med_f9(1,:),NUM_OF_STRATEGIES-1,1)./T_med_f9(2:NUM_OF_STRATEGIES,:));
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

end