% Objective [update step size C1,C2,C3] add hist of objFunCalls [over lambda]
% 1 Plot speed-ups of over (1+1)-ES
%      GP-(1+1)-ES
%      GP-(3/3,10)-ES
%      GP-(5/5,20)-ES
%      GP-(10/10,40)-ES
%    - over test functions: sphere, quartic, ellipsoid, and schewefel
%    - over dimension of data: [2,4,8,16]
% 2. Save
%      Data 
%      plots
%
% difficulty: 
%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.04], [0.08 0.04], [0.08 0.04]);
if ~make_it_tight,  clear subplot;  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIGURE_index = 101;

n_array = [2,4,8,16];
for n = n_array
    if n==2
        load('data_dim=2.mat');
        fig_col_index = 1;
    elseif n==4
        load('data_dim=4.mat')
        fig_col_index = 2;
    elseif n==8
        load('data_dim=8.mat')
        fig_col_index = 3;
    elseif n==16
        load('data_dim=16.mat')
        fig_col_index = 4;
    end
        
FIGURE_NUM = FIGURE_index;
figure(FIGURE_NUM);
subplot_COL = 4;
subplot_ROW = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig. 1 [shperes]
subplot(subplot_ROW,subplot_COL,fig_col_index+subplot_COL*1)
for i = 2:1:NUM_OF_STRATEGIES
    plot(f6_range, T_med_f6(1,:)./T_med_f6(i,:));hold on;
end
% legend({'(1,1)','(3/3,10)','(5/5,20)','(10/10,40)'},'Fontsize',15,'Interpreter','latex','NumColumns',2);
title(sprintf('dimension $n=%d$',n),'Fontsize',15,'Interpreter','latex');
if fig_col_index==1
    ylabel({'speed-up' ;'sphere functions'},'Fontsize',15,'Interpreter','latex');
end
set(gca, 'YScale', 'log', 'XScale', 'log', 'Fontsize',15);
xtickformat('%.1f');
ylim([1 8]);
if fig_col_index==4
    legend({'(1/1,1)','(3/3,10)','(5/5,20)','(10/10,40)'},'Fontsize',15,'Interpreter','latex','NumColumns',4);
end
xlabel('parameter $\alpha$','Fontsize',15,'Interpreter','latex');
yticks([1 2 4 8]);
yticklabels({'1.0','2.0','4.0', '8.0'});
xticks([0.25, 0.5, 1.0, 2.0, 4.0]);
xticklabels({'0.25','0.50', '1.00','2.00', '4.00'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig. 2 [elliposids]
subplot(subplot_ROW,subplot_COL,fig_col_index+subplot_COL*2)
for i = 2:1:NUM_OF_STRATEGIES
    plot(f8_range, T_med_f8(1,:)./T_med_f8(i,:));hold on;
end
ylim([1 32]);
% legend({'(1,1)','(3/3,10)','(5/5,20)','(10/10,40)'},'Fontsize',15,'Interpreter','latex','NumColumns',2);
if fig_col_index==1
    ylabel({'speed-up'; 'ellipsoid functions'},'Fontsize',15,'Interpreter','latex');
end
ytickformat('%.1f');
set(gca, 'YScale', 'log', 'XScale', 'log', 'Fontsize',15);
xlabel('parameter $\beta$','Fontsize',15,'Interpreter','latex');
yticks([1 2 4 8 16 32]);
yticklabels({'1.0','2.0','4.0','8.0','16.0','32.0'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig. 3 [quartic]
subplot(subplot_ROW,subplot_COL,fig_col_index+subplot_COL*3)
for i = 2:1:NUM_OF_STRATEGIES
    plot(f7_range, T_med_f7(1,:)./T_med_f7(i,:));hold on;
end
% legend({'(1,1)','(3/3,10)','(5/5,20)','(10/10,40)'},'Fontsize',15,'Interpreter','latex','NumColumns',2);
if fig_col_index==1
    ylabel({'speed-up'; 'quartic functions'},'Fontsize',15,'Interpreter','latex');
end
set(gca, 'YScale', 'log', 'Fontsize',15);
ylim([1 32]);
xlim([1 5]);
xtickformat('%.1f');
ytickformat('%.1f');
xlabel('parameter $\gamma$','Fontsize',15,'Interpreter','latex');
yticks([1 2 4 8 16 32]);
yticklabels({'1.0','2.0','4.0','8.0','16.0','32.0'});


end

fig_name = sprintf('speed-up1.fig');
saveas(gcf,fig_name); 

