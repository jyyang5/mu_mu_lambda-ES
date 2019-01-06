% test file [for CSA with emergency store emergency rate]
% using 5 different obejctive functions plot 5 seperate graphs
% NOTE: DO NOT CHANGE THE NUMBER WHEN CALL fun_graph_merged
%       OR CANNOT SAVE ITERTAION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = @(x) (x'*x)^(1/2);
f2 = @(x) (x'*x);
f3 = @(x) (x'*x)^(3/2);

close all;
FIGURE_NUM = 1;
NUM_OF_RUNS = 5;
% NUM_OF_RUNS = 2;
TRAINING_SIZE = 40;
LENGTH_SCALE = 16;



t_array = zeros(5,4,NUM_OF_RUNS,1);                                           % # of iterations for the stop creteria
sigma_matrix = zeros(5,4,NUM_OF_RUNS,10000);                                  % store all sigma
T_array = zeros(5,4,NUM_OF_RUNS,1);                                           % # of objective function evaluations for the stop creteria
f_x_matrix = zeros(5,4,NUM_OF_RUNS,10000);                                    % store all fx
convergence_rate_array = zeros(5,4,NUM_OF_RUNS,1);                            % convergence rate 
GP_error_matrix = zeros(5,4,NUM_OF_RUNS,10000);                               % store similar to noise-to-signal ratio
sigma_star_matrix = zeros(5,4,NUM_OF_RUNS,10000);                             % normalized step size 
success_rate_array = zeros(5,4,NUM_OF_RUNS,1);                                % success rate 
delta_matrix = zeros(5,4,NUM_OF_RUNS,10000);                                  % each [i,j] stores a delta array 
% emergency_rate_matrix = zeros(5,4,NUM_OF_RUNS);                             % proportion of emergency [:,1,:] all 0 


% lambda_array = [0,10,20,40];
% success_rate = zeros(3,length(lambda_array),2,LEN_SIGMA_STAR);
lambda_array = [10,20,40];
SUCCESS_RATIO_array = 0.05:0.15:0.95;
subplot_ROW = length(lambda_array);
subplot_COL = 5;

T_matrix = zeros(length(SUCCESS_RATIO_array),length(lambda_array),5);      % [S,lambda,testFunction]

for j=1:1:length(SUCCESS_RATIO_array)
    SUCCESS_RATE = SUCCESS_RATIO_array(j);
    for fname = 1:1:5
        for i = 1:1:length(lambda_array)
            lambda = lambda_array(i);
            temp = fun_multi_run_change_success_rate_FUNCALLS(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE,SUCCESS_RATE,FIGURE_NUM,subplot_ROW,subplot_COL,i);         
            T_matrix(j,i,fname) = cell2mat(temp);
    %         success_rate(fname,i,1,:) = cell2mat(temp(1));
    %         success_rate(fname,i,2,:) = cell2mat(temp(2));
        end   
    end
end
save('med_T.mat','lambda_array','SUCCESS_RATIO_array','T_matrix');


% 
% subplot_ROW = 4;
% subplot_COL = 3;
% figure(FIGURE_NUM+1);
% legend('-DynamicLegend'); 
% hold on;
% mu = ceil(lambda/4);
% d = sprintf('GP-(%d/%d,%d)-ES',mu,mu,lambda);
% d_noGP = sprintf('(%d/%d,%d)-ES',mu,mu,lambda);
% c = categorical({'(3/3,10)-ES','(5/5,20)-ES','(10/10,40)-ES'});
% for fname = 1:1:3
%     for i=1:1:LEN_SIGMA_STAR
% 
%         subplot(subplot_ROW,subplot_COL,(i-1)*subplot_COL+fname);
%         bar(squeeze(success_rate(fname,:,:,i)));
%         legend('with GP','no GP');
%         set(gca,'xticklabel',{'(3/3,10)-ES','(5/5,20)-ES','(10/10,40)-ES'});
%         if i==1
%             if(fname == 1)
%                 dt =sprintf('linear sphere');
%                 title(dt,'fontsize',15);
%             elseif(fname == 2)
%                 dt =sprintf('quadratic sphere');
%                 title(dt,'fontsize',15);
%             elseif(fname == 3)
%                 dt =sprintf('cubic sphere');
%                 title(dt,'fontsize',15);
%             elseif(fname == 4)
%                 dt =sprintf('Schwefel function');
%                 title(dt,'fontsize',15);
%             elseif(fname == 5)
%                 dt =sprintf('quartic function');
%                 title(dt,'fontsize',15); 
%             end
%         end
% %         plot(1:t_med(i), f_x_med(i,1:t_med(i)),'DisplayName',d);hold on; % mml with GP
% %         plot(1:t_med_noGP(i), f_x_med_noGP(i,1:t_med_noGP(i)),'DisplayName',d_noGP);hold on; % mml with GP
%         if(fname==1)
%             ylabel('Success rate','FontSize',15);%
%         end
%         legend('-DynamicLegend'); 
%         legend('show');
%         ylim([0 1]);
%     end
% end
% saveas(gcf,'success_rate_SIGMA_STAR.fig'); 

% save('data_SIGMA_STAR.mat','NUM_OF_RUNS','TRAINING_SIZE','LENGTH_SCALE','SIGMA_STAR',...
%     't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array','GP_error_matrix',...
%     'sigma_star_matrix','success_rate_array','delta_matrix');



