% test file [each lambda is a fig.]
% using 5 different obejctive functions plot 5 seperate graphs
% NOTE: DO NOT CHANGE THE NUMBER WHEN CALL fun_graph_merged
%       OR CANNOT SAVE ITERTAION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = @(x) (x'*x)^(1/2);
f2 = @(x) (x'*x);
f3 = @(x) (x'*x)^(3/2);

close all;
FIGURE_NUM = 1;
NUM_OF_RUNS =2;
% NUM_OF_RUNS = 2;
TRAINING_SIZE = 40;
LENGTH_SCALE = 8;



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
% lambda_array = [10,20,40];
lambda_array = [10,20,40];

SIGMA_STAR_array = 2.^(-2:1:2);
C1 = 0.6;
C2 = 0.6;
C3_array = 0.05:0.05:0.3;
PROB_array = C3_array;

subplot_ROW = 5;
subplot_COL = 5;
str_cell_PROB_RATE = cellstr(num2str(PROB_array'));

T_matrix = zeros(length(PROB_array),length(lambda_array),5);             % dim = [sigmaStar, lambda, fname]
% str_cell_SIGMA_STAR = [];       % store sigmaStar in a cell
% for i = 1:1:length(2.^(0:1:3))
%     temp = SIGMA_STAR_array(i);
%     str_cell_SIGMA_STAR = [str_cell_SIGMA_STAR; num2str(temp)];
% end
% str_cell_SIGMA_STAR = {str_cell_SIGMA_STAR};
for i=1:1:length(lambda_array)
    lambda = lambda_array(i);
	FIGURE_NUM = i;
    for fname = 1:1:5
        temp = fun_multi_change_C3_FUNCALLS(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE,C1,C2,PROB_array,FIGURE_NUM,subplot_ROW,subplot_COL,str_cell_PROB_RATE);
        T_matrix(:,i,fname) = cell2mat(temp);
    end
end

save('med_T.mat','lambda_array','SIGMA_STAR_array','T_matrix','LENGTH_SCALE','TRAINING_SIZE','NUM_OF_RUNS','C1','C2','C3_array','NUM_OF_RUNS');

