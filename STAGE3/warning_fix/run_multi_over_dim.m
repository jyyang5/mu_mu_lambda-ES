% Plot same dim different setting of TRAIN_SIZE and LENGTH_SCALE together
%%%%%%%%%%%%%%%%%%
% N = 4 
% n_dim_array = [4,4,4,4];
% LS_array = [16,10,8,6];
% TRAIN_SIZE_array = [9,10,11,12];
% FIGURE_NUM = 1;


% n_dim_array = [8,8,8,8];
% LS_array = [16,10,8,6];
% % TRAIN_SIZE_array = [15,16,18,20];
% FIGURE_NUM = 2;

n_dim_array = [16,16,16,16];
LS_array = [16,10,8,6];
TRAIN_SIZE_array = [24,28,30,34];
FIGURE_NUM = 3;

% close all;
NUM_OF_RUNS = 5; 

C1 = 1.0;
C2 = 1.0;
C3 = 0.2;
% LS_mml = 6;
% LS_onePlusOne = LS_mml;

% Standard [renew]
% f6_range = [10.^(-1:1/5:0),2,3, 4,6,10 ];
% f7_range = 1:0.4:5;
% f8_range = 10.^(-2:0.4:2);



f6_range =  [0.25 10.^(-0.6+0.12:0.12:0.6-0.12*3) 2 3 4];%[0.25,2];%
f7_range = 1:0.4:5;%[];%
f8_range = 10.^(-2:0.4:2);% [];%

% TRAINING_SIZE = 40;
NUM_OF_ITERATIONS = 50000;
subplot_COL = length(n_dim_array);
subplot_ROW = 3;

fig_row_index = 1;
clc
for i = 1:1:length(n_dim_array)
    n = n_dim_array(i);
    TRAINING_SIZE = TRAIN_SIZE_array(i);
    
    LS_mml = LS_array(i);
    LS_onePlusOne = LS_array(i);
    fun_multi_over_dim(n,NUM_OF_RUNS,f6_range, f7_range, f8_range,TRAINING_SIZE,...
        LS_onePlusOne,LS_mml,NUM_OF_ITERATIONS,FIGURE_NUM,subplot_ROW,subplot_COL,...
        i,C1,C2,C3);
end

