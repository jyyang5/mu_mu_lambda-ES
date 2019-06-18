

n_dim_array = [4,8,12];
close all;
n = 10;
NUM_OF_RUNS = 3; 

C1 = 1.0;
C2 = 1.0;
C3 = 0.1;
LS_mml = 16;

% f6_range = [1,2,3];%10.^(-1:1:1);
% f7_range = [1];%1:5:10;
% f8_range = [2];%10.^(-2:2:2);

% f6_range = 10.^(-1:1/4:1);
% f7_range = 1:0.5:5;
% f8_range = 10.^(-2:0.5:2);
f6_range = 10.^(-0.5:1/4:1);
f7_range = 1;
f8_range = 1;
TRAINING_SIZE = 40;
LS_onePlusOne = 8;
NUM_OF_ITERATIONS = 50000;
FIGURE_NUM = 100;
subplot_ROW = length(n_dim_array);
subplot_COL = 4;
kappa = 3;
fig_row_index = 1;
for i = 1:1:length(n_dim_array)
    n = n_dim_array(i);
    fun_multi_over_dim_w(n,NUM_OF_RUNS,f6_range, f7_range, f8_range,TRAINING_SIZE,...
        LS_onePlusOne,LS_mml,NUM_OF_ITERATIONS,FIGURE_NUM,subplot_ROW,subplot_COL,...
        i,C1,C2,C3,kappa);
end