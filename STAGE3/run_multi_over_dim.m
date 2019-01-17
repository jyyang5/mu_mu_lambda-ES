

n_dim_array = [4,8,10,16];
close all;
n = 10;
NUM_OF_RUNS = 21; 
% f6_range = 10.^(-1:1:1);
% f7_range = 1:5:10;
% f8_range = 10.^(-2:2:2);
f6_range = 10.^(-1:1/4:1);
f7_range = 1:0.5:5;
f8_range = 10.^(-2:0.5:2);
TRAINING_SIZE = 40;
LS_onePlusOne = 8;
LS_mml = 20;
NUM_OF_ITERATIONS = 50000;
FIGURE_NUM = 100;
subplot_ROW = length(n_dim_array);
subplot_COL = 4;

fig_row_index = 1;
for i = 1:1:length(n_dim_array)
    n = n_dim_array(i);
    fun_multi_over_dim(n,NUM_OF_RUNS,f6_range, f7_range, f8_range,TRAINING_SIZE,LS_onePlusOne,LS_mml,NUM_OF_ITERATIONS,FIGURE_NUM,subplot_ROW,subplot_COL,i);
end