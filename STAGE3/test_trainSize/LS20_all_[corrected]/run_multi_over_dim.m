

n_dim_array = [2, 4, 8, 16];
% close all;
NUM_OF_RUNS = 51; 

C1 = 1.0;
C2 = 1.0;
C3 = 0.2;
LS_mml = 20;

% Standard [renew]
% f6_range = [10.^(-1:1/5:0),2,3, 4,6,10 ];
% f7_range = 1:0.4:5;
% f8_range = 10.^(-2:0.4:2);


f6_range =  [0.25 10.^(-0.6+0.12:0.12:0.6-0.12*3) 2 3 4];
f7_range = 1:0.4:5;
f8_range = 10.^(-2:0.4:2);

% f6_range =  [0.25 0.5 1 2 3 4];
% f7_range = [1 2];
% f8_range = [1 2];

% TRAINING_SIZE = 40;
LS_onePlusOne = 20;
NUM_OF_ITERATIONS = 50000;
FIGURE_NUM = 2;
subplot_ROW =3;
subplot_COL =  length(n_dim_array);

fig_row_index = 1;

for i = 1:1:length(n_dim_array)
    n = n_dim_array(i);
    TRAINING_SIZE = 20+2*n;
    fun_multi_over_dim(n,NUM_OF_RUNS,f6_range, f7_range, f8_range,TRAINING_SIZE,...
        LS_onePlusOne,LS_mml,NUM_OF_ITERATIONS,FIGURE_NUM,subplot_ROW,subplot_COL,...
        i,C1,C2,C3);
end


