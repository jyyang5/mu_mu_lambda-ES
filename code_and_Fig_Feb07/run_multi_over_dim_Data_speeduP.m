% save data and speed-up Fig.
n_dim_array = [2,4,8,16];
% close all;
NUM_OF_RUNS = 51; % Number of replicates

C1 = 1.0;
C2 = 1.0;
C3 = 0.2;
LS_mml = 8;
LS_onePlusOne = LS_mml;




% \alpha for spheres
f6_range =  2.^(-2:0.4:2);
% % \gamma for quartic functions 
f7_range = 1:0.4:5;
% % \beta for sellipsoids  
f8_range = 10.^(-2:0.4:2);
% f7_range = [1];
% f8_range = [0.01 100];

% TRAINING_SIZE = 40;
NUM_OF_ITERATIONS = 50000;                          % max nunber of iterations
FIGURE_NUM = 3;                                               % fig. number 
subplot_COL = length(n_dim_array);                  % number of columns in Fig. 
subplot_ROW = 4;                                               % number of test functions + 1 (1st row for legend) NOTE: need to remove the lengend from second row 4th col to top middle

for i = 1:1:length(n_dim_array)
    n = n_dim_array(i);
    TRAINING_SIZE = 8*n;
    LS_onePlusOne = 8*n;
    LS_mml = 8*n;
    fun_multi_over_dim(n,n_dim_array,NUM_OF_RUNS,f6_range, f7_range, f8_range,TRAINING_SIZE,...
        LS_onePlusOne,LS_mml,NUM_OF_ITERATIONS,FIGURE_NUM,subplot_ROW,subplot_COL,...
        C1,C2,C3);
end


