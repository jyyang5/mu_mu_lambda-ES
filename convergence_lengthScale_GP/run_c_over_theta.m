% test file
% using 5 different obejctive functions plot 5 seperate graphs
% NOTE: DO NOT CHANGE THE NUMBER WHEN CALL fun_graph_merged
%       OR CANNOT SAVE ITERTAION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = @(x) (x'*x)^(1/2);
f2 = @(x) (x'*x);
f3 = @(x) (x'*x)^(3/2);


strategy_name = '[theta]';
lambda_start = 1;
lambda_end = 50;
lambda_increment = 20;

NUM_OF_RUNS = 2;
sigma_star_start = 1;
sigma_star_end = 3;
sigma_increment = 0.5;

% sigma_star_array = sigma_star_start:sigma_increment:sigma_star_end;
sigma_star_array = 1:0.5:3.0;
theta_array = 2:2:32;
lambda_array = [30 50 70];
[temp cols] = size(lambda_array);

% NUM_OF_RUNS = 2;
% sigma_star_array = [0.1 0.2 0.5 1 2 4 5 10];

disp("===========================================================");
disp("linear sphere");
disp("---------------");
TRAINING_FACTOR = 0;
lambda = 30;
subplotIndex = [1 cols 1];
T_f11=fun_c_over_theta(f1,6,NUM_OF_RUNS,sigma_star_array,lambda,theta_array,TRAINING_FACTOR,strategy_name,subplotIndex);
lambda = 50;
subplotIndex = [1 cols 2];
T_f12=fun_c_over_theta(f1,6,NUM_OF_RUNS,sigma_star_array,lambda,theta_array,TRAINING_FACTOR,strategy_name,subplotIndex);
lambda = 70;
subplotIndex = [1 cols 3];
T_f13=fun_c_over_theta(f1,6,NUM_OF_RUNS,sigma_star_array,lambda,theta_array,TRAINING_FACTOR,strategy_name,subplotIndex);


% disp("===========================================================");
% disp("quadratic sphere");
% disp("---------------");
% lambda = 30;
% subplotIndex = [1 cols 1];
% T_f21=fun_c_over_theta(f2,7,NUM_OF_RUNS,sigma_star_array,lambda,theta_array,TRAINING_FACTOR,strategy_name,subplotIndex);
% lambda = 50;
% subplotIndex = [1 cols 2];
% T_f22=fun_c_over_theta(f2,7,NUM_OF_RUNS,sigma_star_array,lambda,theta_array,TRAINING_FACTOR,strategy_name,subplotIndex);
% lambda = 70;
% subplotIndex = [1 cols 3];
% T_f23=fun_c_over_theta(f2,7,NUM_OF_RUNS,sigma_star_array,lambda,theta_array,TRAINING_FACTOR,strategy_name,subplotIndex);
% 
% 
% disp("===========================================================");
% disp("cubic sphere");
% disp("---------------");
% lambda = 30;
% subplotIndex = [1 cols 1];
% T_f31=fun_c_over_theta(f3,8,NUM_OF_RUNS,sigma_star_array,lambda,theta_array,TRAINING_FACTOR,strategy_name,subplotIndex);
% lambda = 50;
% subplotIndex = [1 cols 2];
% T_f32=fun_c_over_theta(f3,8,NUM_OF_RUNS,sigma_star_array,lambda,theta_array,TRAINING_FACTOR,strategy_name,subplotIndex);
% lambda = 70;
% subplotIndex = [1 cols 3];
% T_f33=fun_c_over_theta(f3,8,NUM_OF_RUNS,sigma_star_array,lambda,theta_array,TRAINING_FACTOR,strategy_name,subplotIndex);

% 
% disp("===========================================================");
% disp("schwefel's function");
% disp("---------------");
% TRAINING_FACTOR = 3;
% T_f41=fun_mergedGraph_fx_funCalls(@f4,9,NUM_OF_RUNS,sigma_star_array,lambda_array,TRAINING_FACTOR);
% TRAINING_FACTOR = 4;
% T_f42=fun_mergedGraph_fx_funCalls(@f4,9,NUM_OF_RUNS,sigma_star_array,lambda_array,TRAINING_FACTOR);
% 
% 
% disp("===========================================================");
% disp("quartic function");
% disp("---------------");
% TRAINING_FACTOR = 3;
% T_f51=fun_graph_funCall_merged(@f5,10,NUM_OF_RUNS,sigma_star_array,lambda_array);
% TRAINING_FACTOR = 4;
% T_f52=fun_graph_funCall_merged(@f5,10,NUM_OF_RUNS,sigma_star_array,lambda_array);
% 
% 
% disp("===========================================================");


function [val] = f4(x)

val = 0;
for i = 1:length(x)
    temp = 0;
    for j = 1:i
        temp = temp + x(j);
    end
    val = val + temp^2;    
        
end


end

function [val] = f5(x)

val = 0;
beta = 1;
for i = 1:length(x)-1
    val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
end
end