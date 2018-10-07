% test file
% using 5 different obejctive functions plot 5 seperate graphs
% NOTE: DO NOT CHANGE THE NUMBER WHEN CALL fun_graph_merged
%       OR CANNOT SAVE ITERTAION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = @(x) (x'*x)^(1/2);
f2 = @(x) (x'*x);
f3 = @(x) (x'*x)^(3/2);


strategy_name = '[centroidQuadratic]';
lambda_start = 1;
lambda_end = 50;
lambda_increment = 20;

NUM_OF_RUNS = 1;
sigma_star_start = 0.2;
sigma_star_end = 4.2;
sigma_increment = 0.8;
sigma_star_array = sigma_star_start:sigma_increment:sigma_star_end;
lambda_array = [30 50];
theta_array = [2,4,8];
% theta_array = [2,4,8,16,32]



disp("===========================================================");
disp("linear sphere");
disp("---------------");
%graph_fun(f1,1);
TRAINING_FACTOR = 3;
T_f11=fun_c_over_thetaAndSigmaStar(f1,6,NUM_OF_RUNS,sigma_star_array,lambda_array,theta_array,TRAINING_FACTOR,strategy_name);

% 
% disp("===========================================================");
% disp("quadratic sphere");
% disp("---------------");
% TRAINING_FACTOR = 3;
% T_f21=fun_mergedGraph_fx_funCalls(f2,7,NUM_OF_RUNS,sigma_star_array,lambda_array,TRAINING_FACTOR,strategy_name);
% 
% 
% disp("===========================================================");
% disp("cubic sphere");
% disp("---------------");
% TRAINING_FACTOR = 3;
% T_f31=fun_mergedGraph_fx_funCalls(f3,8,NUM_OF_RUNS,sigma_star_array,lambda_array,TRAINING_FACTOR,strategy_name);
% 
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