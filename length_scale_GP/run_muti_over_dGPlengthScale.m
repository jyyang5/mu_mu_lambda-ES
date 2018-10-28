% test file
% using 5 different obejctive functions plot 5 seperate graphs
% NOTE: DO NOT CHANGE fname each corespond to a test function
% 
% Input:
%       NUM_OF_RUNS:        number of replicates
%       TRAINING_SIZE:      GP training size
%       LENGTH_SCALE:       length scale factor for GP 
%       lambda_array:       store an array of lambda for plot
%       fname:
% Output:
%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = @(x) (x'*x)^(1/2);
f2 = @(x) (x'*x);
f3 = @(x) (x'*x)^(3/2);

DECREASE_FACTOR = 0.72;
STRATEGY_NAME = @mml_GP_final_emergency;

NUM_OF_RUNS = 1;
% NUM_OF_RUNS = 2;
TRAINING_SIZE = 40;
% LENGTH_SCALE = 8;
% lambda_array = [5 10 15 20 25 40 50 60 80];
lambda=40;
LENGTH_SCALE_array = [2 4 8 16 32 46 64 80 96 112 128 144 256];

fname = 1;
temp1 = fun_multi_over_GPlengthScale(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE_array,DECREASE_FACTOR,STRATEGY_NAME);

fname = 2;
temp2 = fun_multi_over_GPlengthScale(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE_array,DECREASE_FACTOR,STRATEGY_NAME);

fname = 3;
temp3 = fun_multi_over_GPlengthScale(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE_array,DECREASE_FACTOR,STRATEGY_NAME);
% 
% fname = 4;
% temp4 = fun_multi_over_decreaseFactor(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE,DECREASE_FACTOR_array);
% 
% fname = 5;
% temp5 = fun_multi_over_decreaseFactor(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE,DECREASE_FACTOR_array);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get data
a = cell2mat(temp1);
b = cell2mat(temp2);
c = cell2mat(temp3);
% d = cell2mat(temp4);
% e = cell2mat(temp5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot summary
figure(100)

% plot objective function evaluation
subplot(3,1,3)
hold on;
plot(LENGTH_SCALE_array,a(:,1));hold on;
plot(LENGTH_SCALE_array,b(:,1));hold on;
plot(LENGTH_SCALE_array,c(:,1));hold on;
% plot(DECREASE_FACTOR_array,d(:,1));hold on;
% plot(DECREASE_FACTOR_array,e(:,1));hold on;
xlabel('DECREASE_FATOR (when inferior offspring)','fontsize',15);
ylabel('Objective function evaluation','fontsize',15);
legend('Linear sphere','Quadratic sphere','Cubic sphere')%,'Schwefel?s function','Quartic function')
d1 =sprintf('Objective function evaluation');
title(d1,'fontsize',20);

% plot convergence rate
subplot(3,1,1)
hold on;
plot(LENGTH_SCALE_array,a(:,2));hold on;
plot(LENGTH_SCALE_array,b(:,2));hold on;
plot(LENGTH_SCALE_array,c(:,2));hold on;
% plot(DECREASE_FACTOR_array,d(:,2));hold on;
% plot(DECREASE_FACTOR_array,e(:,2));hold on;
xlabel('GP length scale factor','fontsize',15);
ylabel('Convergence rate','fontsize',15);
legend('Linear sphere','Quadratic sphere','Cubic sphere')%,'Schwefel?s function','Quartic function')
d1 =sprintf('Convergence rate');
title(d1,'fontsize',20);

% plot success rate
subplot(3,1,2)
hold on;
plot(LENGTH_SCALE_array,a(:,3));hold on;
plot(LENGTH_SCALE_array,b(:,3));hold on;
plot(LENGTH_SCALE_array,c(:,3));hold on;
% plot(DECREASE_FACTOR_array,d(:,3));hold on;
% plot(DECREASE_FACTOR_array,e(:,3));hold on;
xlabel('GP length scale factor','fontsize',15);
ylabel('Objective function evaluation','fontsize',15);
legend('Linear sphere','Quadratic sphere','Cubic sphere')%,'Schwefel?s function','Quartic function')
d1 =sprintf('Success rate');
title(d1,'fontsize',20);

% Saves fig
saveas(gcf,'summary.fig');

% 
% % linear sphere
% disp("===========================================================");
% disp("linear sphere");
% disp("---------------");
% name = 6;
% result_f1=fun_graph_funCall_NEW(f1,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% 
% % quadratic sphere
% disp("===========================================================");
% disp("quadratic sphere");
% disp("---------------");
% name = 7;
% result_f2 = fun_graph_funCall_NEW(f2,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% 
% % cubic sphere
% disp("===========================================================");
% disp("cubic sphere");
% disp("---------------");
% name = 8;
% result_f3=fun_graph_funCall_NEW(f3,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);

% disp("===========================================================");
% disp("schwefel's function");
% disp("---------------");
% name = 9;
% result_f4=fun_graph_funCall_NEW(@f4,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% 
% disp("===========================================================");
% disp("quartic function");
% disp("---------------");
% name = 10;
% result_f5=fun_graph_funCall_NEW(@f5,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% disp("===========================================================");

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display statics
% plotValue(result_f1,6);
% plotValue(result_f2,7);
% plotValue(result_f3,8);
% plotValue(result_f4,9);
% plotValue(result_f5,10);
% % save('final_result.mat','result_f1','result_f2','result_f3','NUM_OF_RUNS','TRAINING_SIZE','mu','lambda');
% save('final_result.mat','result_f1','result_f2','result_f3','result_f4','result_f5','NUM_OF_RUNS','TRAINING_SIZE','mu','lambda');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotValue(temp,name)
    if name == 6
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('Linear sphere');
    elseif name == 7
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('Quadratic sphere');    
    elseif name == 8
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('Cubic sphere'); 
    elseif name == 9 
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('Schwefel function');    
    elseif name == 10
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('Quartic function');
    end
        
        
    T = cell2mat(temp(1));
    T1 = cell2mat(temp(2));
    convergence = cell2mat(temp(3));
    convergence1 = cell2mat(temp(4));
    success = cell2mat(temp(5));
    success1 = cell2mat(temp(6));
    disp('# objective function evaluations (with model assistance)');
    d = sprintf('mml-ES: %d',T);
    disp(d);
    d = sprintf('(1+1)-ES: %d',T1);
    disp(d);
    disp('Convergence rate (with model assistance)');
    d = sprintf('mml-ES: %.3f',convergence);
    disp(d);
    d = sprintf('(1+1)-ES: %.3f',convergence1);
    disp(d);
    disp('Success rate (with model assistance)');
    d = sprintf('mml-ES: %.3f',success);
    disp(d);
    d = sprintf('(1+1)-ES: %.3f',success1);
    disp(d);

end