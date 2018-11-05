% test file
% using 5 different obejctive functions plot 5 seperate graphs
% NOTE: DO NOT CHANGE THE NUMBER WHEN CALL fun_graph_merged
%       OR CANNOT SAVE ITERTAION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = @(x) (x'*x)^(1/2);
f2 = @(x) (x'*x);
f3 = @(x) (x'*x)^(3/2);

close all;
NUM_OF_RUNS = 200;
% NUM_OF_RUNS = 2;
TRAINING_SIZE = 40;
LENGTH_SCALE = 8;
% lambda_array = [5 10 15 20 25 30 35 40 50 60 70 80 90 100];

% 
% fname = 1;
% data11 = fun_multi_run(fname,NUM_OF_RUNS,0,TRAINING_SIZE,LENGTH_SCALE);
% data12 = fun_multi_run(fname,NUM_OF_RUNS,10,TRAINING_SIZE,LENGTH_SCALE);
% data13 = fun_multi_run(fname,NUM_OF_RUNS,20,TRAINING_SIZE,LENGTH_SCALE);
% data14 = fun_multi_run(fname,NUM_OF_RUNS,40,TRAINING_SIZE,LENGTH_SCALE);
% 
% 
% fname = 2;
% data21 = fun_multi_run(fname,NUM_OF_RUNS,0,TRAINING_SIZE,LENGTH_SCALE);
% data22 = fun_multi_run(fname,NUM_OF_RUNS,10,TRAINING_SIZE,LENGTH_SCALE);
% data23 = fun_multi_run(fname,NUM_OF_RUNS,20,TRAINING_SIZE,LENGTH_SCALE);
% data24 = fun_multi_run(fname,NUM_OF_RUNS,40,TRAINING_SIZE,LENGTH_SCALE);
% 
% fname = 3;
% data31 = fun_multi_run(fname,NUM_OF_RUNS,0,TRAINING_SIZE,LENGTH_SCALE);
% data32 = fun_multi_run(fname,NUM_OF_RUNS,10,TRAINING_SIZE,LENGTH_SCALE);
% data33 = fun_multi_run(fname,NUM_OF_RUNS,20,TRAINING_SIZE,LENGTH_SCALE);
% data34 = fun_multi_run(fname,NUM_OF_RUNS,40,TRAINING_SIZE,LENGTH_SCALE);
% 
% fname = 4;
% data41 = fun_multi_run(fname,NUM_OF_RUNS,0,TRAINING_SIZE,LENGTH_SCALE);
% data42 = fun_multi_run(fname,NUM_OF_RUNS,10,TRAINING_SIZE,LENGTH_SCALE);
% data43 = fun_multi_run(fname,NUM_OF_RUNS,20,TRAINING_SIZE,LENGTH_SCALE);
% data44 = fun_multi_run(fname,NUM_OF_RUNS,40,TRAINING_SIZE,LENGTH_SCALE);
% 
% fname = 5;
% data51 = fun_multi_run(fname,NUM_OF_RUNS,0,TRAINING_SIZE,LENGTH_SCALE);
% data52 = fun_multi_run(fname,NUM_OF_RUNS,10,TRAINING_SIZE,LENGTH_SCALE);
% data53 = fun_multi_run(fname,NUM_OF_RUNS,20,TRAINING_SIZE,LENGTH_SCALE);
% data54 = fun_multi_run(fname,NUM_OF_RUNS,40,TRAINING_SIZE,LENGTH_SCALE);

t_array = zeros(5,4,NUM_OF_RUNS,1);                                           % # of iterations for the stop creteria
sigma_matrix = zeros(5,4,NUM_OF_RUNS,10000);                                  % store all sigma
T_array = zeros(5,4,NUM_OF_RUNS,1);                                           % # of objective function evaluations for the stop creteria
f_x_matrix = zeros(5,4,NUM_OF_RUNS,10000);                                    % store all fx
convergence_rate_array = zeros(5,4,NUM_OF_RUNS,1);                            % convergence rate 
GP_error_matrix = zeros(5,4,NUM_OF_RUNS,10000);                               % store similar to noise-to-signal ratio
sigma_star_matrix = zeros(5,4,NUM_OF_RUNS,10000);                             % normalized step size 
success_rate_array = zeros(5,4,NUM_OF_RUNS,1);                                % success rate 
delta_matrix = zeros(5,4,NUM_OF_RUNS,10000);                                  % each [i,j] stores a delta array 



lambda_array = [0,10,20,40];
for fname = 1:1:5
    for i = 1:1:length(lambda_array)
        lambda = lambda_array(i);
        temp = fun_multi_run(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE);
        t_array(fname,i,:) =  cell2mat(temp(1));
        sigma_matrix(fname,i,:,:) = cell2mat(temp(2));
        T_array(fname,i,:) = cell2mat(temp(3));                              % # of objective function evaluations for the stop creteria
        f_x_matrix(fname,i,:,:) = cell2mat(temp(4));                         % store all fx
        convergence_rate_array(fname,i,:) = cell2mat(temp(5));               % convergence rate 
        GP_error_matrix(fname,i,:,:) = cell2mat(temp(6));                    % store similar to noise-to-signal ratio
        sigma_star_matrix(fname,i,:,:) = cell2mat(temp(7));                  % normalized step size 
        success_rate_array(fname,i,:) = cell2mat(temp(8));                   % success rate 
        delta_matrix(fname,i,:,:) = cell2mat(temp(9));                       % each [i,j] stores a delta array 
    end   
end





% 
% %Get data
% a = cell2mat(temp1);
% b = cell2mat(temp2);
% c = cell2mat(temp3);
% % d = cell2mat(temp4);
% % e = cell2mat(temp5);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot summary
% figure(100)
% 
% % plot objective function evaluation
% subplot(3,1,3)
% hold on;
% plot(lambda_array,a(:,1));hold on;
% plot(lambda_array,b(:,1));hold on;
% plot(lambda_array,c(:,1));hold on;
% % plot(lambda_array,d(1));hold on;
% % plot(lambda_array,e(1));hold on;
% xlabel('Population size \lambda','fontsize',15);
% ylabel('Objective function evaluation','fontsize',15);
% legend('Linear sphere','Quadratic sphere','Cubic sphere');
% % legend('Linear sphere','Quadratic sphere','Cubic sphere','Schwefel?s function','Quartic function')
% d =sprintf('Objective function evaluation');
% title(d,'fontsize',20);
% 
% % plot convergence rate
% subplot(3,1,1)
% hold on;
% plot(lambda_array,a(:,2));hold on;
% plot(lambda_array,b(:,2));hold on;
% plot(lambda_array,c(:,2));hold on;
% % plot(lambda_array,d(2));hold on;
% % plot(lambda_array,e(2));hold on;
% xlabel('Population size \lambda','fontsize',15);
% ylabel('Convergence rate','fontsize',15);
% legend('Linear sphere','Quadratic sphere','Cubic sphere');%,'Schwefel?s function','Quartic function')
% d =sprintf('Convergence rate');
% title(d,'fontsize',20);
% 
% % plot success rate
% subplot(3,1,2)
% hold on;
% plot(lambda_array,a(:,3));hold on;
% plot(lambda_array,b(:,3));hold on;
% plot(lambda_array,c(:,3));hold on;
% % plot(lambda_array,d(3));hold on;
% % plot(lambda_array,e(3));hold on;
% xlabel('Population size \lambda','fontsize',15);
% ylabel('Objective function evaluation','fontsize',15);
% legend('Linear sphere','Quadratic sphere','Cubic sphere');
% % legend('Linear sphere','Quadratic sphere','Cubic sphere','Schwefel?s function','Quartic function')
% d =sprintf('Objective function evaluation');
% title(d,'fontsize',20);
% 
% % Saves fig
% saves(gcf,'summary.fig');
% % 
% % % linear sphere
% % disp("===========================================================");
% % disp("linear sphere");
% % disp("---------------");
% % name = 6;
% % result_f1=fun_graph_funCall_NEW(f1,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% % 
% % % quadratic sphere
% % disp("===========================================================");
% % disp("quadratic sphere");
% % disp("---------------");
% % name = 7;
% % result_f2 = fun_graph_funCall_NEW(f2,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% % 
% % % cubic sphere
% % disp("===========================================================");
% % disp("cubic sphere");
% % disp("---------------");
% % name = 8;
% % result_f3=fun_graph_funCall_NEW(f3,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% 
% % disp("===========================================================");
% % disp("schwefel's function");
% % disp("---------------");
% % name = 9;
% % result_f4=fun_graph_funCall_NEW(@f4,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% % 
% % disp("===========================================================");
% % disp("quartic function");
% % disp("---------------");
% % name = 10;
% % result_f5=fun_graph_funCall_NEW(@f5,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% % disp("===========================================================");
% 
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Display statics
% % plotValue(result_f1,6);
% % plotValue(result_f2,7);
% % plotValue(result_f3,8);
% % plotValue(result_f4,9);
% % plotValue(result_f5,10);
% % % save('final_result.mat','result_f1','result_f2','result_f3','NUM_OF_RUNS','TRAINING_SIZE','mu','lambda');
% % save('final_result.mat','result_f1','result_f2','result_f3','result_f4','result_f5','NUM_OF_RUNS','TRAINING_SIZE','mu','lambda');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [val] = f4(x)
%     val = 0;
%     for i = 1:length(x)
%         temp = 0;
%         for j = 1:i
%             temp = temp + x(j);
%         end
%         val = val + temp^2;          
%     end
% end
% 
% function [val] = f5(x)
%     val = 0;
%     beta = 1;
%     for i = 1:length(x)-1
%     val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function plotValue(temp,name)
%     if name == 6
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%         disp('Linear sphere');
%     elseif name == 7
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%         disp('Quadratic sphere');    
%     elseif name == 8
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%         disp('Cubic sphere'); 
%     elseif name == 9 
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%         disp('Schwefel function');    
%     elseif name == 10
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%         disp('Quartic function');
%     end
%         
%         
%     T = cell2mat(temp(1));
%     T1 = cell2mat(temp(2));
%     convergence = cell2mat(temp(3));
%     convergence1 = cell2mat(temp(4));
%     success = cell2mat(temp(5));
%     success1 = cell2mat(temp(6));
%     disp('# objective function evaluations (with model assistance)');
%     d = sprintf('mml-ES: %d',T);
%     disp(d);
%     d = sprintf('(1+1)-ES: %d',T1);
%     disp(d);
%     disp('Convergence rate (with model assistance)');
%     d = sprintf('mml-ES: %.3f',convergence);
%     disp(d);
%     d = sprintf('(1+1)-ES: %.3f',convergence1);
%     disp(d);
%     disp('Success rate (with model assistance)');
%     d = sprintf('mml-ES: %.3f',success);
%     disp(d);
%     d = sprintf('(1+1)-ES: %.3f',success1);
%     disp(d);
% 
% end