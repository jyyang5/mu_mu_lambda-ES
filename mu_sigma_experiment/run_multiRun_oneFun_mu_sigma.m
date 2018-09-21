% Test parameter E.G. mu, lambda for mml-ES
% run several replicates and take the median as result
% print and save the result to a txt file.

% initialization
% f:                  objective function value
% x0:                 mu initial point
% mu:                 population size
% sigma0:             initial muttaion strength
% NUM_OF_ITERATIONS:  number of maximum iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test functions
f1 = @(x) (x'*x)^(1/2);  % linear sphere
f2 = @(x) (x'*x);        % quadratic sphere
f3 = @(x) (x'*x)^(3/2);  % cubic sphere


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open the file to write data 
fileID = fopen('oneFun_mu_lambda_result.txt','w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialization
n = 10;
NUM_OF_RUNS = 5;
sigma0 = 1;
mu_start = 2;
mu_end = 5;
mu_increment = 1;

lambda_start = 10;
lambda_end = 40;
% lambda_end = 50;
lambda_increment = 10;

f=@f4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% (mu/mu,lambda)-ES with GP
                    % # of iterations for the stop creteria
% f_x_matrix = zeros(NUM_OF_RUNS,10000);             % store all fx
% sigma_matrix = zeros(NUM_OF_RUNS,10000);           % store all sigma
% sigma_final = zeros(NUM_OF_RUNS,1);                % the last sigma
% x_final = zeros(NUM_OF_RUNS,n);                    % last x
% sigma_star_matrix = zeros(NUM_OF_RUNS,10000);      % normalized step size           
% GP_error_matrix = zeros(NUM_OF_RUNS,10000);        % relative error of GP            


mu_length = (mu_end-mu_start)/mu_increment;
lambda_length = (lambda_start-lambda_end)/lambda_increment;

T_array = zeros(mu_length,lambda_length,1);

i = 1;
j = 1;
for mu_temp = mu_start:mu_increment:mu_end
    for lambda_temp = lambda_start:lambda_increment:lambda_end
        x0 = randn(n,mu_temp);
        % (mu/mu,lambda)-ES with GP
        a = fun_multiRun_strategy(f, mu_temp, lambda_temp, NUM_OF_RUNS, n, sigma0);
        T_array(i,j) = cell2mat(a(1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Print result
        d = sprintf('(%d/%d,%d)-ES with GP, T_array[%d,%d]numOfFunCalls = %d',mu_temp,mu_temp,lambda_temp,i,j,T_array(i,j));
        disp(d);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Write data 
        formatSpec = '(%d/%d,%d)-ES with GP, T_array[%d,%d]numOfFunCalls = %d\n';
        fprintf(fileID,formatSpec,mu_temp,mu_temp,lambda_temp,i,j,T_array(i,j));
        
%         disp('number of objective function calls');
%         disp(T_array(i,j));
        j = j + 1;
    end
    i = i+1;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close the file to put data
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save data    
save('oneFun_mu_lambda_data.mat','a','T_array','mu_start','mu_end','mu_increment','lambda_start','lambda_end','lambda_increment');

    
%     % graph
%     figure(name);
%     subplot(1,3,1)
%    
%     hold on;
%     % plot first 40 iterations seperately
%     t_range1 = 1:10:40;
%     t_range2 = 5:T_max;
%     sigma_med_range1 = sigma_med(1:4);
%     sigma_med_range2 = sigma_med(5:T_max);
%     
%     semilogy([t_range1 t_range2], [sigma_med_range1 sigma_med_range2]);% mml with GP
%     semilogy(1:10:10*T_max1, sigma_med1(1:T_max1));  % MML
%     semilogy(1:T_max2, sigma_med2(1:T_max2));        % (1+1)-ES with GP
%     semilogy(1:T_max3, sigma_med3(1:T_max3));        % (1+1)-ES
% 
%     hold off;
%     legend({'mml with GP','mml','(1+1)-ES with GP','(1+1)-ES'},'fontsize',10);
%     xlabel('number of function evaluations','fontsize',15);
%     ylabel('log(\sigma)','fontsize',15);
%     set(gca,'yscale','log')
%     title('step size \sigma','fontsize',20);
%     
%     
%     
%     subplot(1,3,2)
%     hold on;
%     semilogy(1:T_max, f_x_med(1:T_max));             % mml with GP
%     semilogy(1:10:10*T_max1, f_x_med1(1:T_max1));    % mml 
%     semilogy(1:T_max2, f_x_med2(1:T_max2));          % (1+1)-ES with GP
%     semilogy(1:T_max3, f_x_med3(1:T_max3));          % (1+1)-ES
%     
%     hold off;
%     legend({'mml with GP','mml','(1+1)-ES with GP','(1+1)-ES'},'fontsize',10);
%     xlabel('number of function evaluations','fontsize',15);
%     ylabel('log( f(x) )','fontsize',15);
%     set(gca,'yscale','log')
%     title('objective function value f(x)','fontsize',20);
%     
%     
%     subplot(1,3,3)
%     hold on;
%     plot(1:T_max, sigma_star_med(1:T_max));             % mml with GP
%     %plot(1:10:10*T_max1, sigma_star_med1(1:T_max1));    % mml 
%     plot(1:T_max2, sigma_star_med2(1:T_max2));          % (1+1)-ES with GP
%     %plot(1:T_max3, sigma_star_med3(1:T_max3));          % (1+1)-ES
%     
%     hold off;
%     %legend({'mml with GP','mml','(1+1)-ES with GP','(1+1)-ES'},'fontsize',10);
%     legend({'mml with GP','(1+1)-ES with GP'},'fontsize',10);
%     xlabel('number of function evaluations','fontsize',15);
%     ylabel('normalized step size \sigma*','fontsize',15);
%     set(gca,'yscale','log')
%     title('normalized step size \sigma*','fontsize',20);
    
    
% Save data    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if name == 6
%        save('T_linear_sphere.mat','T_array','f_x_med','sigma_med','sigma_star_med','t','GP_error_med');
%        saveas(gcf,'linear_sphere.fig');
%     elseif name == 7
%        save('T_quadratic_sphere.mat','T_array','f_x_med','sigma_med','sigma_star_med','t','GP_error_med');
%        saveas(gcf,'quadratic_sphere.fig');
%     elseif name == 8
%        save('T_cubic_sphere.mat','T_array','f_x_med','sigma_med','sigma_star_med','t','GP_error_med');
%        saveas(gcf,'cubic_sphere.fig');
%     elseif name == 9
%        save('T_schwefel_function.mat','T_array','f_x_med','sigma_med','sigma_star_med','t','GP_error_med');
%        saveas(gcf,'schwefel_function.fig');
%     elseif name == 10
%        save('T_quartic_function.mat','T_array','f_x_med','sigma_med','sigma_star_med','t','GP_error_med');
%        saveas(gcf,'quartic_function.fig');
%     end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schwefel?s Problem 1.2
function val = f4(x)
    val = 0;
    for i = 1:1:length(x)
        temp = 0;
        for j = 1:1:i
            temp = temp + x(j);
        end
        val = val + temp^2;
    end
end
% quartic function
function val = f5(x)
    beta = 1;
    val = 0;
    for i = 1:1:length(x)-1
        val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
    end
end
