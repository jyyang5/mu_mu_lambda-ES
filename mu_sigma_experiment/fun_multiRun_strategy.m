% Test parameter E.G. mu, lambda for mml-ES
% run several replicates and take the median as result
function val = fun_multiRun_strategy(f,mu,lambda,NUM_OF_RUNS)
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
% initialization
% f = @f5;
% name = 10;
% NUM_OF_RUNS = 20;
% mu = 4;
% lambda = 15;

n = 10;
sigma0 = 1;
NUM_OF_ITERATIONS = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% (mu/mu,lambda)-ES with GP
T_array = zeros(NUM_OF_RUNS,1);                    % # of iterations for the stop creteria
f_x_matrix = zeros(NUM_OF_RUNS,10000);             % store all fx
sigma_matrix = zeros(NUM_OF_RUNS,10000);           % store all sigma
sigma_final = zeros(NUM_OF_RUNS,1);                % the last sigma
x_final = zeros(NUM_OF_RUNS,n);                    % last x
sigma_star_matrix = zeros(NUM_OF_RUNS,10000);      % normalized step size           
GP_error_matrix = zeros(NUM_OF_RUNS,10000);        % relative error of GP            


for i = 1:NUM_OF_RUNS
    
    x0 = randn(n,mu);
    % (mu/mu,lambda)-ES with GP
    a = mml_GP(f,x0,sigma0,lambda,NUM_OF_ITERATIONS);
%     a = withGP(f,x1,sigma0,NUM_OF_ITERATIONS);

    T_array(i) = cell2mat(a(1));
    f_x_matrix(i,:) = cell2mat(a(6));
    sigma_matrix(i,:) = cell2mat(a(4));         % last sigma
    temp = cell2mat(a(4));
    sigma_final(i) = temp(T_array(i));
    temp = cell2mat(a(2));
    x_final(i,:) = temp';
    sigma_star_matrix(i,:) = cell2mat(a(9));
    GP_error_matrix(i,:) = cell2mat(a(8));
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take median
    T_med = median(T_array);
    T_max = max(T_array);
    T_min = min(T_array);
    f_x_med = median(f_x_matrix);
    sigma_med = median(sigma_matrix);
    sigma_final_med = median(sigma_final);
    x_final_med = median(x_final);
    sigma_star_med = median(sigma_star_matrix);
    GP_error_med = median(GP_error_matrix);
    
t = T_med+40;
    
    
    val = {t};
end
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
