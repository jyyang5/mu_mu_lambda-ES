% Plot relative error for GP 
% (1+1)-ES with GP and mml-ES with GP
% use 5 test functions 
% plot is merged into one
% plot and data automatically saved

% initialization
% NUM_OF_RUNS:        # of reolicates to take dedian
% mu and lambda       initialization for (mu/mu,lambda)-ES
% n:                  dimension of the problem(data)
% f:                  objective function value
% sigma0:             initial muttaion strength
% NUM_OF_ITERATIONS:  number of maximum iterations

 

% Test functions
f1 = @(x) (x'*x)^(1/2);  % linear sphere
f2 = @(x) (x'*x);        % quadratic sphere
f3 = @(x) (x'*x)^(3/2);  % cubic sphere


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up
NUM_OF_RUNS = 2;
mu = 4;
lambda = 15;
n = 10;
sigma0 = 1;
NUM_OF_ITERATIONS = 10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lambda = 10;
% mu = 3;



% (mu/mu,lambda)-ES with GP
T_array1 = zeros(5,NUM_OF_RUNS,1);                    % # of iterations for the stop creteria
f_x_matrix1 = zeros(5,NUM_OF_RUNS,10000);             % store all fx
sigma_matrix1 = zeros(5,NUM_OF_RUNS,10000);           % store all sigma
GP_error_matrix1 = zeros(5,NUM_OF_RUNS,10000);        % relative error of GP 
T_med1 = zeros(5,1);
T_max1 = zeros(5,1);
T_min1 = zeros(5,1);
f_x_med1 = zeros(5,10000);
sigma_med1 = zeros(5,10000);
GP_error_med1 = zeros(5,10000);

% (1+1)-ES with GP
T_array2 = zeros(5,NUM_OF_RUNS,1);                    % # of iterations for the stop creteria
f_x_matrix2 = zeros(5,NUM_OF_RUNS,10000);             % store all fx
sigma_matrix2 = zeros(5,NUM_OF_RUNS,10000);           % store all sigma       
GP_error_matrix2 = zeros(5,NUM_OF_RUNS,10000);        % relative error of GP            
T_med2 = zeros(5,1);
T_max2 = zeros(5,1);
T_min2 = zeros(5,1);
f_x_med2 = zeros(5,10000);
sigma_med2 = zeros(5,10000);
GP_error_med2 = zeros(5,10000);


for i = 1:NUM_OF_RUNS
    
    x0 = randn(n,mu);
    % (mu/mu,lambda)-ES with GP
    a1 = mml_GP(f1,x0,sigma0,lambda,NUM_OF_ITERATIONS);
    a2 = mml_GP(f2,x0,sigma0,lambda,NUM_OF_ITERATIONS);
    a3 = mml_GP(f3,x0,sigma0,lambda,NUM_OF_ITERATIONS);
    a4 = mml_GP(@f4,x0,sigma0,lambda,NUM_OF_ITERATIONS);
    a5 = mml_GP(@f5,x0,sigma0,lambda,NUM_OF_ITERATIONS);
    
    T_array1(1,i) = cell2mat(a1(1));
    f_x_matrix1(1,i,:) = cell2mat(a1(6));
    sigma_matrix1(1,i,:) = cell2mat(a1(4));         
    GP_error_matrix1(1,i,:) = cell2mat(a1(8));
    
    T_array1(2,i) = cell2mat(a2(1));
    f_x_matrix1(2,i,:) = cell2mat(a2(6));
    sigma_matrix1(2,i,:) = cell2mat(a2(4));         
    GP_error_matrix1(2,i,:) = cell2mat(a2(8));
    
    T_array1(3,i) = cell2mat(a3(1));
    f_x_matrix1(3,i,:) = cell2mat(a3(6));
    sigma_matrix1(3,i,:) = cell2mat(a3(4));         
    GP_error_matrix1(3,i,:) = cell2mat(a3(8));
    
    T_array1(4,i) = cell2mat(a4(1));
    f_x_matrix1(4,i,:) = cell2mat(a4(6));
    sigma_matrix1(4,i,:) = cell2mat(a4(4));         
    GP_error_matrix1(4,i,:) = cell2mat(a4(8));
    
    T_array1(5,i) = cell2mat(a5(1));
    f_x_matrix1(5,i,:) = cell2mat(a5(6));
    sigma_matrix1(5,i,:) = cell2mat(a5(4));         
    GP_error_matrix1(5,i,:) = cell2mat(a5(8));
    
    x0 = randn(n,1);
    b1 = withGP(f1,x0,sigma0,NUM_OF_ITERATIONS);
    b2 = withGP(f2,x0,sigma0,NUM_OF_ITERATIONS);
    b3 = withGP(f3,x0,sigma0,NUM_OF_ITERATIONS);
    b4 = withGP(@f4,x0,sigma0,NUM_OF_ITERATIONS);
    b5 = withGP(@f5,x0,sigma0,NUM_OF_ITERATIONS);
    
    T_array2(1,i) = cell2mat(b1(1));
    f_x_matrix2(1,i,:) = cell2mat(b1(6));
    sigma_matrix2(1,i,:) = cell2mat(b1(4));         
    GP_error_matrix2(1,i,:) = cell2mat(b1(8));
    
    T_array2(2,i) = cell2mat(b2(1));
    f_x_matrix2(2,i,:) = cell2mat(b2(6));
    sigma_matrix2(2,i,:) = cell2mat(b2(4));         
    GP_error_matrix2(2,i,:) = cell2mat(b2(8));
    
    T_array2(3,i) = cell2mat(b3(1));
    f_x_matrix2(3,i,:) = cell2mat(b3(6));
    sigma_matrix2(3,i,:) = cell2mat(b3(4));         
    GP_error_matrix2(3,i,:) = cell2mat(b3(8));
    
    T_array2(4,i) = cell2mat(b4(1));
    f_x_matrix2(4,i,:) = cell2mat(b4(6));
    sigma_matrix2(4,i,:) = cell2mat(b4(4));         
    GP_error_matrix2(4,i,:) = cell2mat(b4(8));
    
    T_array2(5,i) = cell2mat(b5(1));
    f_x_matrix2(5,i,:) = cell2mat(b5(6));
    sigma_matrix2(5,i,:) = cell2mat(b5(4));         
    GP_error_matrix2(5,i,:) = cell2mat(b5(8));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take median

T_med1 = median(T_array1,2);
T_max1 = max(T_array1,[],2);
T_min1 = min(T_array1,[],2);
f_x_med1 = median(f_x_matrix1,2);
sigma_med1 = median(sigma_matrix1,2);
GP_error_med1 = median(GP_error_matrix1,2);

T_med2 = median(T_array2,2);
T_max2 = max(T_array2,[],2);
T_min2 = min(T_array2,[],2);
f_x_med2 = median(f_x_matrix2,2);
sigma_med2 = median(sigma_matrix2,2);
GP_error_med2 = median(GP_error_matrix2,2);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph

figure(1);
for i=1:1:5
    subplot(1,5,i)
    hold on;
    % plot first 40 iterations seperately
    temp1 = GP_error_med1(i,:);
    temp2 = GP_error_med2(i,:);
    plot(temp1(temp1~=0));
    plot(temp2(temp2~=0));
%     plot(GP_error_med1(i,1:floor(T_max1(i))));
%     plot(GP_error_med2(i,1:floor(T_max2(i))));
    set(gca,'yscale','log')
    xlabel('iterations','fontsize',12);
    legend({'mml with GP','(1+1)-ES with GP'},'fontsize',10);
    if(i==1)
        ylabel('relative model error of GP','fontsize',12);
        title('linear sphere','fontsize',15);
    elseif(i==2)
        title('quadratic sphere','fontsize',15);
    elseif(i==3)
        title('cubic sphere','fontsize',15);
    elseif(i==4)
        title('schwefel function','fontsize',15);
    elseif(i==5)
        title('quartic function','fontsize',15);
    end
    hold off;
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save data

save('relative_error_GP.mat','mu','lambda','sigma0','a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','T_array1','T_array2','GP_error_med1','GP_error_med2','f_x_matrix1','f_x_matrix2','GP_error_matrix1','GP_error_matrix2');
saveas(gcf,'relative_error_GP.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schwefel's Problem 1.2
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
