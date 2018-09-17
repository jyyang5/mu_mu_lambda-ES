function val = fun_graph_merged(f,name,NUM_OF_RUNS)


sigma0 = 1;
NUM_OF_ITERATIONS = 10000;
lambda = 10;
mu = 3;
n = 10;


% (mu/mu,lambda)-ES with GP
T_array = zeros(NUM_OF_RUNS,1);                   % # of iterations for the stop creteria
f_x_matrix = zeros(NUM_OF_RUNS,10000);             % store all fx
sigma_matrix = zeros(NUM_OF_RUNS,10000);           % store all sigma
sigma_final = zeros(NUM_OF_RUNS,1);               % the last sigma
x_final = zeros(NUM_OF_RUNS,n);                   % last x
% (mu/mu,lambda)-ES usual 
T_array1 = zeros(NUM_OF_RUNS,1);                   % # of iterations for the stop creteria
f_x_matrix1 = zeros(NUM_OF_RUNS,10000);             % store all fx
sigma_matrix1 = zeros(NUM_OF_RUNS,10000);           % store all sigma
sigma_final1 = zeros(NUM_OF_RUNS,1);               % the last sigma
x_final1 = zeros(NUM_OF_RUNS,n);                   % last x

% (1+1)-ES with GP
T_array2 = zeros(NUM_OF_RUNS,1);                   % # of iterations for the stop creteria
f_x_matrix2 = zeros(NUM_OF_RUNS,10000);             % store all fx
sigma_matrix2 = zeros(NUM_OF_RUNS,10000);           % store all sigma
sigma_final2 = zeros(NUM_OF_RUNS,1);               % the last sigma
x_final2 = zeros(NUM_OF_RUNS,n);                   % last x

% (1+1)-ES 
T_array3 = zeros(NUM_OF_RUNS,1);                   % # of iterations for the stop creteria
f_x_matrix3 = zeros(NUM_OF_RUNS,10000);             % store all fx
sigma_matrix3 = zeros(NUM_OF_RUNS,10000);           % store all sigma
sigma_final3 = zeros(NUM_OF_RUNS,1);               % the last sigma
x_final3 = zeros(NUM_OF_RUNS,n);                   % last x


for i = 1:NUM_OF_RUNS
    x0 = randn(n,mu);
    a = mml_GP(f,x0,sigma0,lambda,NUM_OF_ITERATIONS);
    T_array(i) = cell2mat(a(1));
    %x_matrix(i,:) = cell2mat(a(2));
    f_x_matrix(i,:) = cell2mat(a(6));
    sigma_matrix(i,:) = cell2mat(a(4));
    % last sigma
    temp = cell2mat(a(4));
    sigma_final(i) = temp(T_array(i));
    temp = cell2mat(a(2));
    x_final(i,:) = temp';%temp(T_array(i),:);
    
    b = mml(f,x0,sigma0,lambda,NUM_OF_ITERATIONS);
    T_array1(i) = cell2mat(b(1));
    f_x_matrix1(i,:) = cell2mat(b(6));
    sigma_matrix1(i,:) = cell2mat(b(4));
    % last sigma
    temp = cell2mat(b(4));
    sigma_final1(i) = temp(T_array1(i));
    temp = cell2mat(b(2));
    x_final1(i,:) = temp';%temp(T_array(i),:);
    
    x0 = randn(n,1);
    c = withGP(f,x0,sigma0,NUM_OF_ITERATIONS);
    T_array2(i) = cell2mat(c(1));
    f_x_matrix2(i,:) = cell2mat(c(3));
    sigma_matrix2(i,:) = cell2mat(c(4));
    % last sigma
    temp = cell2mat(c(4));
    sigma_final2(i) = temp(T_array2(i));
    temp = cell2mat(c(2));
    x_final2(i,:) = temp';%temp(T_array(i),:);
    
    d = noGP(f,x0,sigma0,NUM_OF_ITERATIONS);
    T_array3(i) = cell2mat(d(1));
    f_x_matrix3(i,:) = cell2mat(d(3));
    sigma_matrix3(i,:) = cell2mat(d(4));
    % last sigma
    temp = cell2mat(d(4));
    sigma_final3(i) = temp(T_array3(i));
    temp = cell2mat(d(2));
    x_final3(i,:) = temp';%temp(T_array(i),:);
    
    
    
end
    T_med = median(T_array);
    T_max = max(T_array);
    T_min = min(T_array);
    f_x_med = median(f_x_matrix);
    sigma_med = median(sigma_matrix);
    sigma_final_med = median(sigma_final);
    x_final_med = median(x_final);
    
    T_med1 = median(T_array1);
    T_max1 = max(T_array1);
    T_min1 = min(T_array1);
    f_x_med1 = median(f_x_matrix1);
    sigma_med1 = median(sigma_matrix1);
    sigma_final_med1 = median(sigma_final1);
    x_final_med1 = median(x_final1);
    
    T_med2 = median(T_array2);
    T_max2 = max(T_array2);
    T_min2 = min(T_array2);
    f_x_med2 = median(f_x_matrix2);
    sigma_med2 = median(sigma_matrix2);
    sigma_final_med2 = median(sigma_final2);
    x_final_med2 = median(x_final2);
    
    T_med3 = median(T_array3);
    T_max3 = max(T_array3);
    T_min3 = min(T_array3);
    f_x_med3 = median(f_x_matrix3);
    sigma_med3 = median(sigma_matrix3);
    sigma_final_med3 = median(sigma_final3);
    x_final_med3 = median(x_final3);
    
    
    % graph
    figure(name);
    subplot(1,2,1)
   
    hold on;
    % plot first 40 iterations seperately
    t_range1 = 1:10:40;
    t_range2 = 5:T_max;
    sigma_med_range1 = sigma_med(1:4);
    sigma_med_range2 = sigma_med(5:T_max);
    
    semilogy([t_range1 t_range2], [sigma_med_range1 sigma_med_range2]);% mml with GP
    semilogy(1:10:10*T_max1, sigma_med1(1:T_max1));  % MML
    semilogy(1:T_max2, sigma_med2(1:T_max2));        % (1+1)-ES with GP
    semilogy(1:T_max3, sigma_med3(1:T_max3));        % (1+1)-ES

    hold off;
    legend('mml with GP','mml','(1+1)-ES with GP','(1+1)-ES');
    xlabel('number of function evaluations');
    ylabel('log(sigma)');
    set(gca,'yscale','log')
    title('sigma');
    
    
    
    subplot(1,2,2)
    hold on;
    semilogy(1:T_max, f_x_med(1:T_max));             % mml with GP
    semilogy(1:10:10*T_max1, f_x_med1(1:T_max1));    % mml 
    semilogy(1:T_max2, f_x_med2(1:T_max2));          % (1+1)-ES with GP
    semilogy(1:T_max3, f_x_med3(1:T_max3));          % (1+1)-ES
    
    hold off;
    legend('mml with GP','mml','(1+1)-ES with GP','(1+1)-ES');
    xlabel('number of function evaluations');
    ylabel('log( f(x) )');
    set(gca,'yscale','log')
    title('function evaluation');
    
    if name == 6
       save('T_linear_sphere.mat','T_array','T_array1','T_array2','T_array3');
    elseif name == 7
       save('T_quadratic_sphere.mat','T_array','T_array1','T_array2','T_array3');
    elseif name == 8
       save('T_cubic_sphere.mat','T_array','T_array1','T_array2','T_array3');
    elseif name == 9
       save('T_schwefel_function.mat','T_array','T_array1','T_array2','T_array3');
    elseif name == 10
       save('T_quartic_function.mat','T_array','T_array1','T_array2','T_array3');
    end
    
    
    val = {T_med,T_med1,T_med2,T_med3};
%     % final evaluated value 
%     disp('Last iteration output(noGP)');
%     d0 = sprintf('Number of function evaluations(median) of %d runs: %d',NUM_OF_RUNS, T_med);
%     d1 = sprintf("mutation strength(last function evaluation): %f", sigma_final_med);
%     disp(d0);
%     disp(d1);
%     disp(x_final_med);
%     
%     disp('Last iteration output(withGP)');
%     d2 = sprintf('Number of function evaluations(median) of %d runs: %d',NUM_OF_RUNS, T_med1);
%     d3 = sprintf("mutation strength(last function evaluation): %f", sigma_final_med1);
%     disp(d2);
%     disp(d3);
%     disp(x_final_med1);
    
    
    
%     disp('Last iteration output(Strategy)');
%     d4 = sprintf('Number of function evaluations(median) of %d runs: %d',NUM_OF_RUNS, T_med2);
%     d5 = sprintf("mutation strength(last function evaluation): %f", sigma_final_med2);
%     disp(d4);
%     disp(d5);
%     disp(x_final_med2);
    
end