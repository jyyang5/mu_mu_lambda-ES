% Refactored simplify sigma_star and lambda array input
% plot objective function value over funCalls
% For a combination of sigma* and lambda
% 
% save objective function calls data & convergence plot
%
% difficulty: 
%            first several iterations lambda objective function calls
%            add one legned when doing a plot
%            use a 4-dim matrix to save the data for different lambda and
%            plot replaced in the first loop
%            enable adding new plot using different TRAINING_FACTOR to
%            original
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = fun_mergedGraph_fx_funCalls(f,name,NUM_OF_RUNS,sigma_star_array,lambda_array,TRAINING_FACTOR)
%Input:
%    f:                 objective function 
%    name:              a number specify f  
%    NUM_OF_RUNS        # of runs to average
%    sigma_star_array   an array of sigma_star's
%    lambda_array       an array of lambda's
%    TRAINING_FACTOR    number of iterations needed to build the GP model 
%Return:
%    iteration number for [mmlWithGP,mmlNoGP,1+1WithGP,1+1NoGP]  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NUM_OF_ITERATIONS = 1000;
n = 10;


[index ,SIGMA_LENGTH] = size(sigma_star_array);
[index ,LAMBDA_LENGTH] = size(lambda_array);

t_array = zeros(SIGMA_LENGTH, LAMBDA_LENGTH,NUM_OF_RUNS,1);                 % # of iterations 
T_array = zeros(SIGMA_LENGTH, LAMBDA_LENGTH,NUM_OF_RUNS,1);                 % # of objective function calls
f_x_matrix = zeros(SIGMA_LENGTH, LAMBDA_LENGTH,NUM_OF_RUNS,10000);          % store all fx
T_med = zeros(SIGMA_LENGTH,LAMBDA_LENGTH);
t_med = zeros(SIGMA_LENGTH,LAMBDA_LENGTH);
f_x_med = zeros(SIGMA_LENGTH,LAMBDA_LENGTH,10000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open the file to write data 
if name == 6
    fileID = fopen('linear_sphere_funCalls.txt','a');
    d = sprintf('The file is obtained by taking the median of %d runs of mml-ES with GP (dim(data)=%d)\n \n',NUM_OF_RUNS,n);
    fprintf(fileID,d);
    fprintf(fileID,'%10s %5s %4s %9s\n','(mu/mu,lambda)','\sigma^*','F','linear'); 
elseif name == 7
    fileID = fopen('quadratic_sphere_funCalls.txt','a');
    d = sprintf('The file is obtained by taking the median of %d runs of mml-ES with GP (dim(data)=%d)\n \n',NUM_OF_RUNS,n);
    fprintf(fileID,d);
    fprintf(fileID,'%10s %5s %5s %6s\n','(mu/mu,lambda)','\sigma^*','F','quadratic'); 
elseif name == 8
    fileID = fopen('cubic_sphere_funCalls.txt','a');
    d = sprintf('The file is obtained by taking the median of %d runs of mml-ES with GP (dim(data)=%d)\n \n',NUM_OF_RUNS,n);
    fprintf(fileID,d);
    fprintf(fileID,'%10s %5s %5s %6s\n','(mu/mu,lambda)','\sigma^*','F','cubic'); 
elseif name == 9
    fileID = fopen('Schwefel_funCalls.txt','a');
    d = sprintf('The file is obtained by taking the median of %d runs of mml-ES with GP (dim(data)=%d)\n \n',NUM_OF_RUNS,n);
    fprintf(fileID,d);
    fprintf(fileID,'%10s %5s %5s %6s\n','(mu/mu,lambda)','\sigma^*','F','Schwefel'); 
elseif name == 10
    fileID = fopen('quartic_funCalls.txt','a');
    d = sprintf('The file is obtained by taking the median of %d runs of mml-ES with GP (dim(data)=%d)\n \n',NUM_OF_RUNS,n);
    fprintf(fileID,d);
    fprintf(fileID,'%10s %5s %5s %6s\n','(mu/mu,lambda)','\sigma^*','F','quartic'); 
end
d = sprintf('------------------------------------------------------------------------\n');
fprintf(fileID,d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute statics and save
for i = 1:1:SIGMA_LENGTH
    sigma_star_temp = sigma_star_array(i);
    for j = 1:1:LAMBDA_LENGTH
        lambda_temp = lambda_array(j);
        for k = 1:NUM_OF_RUNS
            mu = ceil(0.27*lambda_temp);
            x0 = randn(n,mu);
            % (mu/mu,lambda)-ES with GP
            a = mml_sigmaStar_GP_TrainSize(f,x0,sigma_star_temp,lambda_temp,NUM_OF_ITERATIONS,TRAINING_FACTOR);
            t_array(i,j,k) = cell2mat(a(1));
            T_array(i,j,k) = cell2mat(a(5));
            f_x_matrix(i,j,k,:) = cell2mat(a(6));
        end
    end
end

T_med = median(T_array,3);
t_med = median(t_array,3);
t_max = max(t_array,[],3);
t_min = min(t_array,[],3);
f_x_med = median(f_x_matrix,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graph

figure(name);
legend('-DynamicLegend'); 
hold on;
for i = 1:1:SIGMA_LENGTH
    sigma_star_temp = sigma_star_array(i);
    
    figure(name);
    hold on;
    for j = 1:1:LAMBDA_LENGTH
        lambda_temp = lambda_array(j);
        mu = ceil(0.27*lambda_temp);
        % plot first ceil(TRAINING_SIZE/lambda_temp) iterations seperately
        t_range1 = 1:lambda_temp:(lambda_temp+1)*TRAINING_FACTOR+1;
        t_range2 = (lambda_temp+1)*TRAINING_FACTOR+2:t_max+lambda_temp*TRAINING_FACTOR;
        f_x_temp = squeeze(f_x_med(i,j,:));
        f_x_med_range1 = f_x_temp(1:TRAINING_FACTOR+1);
        f_x_med_range2 = f_x_temp(TRAINING_FACTOR+2:t_max);
        % add legend    
        d =sprintf('(%d/%d,%d) %d,%d',mu,mu,lambda_temp,sigma_star_temp,TRAINING_FACTOR);
        plot([t_range1 t_range2], [f_x_med_range1' f_x_med_range2'],'DisplayName',d);hold on;
        % write file
        d1 =sprintf('(%d/%d,%d)  \t  %d  \t  %d \t %d\n',mu,mu,lambda_temp,sigma_star_temp,TRAINING_FACTOR,T_med(i,j));
        fprintf(fileID,d1); 
    end
end
    legend('-DynamicLegend'); 
    legend('show');
    leg = legend();
    title(leg,'(\mu/\mu,\lambda) \sigma^*,F');
    hold off;
    
    xlabel('number of function evaluations','fontsize',15);
    ylabel('objective function value f(x)','fontsize',15);
    set(gca,'yscale','log','fontsize',15);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close the file to put data
fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot title seperately and save fig
if name == 6
    title('linear sphere: convergence plot ','fontsize',20);
    saveas(gcf,'linear_sphere.fig');
elseif name == 7
    title('quadratic sphere: convergence plot ','fontsize',20);
    saveas(gcf,'quadratic_sphere.fig');
elseif name == 8
    title('cubic sphere: convergence plot ','fontsize',20);
    saveas(gcf,'cubic_sphere.fig');
elseif name == 9
    title('schwefel function: convergence plot ','fontsize',20);
    saveas(gcf,'schwefel_function.fig');
elseif name == 10
    title('quartic function: convergence plot ','fontsize',20);
    saveas(gcf,'quartic_function.fig');
end
    
    val = T_med;

    
end