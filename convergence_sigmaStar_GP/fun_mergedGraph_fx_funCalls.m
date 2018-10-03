% Refactored simplify sigma_star and lambda array input
% plot convergence rate over sigma_star for different lambda (one plot)
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


NUM_OF_ITERATIONS = 1500;
n = 10;


[index ,SIGMA_LENGTH] = size(sigma_star_array);
[index ,LAMBDA_LENGTH] = size(lambda_array);

c_array = zeros(LAMBDA_LENGTH,SIGMA_LENGTH,NUM_OF_RUNS,1);                 % convergence rate
s_array = zeros(LAMBDA_LENGTH,SIGMA_LENGTH,NUM_OF_RUNS,1);                 % success rate
t_array = zeros(LAMBDA_LENGTH,SIGMA_LENGTH,NUM_OF_RUNS,1);                 % # of iterations 
T_array = zeros(LAMBDA_LENGTH,SIGMA_LENGTH,NUM_OF_RUNS,1);                 % # of objective function calls
f_x_matrix = zeros(LAMBDA_LENGTH, SIGMA_LENGTH,NUM_OF_RUNS,10000);          % store all fx
T_med = zeros(LAMBDA_LENGTH,SIGMA_LENGTH);
t_med = zeros(LAMBDA_LENGTH,SIGMA_LENGTH);
f_x_med = zeros(LAMBDA_LENGTH,SIGMA_LENGTH,10000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open the file to write data 
if name == 6
    fileID = fopen('linear_sphere_convergence_sigmaStar.txt','a');
    d = sprintf('The file is obtained by taking the median of %d runs of mml-ES with GP (dim(data)=%d)\n \n',NUM_OF_RUNS,n);
    fprintf(fileID,d);
    fprintf(fileID,'%10s \t','(mu/mu,lambda)'); 
elseif name == 7
    fileID = fopen('quadratic_sphere_convergence_sigmaStar.txt','a');
    d = sprintf('The file is obtained by taking the median of %d runs of mml-ES with GP (dim(data)=%d)\n \n',NUM_OF_RUNS,n);
    fprintf(fileID,d);
    fprintf(fileID,'%10s \t','(mu/mu,lambda)'); 
elseif name == 8
    fileID = fopen('cubic_sphere_convergence_sigmaStar.txt','a');
    d = sprintf('The file is obtained by taking the median of %d runs of mml-ES with GP (dim(data)=%d)\n \n',NUM_OF_RUNS,n);
    fprintf(fileID,d);
    fprintf(fileID,'%10s \t','(mu/mu,lambda)');
elseif name == 9
    fileID = fopen('Schwefel_convergence_sigmaStar.txt','a');
    d = sprintf('The file is obtained by taking the median of %d runs of mml-ES with GP (dim(data)=%d)\n \n',NUM_OF_RUNS,n);
    fprintf(fileID,d);
    fprintf(fileID,'%10s \t','(mu/mu,lambda)');
elseif name == 10
    fileID = fopen('quartic_convergence_sigmaStar.txt','a');
    d = sprintf('The file is obtained by taking the median of %d runs of mml-ES with GP (dim(data)=%d)\n \n',NUM_OF_RUNS,n);
    fprintf(fileID,d);
    fprintf(fileID,'%10s \t','(mu/mu,lambda)');
end
for i = 1:1:SIGMA_LENGTH
    d = sprintf('%.2f \t \t',sigma_star_array(i));
    fprintf(fileID,d);
end
d = sprintf('\n---------------------------------------------------------------------------------------------------------------\n');
fprintf(fileID,d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute statics and save
for i = 1:1:LAMBDA_LENGTH
    lambda_temp = lambda_array(i);
    for j = 1:1:SIGMA_LENGTH
        sigma_star_temp = sigma_star_array(j);
        for k = 1:NUM_OF_RUNS
            mu = floor(lambda_temp/4);
            x0 = randn(n,mu);
            % (mu/mu,lambda)-ES with GP
            a = mml_sigmaStarGP_bestSoFar(f,x0,sigma_star_temp,lambda_temp,NUM_OF_ITERATIONS,TRAINING_FACTOR);
            c_array(i,j,k) = cell2mat(a(7));                                % convergence rate
            s_array(i,j,k) = cell2mat(a(11));                               % success rate
            
            t_array(i,j,k) = cell2mat(a(1));
            T_array(i,j,k) = cell2mat(a(5));
            f_x_matrix(i,j,k,:) = cell2mat(a(6));
        end
    end
end

c_med = median(c_array,3);
s_med = median(s_array,3);

T_med = median(T_array,3);
% t_med = median(t_array,3);
% t_max = max(t_array,[],3);
% t_min = min(t_array,[],3);
% f_x_med = median(f_x_matrix,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graph

figure(name);
legend('-DynamicLegend'); 
hold on;
for i = 1:1:LAMBDA_LENGTH
    lambda_temp = lambda_array(i);
    mu = floor(lambda_temp/4);
    
    d1 =sprintf('(%d/%d,%d)  \t',mu,mu,lambda_temp);
    fprintf(fileID,d1); 
    % plot convergece rate of sigma_temp over sigma_star
    figure(name);
    hold on;
    d =sprintf('(%d/%d,%d)',mu,mu,lambda_temp);
    temp = c_med(i,:);
    temp = temp(temp~=0);
    [index, length] = size(temp);
    scatter(sigma_star_array(1:length),temp,'DisplayName',d); hold on;
    plot(sigma_star_array,c_med(i,:),'DisplayName',d);hold on;
    % write convergence rate for each sigma_star 
    for j = 1:1:SIGMA_LENGTH
        d1 =sprintf('%.4f  \t ',c_med(i,j));
        fprintf(fileID,d1); 
    end
    % change line 
    fprintf(fileID,'\n'); 
%     
%         sigma_star_temp = sigma_star_array(j);        
%         % plot first ceil(TRAINING_SIZE/lambda_temp) iterations seperately
%         t_range1 = 1:lambda_temp:(lambda_temp+1)*TRAINING_FACTOR+1;
%         t_range2 = (lambda_temp+1)*TRAINING_FACTOR+2:t_max+lambda_temp*TRAINING_FACTOR;
%         f_x_temp = squeeze(f_x_med(i,j,:));
%         f_x_med_range1 = f_x_temp(1:TRAINING_FACTOR+1);
%         f_x_med_range2 = f_x_temp(TRAINING_FACTOR+2:t_max);
%     end
    % add legend    
    
%     plot([t_range1 t_range2], [f_x_med_range1' f_x_med_range2'],'DisplayName',d);hold on; 
end
    legend('-DynamicLegend'); 
    legend('show');
    leg = legend();
    title(leg,'(\mu/\mu,\lambda) \sigma^*');
    hold off;
    
    xlabel('\sigma*','fontsize',15);
    ylabel('convergence rate','fontsize',15);
%     set(gca,'yscale','log','fontsize',15);

%     % scatter plots of the graph
%     figure(name);
%     hold on;
% for i = 1:1:LAMBDA_LENGTH
%     scatter(sigma_star_array,c_med(i,:)); hold on;
% end
%     hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close the file to put data
fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot title seperately and save fig
if name == 6
    title('linear sphere: convergence rate','fontsize',20);
    saveas(gcf,'linear_sphere_convergence_sigmaStar.fig');
    save('linear sphere_convergence rate','c_array','s_array','T_array','c_med','sigma_star_array','lambda_array','TRAINING_FACTOR','f');
elseif name == 7
    title('quadratic sphere: convergence rate','fontsize',20);
    saveas(gcf,'quadratic_sphere_convergence_sigmaStar.fig');
    save('quadratic sphere_convergence rate','c_array','s_array','T_array','c_med','sigma_star_array','lambda_array','TRAINING_FACTOR','f');
elseif name == 8
    title('cubic sphere: convergence rate','fontsize',20);
    saveas(gcf,'cubic_sphere_convergence_sigmaStar.fig');
    save('cubic sphere_convergence rate','c_array','s_array','T_array','c_med','sigma_star_array','lambda_array','TRAINING_FACTOR','f');
elseif name == 9
    title('schwefel function: convergence rate','fontsize',20);
    saveas(gcf,'schwefel_function_convergence_sigmaStar.fig');
elseif name == 10
    title('quartic function: convergence rate ','fontsize',20);
    saveas(gcf,'quartic_function_convergence_sigmaStar.fig');
end
    
    val = T_med;

    
end