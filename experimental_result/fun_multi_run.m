% Objective 
% 1 Plot
%     1.1 relative model noise to lambda, mu = ceil(lambda/4) 
%     1.2 then effect of changing theta
% 2. Save
%     Data     
% save objective function calls data & plot into files
%
% difficulty: 
%            first few trainijng iterations
%            histogram does not show data properly for sigmaStar and model
%            error
%            save file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = fun_multi_run(fname,NUM_OF_RUNS,lambda,TRAINING_SIZE,LENGTH_SCALE)
%Input:
%   fname:          an index 
%                       1 for linear
%                       2 for quadratic 
%                       3 for cubic 
%                       4 for schwefel
%                       5 for quartic
%    NUM_OF_RUNS    # of runs to average
%    lambda         0 for (1+1)-ES
%                   10 for (3/3,10)-ES
%                   20 for (5/5,20)-ES
%                   40 for (10/10,40)-ES
%    TRAINING_SIZE  GP training size
%    LENGTH_SCALE   length scale factor for GP
% 
%Return:
% 	1.t_array
%   2.sigma_matrix
%   3.T_array
%   4.f_x_matrix
%   5.convergence_rate_array
%   6.GP_error_matrix
%   7.sigma_star_matrix
%   8.success_rate_array
%   9.delta_matrix
   
% iteration number for [mmlWithGP,mmlNoGP,1+1WithGP,1+1NoGP]    
sigma0 = 1;
NUM_OF_ITERATIONS = 2000;
% lambda = 10;
% mu = 3;
n = 10;

subplot_ROW = 5;
subplot_COL = 5; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Def of variables

% (mu/mu,lambda)-ES with GP
t_array = zeros(NUM_OF_RUNS,1);                    % # of iterations for the stop creteria
sigma_matrix = zeros(NUM_OF_RUNS,10000);           % store all sigma
T_array = zeros(NUM_OF_RUNS,1);                    % # of objective function evaluations for the stop creteria
f_x_matrix = zeros(NUM_OF_RUNS,10000);             % store all fx
convergence_rate_array = zeros(NUM_OF_RUNS,1);     % convergence rate 
GP_error_matrix = zeros(NUM_OF_RUNS,10000);        % store similar to noise-to-signal ratio
sigma_star_matrix = zeros(NUM_OF_RUNS,10000);      % normalized step size 
success_rate_array = zeros(NUM_OF_RUNS,1);         % success rate 
delta_matrix = zeros(NUM_OF_RUNS,10000);           % each [i,j] stores a delta array 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replicates for NUM_OF_RUNS

mu = ceil(lambda/4);
for i = 1:NUM_OF_RUNS

    if(lambda==0)
        x0 = randn(n,1);
        a = withGP(fname,x0,sigma0,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE);
        t_array(i) = cell2mat(a(5));
    else 
        x0 = randn(n,mu);
        a = mml_GP_CSA_Niko(fname,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE);
        t_array(i) = cell2mat(a(1));
    end
        
    
    
    sigma_matrix(i,:) = cell2mat(a(4));
    T_array(i) = cell2mat(a(5));
    f_x_matrix(i,:) = cell2mat(a(6));
    convergence_rate_array(i) = cell2mat(a(7));
    GP_error_matrix(i,:) = cell2mat(a(8));
    sigma_star_matrix(i,:) = cell2mat(a(9));
    success_rate_array(i) = cell2mat(a(10));
    delta_matrix(i,:) = cell2mat(a(11));

    % counter
    fprintf('# of runs %d',i)
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Save data
%     if fname == 1
%         f = @(x) (x'*x)^(1/2);
%         save('linear_over_lambda.mat','fname','f','fname','NUM_OF_RUNS','lambda_array',...
%             'LENGTH_SCALE','TRAINING_SIZE',...
%         't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
%             'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
%             'delta_matrix');
%     elseif fname == 2
%         f = @(x) (x'*x);
%         save('quadratic_over_lambda.mat','fname','f','fname','NUM_OF_RUNS','lambda_array',...
%             'LENGTH_SCALE','TRAINING_SIZE',...
%         't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
%             'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
%             'delta_matrix');
%     elseif fname == 3
%         f = @(x) (x'*x)^(3/2);
%         save('cubic_over_lambda.mat','fname','f','fname','NUM_OF_RUNS','lambda_array',...
%             'LENGTH_SCALE','TRAINING_SIZE',...
%         't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
%             'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
%             'delta_matrix');
%     elseif fname == 4
%         f = @f4;
%         save('schwefel_over_lambda.mat','fname','f','fname','NUM_OF_RUNS','lambda_array',...
%             'LENGTH_SCALE','TRAINING_SIZE',...
%         't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
%             'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
%             'delta_matrix');  
%     elseif fname == 5
%         f = @f5;
%         save('quartic_over_lambda.mat','f','fname','f','fname','NUM_OF_RUNS','lambda_array',...
%             'LENGTH_SCALE','TRAINING_SIZE',...
%         't_array','sigma_matrix','T_array','f_x_matrix','convergence_rate_array',...
%             'GP_error_matrix','GP_error_matrix','sigma_star_matrix','success_rate_array',...
%             'delta_matrix');       
%     end
% 
% end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take median

t_med = median(t_array);
t_max = max(t_array);
t_min = min(t_array);
sigma_matrix_med = squeeze(median(sigma_matrix));
T_med = median(T_array);
f_x_med = squeeze(median(f_x_matrix));
convergence_med = median(convergence_rate_array);
GP_error_matrix_med = squeeze(median(GP_error_matrix));
sigma_star_matrix_med = squeeze(median(sigma_star_matrix));
success_med = median(success_rate_array);
delta_matrix_med = squeeze(median(delta_matrix));
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graph

% For f(x), sigma, sigmaStar, GP_error over objective function evaluations
figure(100);
legend('-DynamicLegend'); 
hold on;
mu = ceil(lambda/4);
if lambda==0
    d = sprintf('(1+1)-ES');
else
    d = sprintf('(%d/%d,%d)-ES',mu,mu,lambda);
end


% 1.objective function [row 1]
subplot(subplot_ROW,subplot_COL,fname);
h = histogram(T_array,'Normalization','probability','DisplayName',d);hold on;
h.BinWidth = 20;
if(lambda==40)
    legend({'(1+1)-ES','(3/3,10)-ES','(5/5,20)-ES','(10/10,40)-ES'})
    if(fname == 1)
        d3 =sprintf('linear sphere');
        title(d3,'fontsize',15);
    elseif(fname == 2)
        d3 =sprintf('quadratic sphere');
        title(d3,'fontsize',15);
    elseif(fname == 3)
        d3 =sprintf('cubic sphere');
        title(d3,'fontsize',15);
    elseif(fname == 4)
        d3 =sprintf('Schwefel;s function');
        title(d3,'fontsize',15);
    elseif(fname == 5)
        d3 =sprintf('quartic function');
        title(d3,'fontsize',15); 
    end
end
% if(fname==1)
%     ylabel('normalized pdf','FontSize',15);%
% end
xlabel('objective function calls','FontSize',15); 


% 2.objective function [row 2]
subplot(subplot_ROW,subplot_COL,fname+5);
if lambda==0
    t_start = TRAINING_SIZE+1;
    plot(1:1:t_med,f_x_med(1:T_med),'DisplayName',d);hold on;
else 
    t_start = ceil(TRAINING_SIZE/lambda);
    fx_range1 = f_x_med(1:t_start);
    fx_range2 = f_x_med(t_start+1:t_med);
    t_range1 = 1:lambda:lambda*t_start;
    t_range2 = lambda*t_start+1:lambda*t_start+length(fx_range2);
    plot([t_range1 t_range2], [fx_range1 fx_range2],'DisplayName',d);hold on;% mml with GP
end
if(fname==1)
    ylabel('objective function value','FontSize',15);%
end
xlabel('objective function calls','FontSize',15); 
set(gca, 'YScale', 'log');
legend('-DynamicLegend'); 
legend('show');

% 3.GP error [row 3]
subplot(subplot_ROW,subplot_COL,fname+10);
if lambda==0
    plot(t_start:1:T_med,GP_error_matrix_med(t_start:T_med),'DisplayName',d);hold on;
    d1 = sprintf('(%d/%d,%d)-ES[smoothed]',mu,mu,lambda);
    plot(t_start:1:T_med,smoothdata(GP_error_matrix_med(t_start:T_med),'gaussian',40),'DisplayName',d1,'LineWidth',2);hold on;
else 
    GP_error_range1 = GP_error_matrix_med(1:t_start);
    GP_error_range2 = GP_error_matrix_med(t_start+1:t_med);
    plot([t_range1 t_range2], [GP_error_range1 GP_error_range2],'DisplayName',d);hold on;
    % GP smoothed
    smoothed_GP_range1 = smoothdata(GP_error_range1,'gaussian',40);
    smoothed_GP_range2 = smoothdata(GP_error_range2,'gaussian',40);
    d1 = sprintf('(%d/%d,%d)-ES[smoothed]',mu,mu,lambda);
    plot([t_range1 t_range2],[smoothed_GP_range1 smoothed_GP_range2],'DisplayName',d1,'LineWidth',2);hold on;
end
if(fname==1)
    ylabel('relative model error','FontSize',15);%
end
xlabel('objective function evaluations','FontSize',15); 
set(gca, 'YScale', 'log');

legend('-DynamicLegend'); 
legend('show');
% title('Logarithmic relative model error','FontSize',20);


% 4.normalized step size [row 4](only makes sense for sphere functions)
if(fname<4)
    subplot(subplot_ROW,subplot_COL,fname+15);
    if lambda==0
        plot(1:1:T_med,sigma_star_matrix_med(1:T_med),'DisplayName',d);hold on;
    else 
        sigma_star_range1 = sigma_star_matrix_med(1:t_start);
        sigma_star_range2 = sigma_star_matrix_med(t_start+1:t_med);
        plot([t_range1 t_range2], [sigma_star_range1 sigma_star_range2],'DisplayName',d);hold on;
    end
    if(fname==1)
        ylabel('normalized step size \sigma*','FontSize',15);%
    end
    xlabel('objective function evaluations','FontSize',15); 
    set(gca, 'YScale', 'log');
%     title('normalized step size \sigma*','FontSize',20);

    legend('-DynamicLegend'); 
    legend('show');
end

% 5. histogram success rate (good step size)
subplot(subplot_ROW,subplot_COL,fname+20);
% No bar colour, bold binwidth
if(T_med == 5000) % early stopping
    h2 = histogram(nonzeros(delta_matrix_med),'Normalization','probability','DisplayName',d,'FaceColor','none','LineWidth',2);hold on;
else
    
    h2 = histogram(nonzeros(delta_matrix_med),'Normalization','probability','DisplayName',d,'FaceColor','none','LineWidth',2);hold on;
end
if(lambda==0)
    h2.EdgeColor= [0  0.4470 0.7410];
    h2.BinEdges=h2.BinEdges-0.015;
elseif(lambda == 10)
    h2.EdgeColor= [0.8500  0.3250  0.0980];
    h2.BinEdges=h2.BinEdges-0.005;
elseif(lambda==20)
    h2.EdgeColor= [0.9290  0.6940  0.1250];
    h2.BinEdges=h2.BinEdges+0.005;
elseif(lambda==40)
    h2.EdgeColor= [0.4940  0.1840  0.5560];
    h2.BinEdges=h2.BinEdges+0.015;
end
h2.LineWidth=1.2;
h2.BinWidth = 0.2;
if(fname==1)
    ylabel('Probability','FontSize',15);%
end
xlabel('Normalized fitness gain','FontSize',15); 
% set(gca, 'YScale', 'log');
% title('step size \sigma','FontSize',20);
legend('-DynamicLegend'); 
legend('show');

saveas(gcf,'merged_plot.fig');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram normalized fitGain [SEPERATE]
% figure(101);
% % 5. histogram success rate (good step size)
% 
% subplot(1,subplot_COL,fname);
% % No bar colour, bold binwidth
% if(T_med == 5000) % early stopping
%     h2 = histogram(nonzeros(delta_matrix_med),'Normalization','probability','DisplayName',d,'FaceColor','none','LineWidth',2);hold on;
% else
%     
%     h2 = histogram(nonzeros(delta_matrix_med),'Normalization','probability','DisplayName',d,'FaceColor','none','LineWidth',2);hold on;
% end
% if(lambda==0)
%     h2.EdgeColor= [0  0.4470 0.7410];
%     h2.BinEdges=h2.BinEdges-0.0125;
% elseif(lambda == 10)
%     h2.EdgeColor= [0.8500  0.3250  0.0980];
%     h2.BinEdges=h2.BinEdges;
% elseif(lambda==20)
%     h2.EdgeColor= [0.9290  0.6940  0.1250];
%     h2.BinEdges=h2.BinEdges+0.0125;
% elseif(lambda==40)
%     h2.EdgeColor= [0.4940  0.1840  0.5560];
%     h2.BinEdges=h2.BinEdges+0.025;
% end
% h2.LineWidth=1.2;
% h2.BinWidth = 0.2;
% if(fname==1)
%     ylabel('Probability','FontSize',15);%
% end
% xlabel('Normalized fitness gain','FontSize',15); 
% % set(gca, 'YScale', 'log');
% % title('step size \sigma','FontSize',20);
% legend('-DynamicLegend'); 
% legend('show');
% saveas(gcf,'fitGain_pdf.fig');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Success rate good step size
figure(102);
% 5. histogram success rate (good step size)
subplot(1,subplot_COL,fname);
% No bar colour, bold binwidth
if(T_med == 5000) % early stopping
    h2 = histogram(success_rate_array,'Normalization','probability','DisplayName',d);hold on;
else
    h2 = histogram(success_rate_array,'Normalization','probability','DisplayName',d);hold on;
end
% if(lambda==0)
%     h2.BinEdges=h.BinEdges-0.0125;
% elseif(lambda == 10)
%     h2.BinEdges=h.BinEdges;
% elseif(lambda==20)
%     h2.BinEdges=h.BinEdges+0.0125;
% elseif(lambda==40)
%     h2.BinEdges=h.BinEdges+0.025;
% end
% h2.LineWidth=1.5;
h2.BinWidth = 0.05;
if(fname==1)
    ylabel('Probability','FontSize',15);%
end
xlabel('Prob good step size','FontSize',15); 
% set(gca, 'YScale', 'log');
% title('step size \sigma','FontSize',20);
legend('-DynamicLegend'); 
legend('show');
saveas(gcf,'prob_good_stepSize.fig');

    
%     % svae plot    
%     if fname == 1
%         saveas(gcf,'linear_sphere_funCall.fig');
%     elseif fname == 2
%         saveas(gcf,'quadratic_sphere_funCall.fig');
%     elseif fname == 3
%         saveas(gcf,'cubic_sphere_funCall.fig');
%     elseif fname == 4 
%         saveas(gcf,'schwefel_function_funCall.fig');
%     elseif fname == 5
%         saveas(gcf,'quartic_function_funCall.fig');
%     end

% svae plot  
  


val = {t_array,sigma_matrix,T_array,f_x_matrix,convergence_rate_array,GP_error_matrix,sigma_star_matrix,success_rate_array,delta_matrix};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each lambda plot pmf fitness gain/iteration hitogram
% assume 3 rows 
%     numOfCol = ceil(LAMBDA_LENGTH/3);
%     figure(10+i+fname);
%     for i = 1:1:LAMBDA_LENGTH     
%         subplot(3,numOfCol,i);
%         lambda = lambda_array(i);
%         mu = ceil(lambda/4);
%     
%         delta_array = delta_matrix_med(i,1:t_med(i));
%         histogram(delta_array(1:t_med(i)),'Normalization','probability');
%         xlabel('number of iterations','fontsize',15);
%         ylabel('prob','fontsize',15);
%         d =sprintf('Normalized fitGain pdf (%d/%d,%d)',mu,mu,lambda);
%         title(d,'fontsize',20);
%    
%     end    
%  
%     if fname == 1
%         d =sprintf('hist_linear_%d.fig',lambda);
%     elseif fname == 2
%         d =sprintf('hist_quadratic_%d.fig',lambda);
%     elseif fname == 3
%         d =sprintf('hist_cubic_%d.fig',lambda);
%     elseif fname == 4 
%         d =sprintf('hist_schwefel_%d.fig',lambda);
%     elseif fname == 5
%         d =sprintf('hist_quartic_%d.fig',lambda);
%         saveas(gcf,'hist_quartic_function.fig');
%     end
%     saveas(gcf,d);
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot by taking the mean of relative model error in median run
%     mean_GP_error_array = zeros(LAMBDA_LENGTH,1);
%     for i = 1:1:LAMBDA_LENGTH   
%         mean_GP_error_array(i) = mean(GP_error_matrix_med(i,1:t_med(i)));
%     end
%     figure(50+fname)
%     plot(lambda_array,mean_GP_error_array);
%     ylabel('mean of relative model error in median run','FontSize',15);%
%     xlabel('\lambda','FontSize',15); 
%     set(gca, 'YScale', 'log');
%     title('Relative model error (mean)','FontSize',20);
%     
%     if fname == 1
%         d3 = sprintf('linear_modelError_lambda.fig');
%     elseif fname == 2
%         d3 = sprintf('quadratic_modelError_lambda.fig');
%     elseif fname == 3
%         d3 = sprintf('cubic_modelError_lambda.fig');
%     elseif fname == 4 
%         d3 = sprintf('schwefel_function_funCall.fig');
%     elseif fname == 5
%         d3 = sprintf('quartic_function_funCall.fig');
%     end
%     saveas(gcf,d3);
%     
%     
%     val = {T_med,convergence_med,success_med,GP_error_matrix_med,delta_matrix_med};
% 
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end