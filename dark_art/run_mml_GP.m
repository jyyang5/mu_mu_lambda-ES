% runner for mml_GP_CSA
% plot the objective function value, step size, relative error
%function val = mml(f,x0,sigma_star,sigma_ep_star,lambda,sigma0,NUM_OF_ITERATIONS)
% initialization
% f:                  objective function value
% x0:                 mu initial point
% mu:                 population size
% sigma0:             initial muttaion strength
% NUM_OF_ITERATIONS:  number of maximum iterations



% OPTIMAL:            global optima
% TARGET_DISTANCE:    target distance to global optima
% example input:      fun = @(x) x' * x
%                     noGP(fun, randn(10,1),1,1000) 

% Test functions
f1 = @(x) (x'*x)^(1/2);  % linear sphere
f2 = @(x) (x'*x);        % quadratic sphere
f3 = @(x) (x'*x)^(3/2);  % cubic sphere


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = 5;
n = 10;

lambda = 40;
mu = ceil(lambda/4);
sigma0 = 1;
NUM_OF_ITERATIONS = 2000;
TRAINING_SIZE = 40;
LENGTH_SCALE = 4;
DECREASE_FACTOR = 0.95;
x0 = randn(n,mu);
% 64 worked for quadratic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if wanna compare and plot different strategy with legends cmp_legend != 0
% sigle plot cmp_legend = 1
cmp_legend =0; 

a = mml_GP_CSA_Niko(fname,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,DECREASE_FACTOR);
t = cell2mat(a(1));
centroid = cell2mat(a(2));
f_centroid = cell2mat(a(3));
sigma_array = cell2mat(a(4));
T = cell2mat(a(5));
fcentroid_array = cell2mat(a(6));
convergence_rate = cell2mat(a(7));
error_array = cell2mat(a(8));
sigma_satr_array = cell2mat(a(9));
success_rate = cell2mat(a(10));
delta_array = cell2mat(a(11));
% emergency_rate = cell2mat(a(12));

% p_array = cell2mat(a(12));
% b = mml(f1,x0,sigma0,lambda,NUM_OF_ITERATIONS);
% 
% t1 = cell2mat(b(1));
% centroid1 = cell2mat(b(2));
% f_centroid1 = cell2mat(b(3));
% sigma_array1 = cell2mat(b(4));
% fcentroid_array1 = cell2mat(b(6));
% convergence_rate1 = cell2mat(b(7));
% fep_centroid_array1 = cell2mat(b(8));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(fname==1)
    disp('linear')
elseif(fname==2)
    disp('quadratic')
elseif(fname==3)
    disp('cubic')
elseif(fname==4)
    disp('Schwefel')
elseif(fname==5)
    disp('quartic')
end

disp("number of iterations");
disp(t);
% disp(t1);

disp("convergence rate");
disp(convergence_rate);

disp("success rate");
disp(success_rate);

disp("number of objective function calls");
disp(T);

% disp("Emergency rate");
% disp(emergency_rate);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot sigma, f(x), sigmaStar
    figure(10);hold on;
    % f(x)
    subplot(1,4,1)          
    semilogy(1:t, fcentroid_array(1:t));hold on;           % mml with GP
    xlabel('number of iterations','fontsize',15);
    ylabel('log(f(x))','fontsize',15);
    set(gca,'yscale','log')
    title('f(centroid)','fontsize',20);
    % sigma    
    subplot(1,4,2)
    semilogy(1:t, sigma_array(1:t));hold on;               % mml with GP
    xlabel('number of iterations','fontsize',15);
    ylabel('log(\sigma)','fontsize',15);
    set(gca,'yscale','log')
    title('step size \sigma','fontsize',20);
    % GP error    
    subplot(1,4,3)
    plot(1:t,error_array(1:t));hold on;
    xlabel('number of iterations','fontsize',15);
    ylabel('relative error','fontsize',15);
    set(gca,'yscale','log')
    title('relative error','fontsize',20);
    xlabel('number of iterations','fontsize',15);
    % sigmaStar
    subplot(1,4,4)
    plot(1:t,sigma_satr_array(1:t));hold on;
    xlabel('number of iterations','fontsize',15);
    ylabel('normalized step size \sigma*','fontsize',15);
    set(gca,'yscale','log')
    title('normalized step size \sigma*','fontsize',20);
    xlabel('number of iterations','fontsize',15);
    
    figure(11)
    
    histogram(delta_array(1:t),'Normalization','probability');hold on;
    title('Delta pdf','FontSize',20);
    
    
    % plot pdf curve (green)[.2 .71 .3]
    data = delta_array(1:t);
    h = histogram(delta_array(1:t),'Normalization','probability','facecolor',[.2 .71 .3]);hold on;
    c = h.BinCounts;
    w = h.BinWidth;
    limit = h.BinLimits;
    plot(limit(1)+w/2:w:limit(2)+w/2,c/sum(c),'color',[.2 .71 .3]);hold on;
    
    % plot pdf curve (dark blue)[.25 .55 .79]
    data = delta_array(1:t);
    h = histogram(delta_array(1:t),'Normalization','probability','facecolor',[.25 .55 .79]);hold on;
    c = h.BinCounts;
    w = h.BinWidth;
    limit = h.BinLimits;
    plot(limit(1)+w/2:w:limit(2)+w/2,c/sum(c),'color',[.9 .1 .14]);hold on;
    
    % dark red [.9 .1 .14]
    data = delta_array(1:t);
    h = histogram(delta_array(1:t),'Normalization','probability','facecolor',[.9 .1 .14]);hold on;
    c = h.BinCounts;
    w = h.BinWidth;
    limit = h.BinLimits;
    plot(limit(1)+w/2:w:limit(2)+w/2,c/sum(c),'color',[.9 .1 .14]);hold on;
    
    % orange [1 0.5 0]
    data = delta_array(1:t);
    h = histogram(delta_array(1:t),'Normalization','probability','facecolor',[1 0.5 0]);hold on;
    c = h.BinCounts;
    w = h.BinWidth;
    limit = h.BinLimits;
    plot(limit(1)+w/2:w:limit(2)+w/2,c/sum(c),'color',[1 0.5 0]);hold on;
    
    
mu = mean(data);
sigma = std(data);

f = exp(-(x-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schwefel's Problem 1.2
function val = f4(x)
    val = 0;
    for i = 1:1:length(x)
        val = val + sum(x(1:i))^2;
    end
end
% quartic function
function val = f5(x)
    beta = 1;
    val = 0;
    for i = 1:1:n-1
        val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
    end
end