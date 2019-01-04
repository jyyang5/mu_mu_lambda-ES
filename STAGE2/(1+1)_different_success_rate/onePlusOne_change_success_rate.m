% Using different success rate to update step size [(1+1)-ES]
% function evaluation using TRUE OBJECTIVE FUNCTION
% In each iteration only centroid is evaluated use true objective Function

function val = onePlusOne_change_success_rate(fname,x0,sigma0,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,SUCCESS_RATE)
% initialization
% fname:              an index 
%                       1 for linear
%                       2 for quadratic 
%                       3 for cubic 
%                       4 for schwefel
%                       5 for quartic
% x0:                 mu initial point size [n, mu]
% sigma0:             initial step size
% lambda:             # of offsprings genenerated in each itertaion  
% mu:                 parent size
% NUM_OF_ITERATIONS:  number of maximum iterations
% TRAINING_SIZE:      surrogate training size
% LENGTH_SCALE:       theta in GP
% SUCCESS_RATE:       success rate used to update step size

% Return 
% 1.t:                  # of objective function calls                    
% 2.x_array(:,t):       last x
% 3.fx:                 last objective function value
% 4.sigma:              simage arrary over # of objective function calls  
% 5.T:                  # of objective function calls
% 6.f_x:                objective function values for parents
% 7.convergence_rate:   rate of convergence
% 8.-1:                 no GP error
% 9.sigma_star_array:   normalized step size
% 10.success_rate       success_rate
% 11.delta_array        normalized delta  (f_centroid(t)-f_centroid(t-1))/factor
%                       linear: factor = R
%                       quadratic: factor = 2*R^2
%                       cubic: factor = 3*R^3
%                       where R=dist(centroid)

% OPTIMAL:            global optima
% example input:      f = @(x) x' * x
%                     mml(f,randn(n,mu),1,40,10,1,4000)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test functions
f1 = @(x) (x'*x)^(1/2);  % linear sphere
f2 = @(x) (x'*x);        % quadratic sphere
f3 = @(x) (x'*x)^(3/2);  % cubic sphere
if(fname==1)
    f=f1;
elseif(fname==2)
    f=f2;
elseif(fname==3)
    f=f3;
elseif(fname==4)
    f=@f4;
elseif(fname==5)
    f=@f5;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization
[n,m] = size(x0);                         % dim of the data
xTrain = zeros(n, 10000);                 % parent solution with dim n               
sigma = sigma0;                           % mutation strength(temp)
% TRAINING_SIZE = 40;

% for graphing
x_array = zeros(n,10000);                 % store all best candiate solution
f_x = zeros(1,10000);                     % store all objective function val of candaite solution
x_array(:,1) = x0;
f_x(1) = f(x0);
sigma_array = zeros(1,10000);
sigma_star_array = zeros(1,10000);        % store normalized step size
fep_x = zeros(1,10000);                   % store GP estimate of parent
GP_error = zeros(1,10000);                % relative error for GP
GP_error_final = zeros(1);
delta_array = zeros(1,10000);
error_array = zeros(1,10000);

D = sqrt(n+1);

x = x0;%randn(n, 1);                      % best candiate solution dim:n
fx = f(x);                                % function value
xTrain(:, 1) = x;                         % parent solutions
fTrain(1) = fx;                           % vector: value of parent solutions

t = 1;                                    % # of iterations
T = 1;                                    % # of distinct parent solution 

% initial sigma_satr 
dist = norm(x);
sigma_star = sigma*n/dist;
sigma_star_array(T) = sigma_star;
sigma_array(1) = sigma;



while((T < NUM_OF_ITERATIONS) && fx > 10^(-8))
    % early stopping 
    if(fx > 500)
        % if diverge -> convergence rate = 0 success rate = 0
        success_rate = 0;
        convergence_rate = 0;
%         val = {9999,mean(x0, 2),9999,sigma_array, 9999, fcentroid_array,-1,error_array,sigma_star_array,success_rate,delta_array}; 
        val = {t,x_array,fx,sigma_array, T, f_x,convergence_rate,error_array,sigma_star_array,success_rate,delta_array};

        return 
    end

    
    % offspring_generation
    f_x(T)=fx;
    y = x + sigma*randn(n,1);%temp;
    fy = f(y);
    T = T+1;
    
    if(fy < fx)
        x = y;
        fx = fy;
        sigma = sigma*exp((1-SUCCESS_RATE)/D);
    else
        sigma = sigma*exp(-SUCCESS_RATE/D);
    end
    
    
%     % update GP iff. there is 40 candiate solutions
%     if T > TRAINING_SIZE
%         theta = sigma*LENGTH_SCALE*sqrt(n);                                 % NOtE: need to be updated every time
%         fy_ep = gp(xTrain(:, T-40:T-1), fTrain(T-40:T-1), y, theta);        % fitness of offspring(use GP)
%         % update GP estimate for parent
%         fep_x(T) = gp(xTrain(:, T-40:T-1), fTrain(T-40:T-1), x, theta);
%     end
%     
%     % update mutation & assign new offspring
%     % if GP already built compare
%     if(T > TRAINING_SIZE && fy_ep >= fx)             % bad offspring
%         sigma = sigma * exp(-c1/D);
%         % relative GP error
%         y_temp = f(y);
%         
%     % GP not built || offspring inferior
%     else
%         fy = f(y);                        % fitness of offspring(use true objective fn)
%         
%         xTrain(:, T) = y;                
%         fTrain(T) = fy;
%         T = T + 1;
%         % GP error
%         if(T>TRAINING_SIZE+1)
%             GP_error(T) = abs(fy_ep-fy)./abs(fy-f(x));                      % relative error of GP |f(y)-fep(y)|/ |f(y)-f(x)| after GP built
%         end
%         if(fy >= fx)                      % bad offspring                      
%             sigma = sigma * exp(-c2/D);   % reduce step size
%         else
%             x = y;
%             fx = fy;
%             sigma = sigma * exp(c3/D);   % increase step size
%         end
%         x_array(:,T) = x;
%         f_x(T) = fx;
%         sigma_funEva_array(T) = sigma;
%         
%         
%     end 
    % new iteration     
    t = t + 1;
    
    sigma_array(T) = sigma;
    
    % store last normalized step size
    sigma_star_array(T) = sigma*n/norm(x);
    
    
end
 

    % convergence rate (overall)
    t_start = ceil(TRAINING_SIZE);
    if(fname==1)
        delta_array(2:t) = -(f_x(2:t)-f_x(1:t-1))./vecnorm(f_x(:,2:t),2,1);
        convergence_rate = -n*sum(log(f_x(t_start+2:t)./f_x(t_start+1:t-1)))/length(f_x(t_start+1:t-1));
    elseif(fname==2)
        delta_array(2:t) = -(f_x(2:t)-f_x(1:t-1))./(vecnorm(f_x(:,2:t),2,1)).^2/2;
        convergence_rate = -n/2*sum(log(f_x(t_start+2:t)./f_x(t_start+1:t-1)))/length(f_x(t_start+1:t-1));
    elseif(fname==3)
        delta_array(2:t) = -(f_x(2:t)-f_x(1:t-1))./(vecnorm(f_x(:,2:t),2,1)).^3/3;        
        convergence_rate = -n/3*sum(log(f_x(t_start+2:t)./f_x(t_start+1:t-1)))/length(f_x(t_start+1:t-1));
    elseif(fname==4)
        delta_array(2:t) = -(f_x(2:t)-f_x(1:t-1))./vecnorm(f_x(:,2:t),2,1);
        convergence_rate = -n/2*sum(log(f_x(t_start+2:t)./f_x(t_start+1:t-1)))/length(f_x(t_start+1:t-1));
    elseif(fname==5)
        delta_array(2:t) = -(f_x(2:t)-f_x(1:t-1))./vecnorm(f_x(:,2:t),2,1);
    	convergence_rate = -n/2*sum(log(f_x(t_start+2:t)./f_x(t_start+1:t-1)))/length(f_x(t_start+1:t-1));
    end
    % success rate
    success_rate = sum(f_x(t_start:T-1)>f_x(t_start+1:T))/length(f_x(t_start:T-1));
    delta_array = -delta_array;
    val = {t,x_array,fx,sigma_array, T, f_x,convergence_rate,error_array,sigma_star_array,success_rate,delta_array};

end



function fTest = gp(xTrain, fTrain, xTest, theta)
% input: 
%       xTrain(40 training pts)
%       fTrain(true objective function value)
%       xTest(1 test pt)   
%       theta length scale of GP    
% return: the prediction of input test data

    [n, m] = size(xTrain);                                       % m:  # of training data

    delta = vecnorm(repmat(xTrain, 1, m)-repelem(xTrain, 1, m)); %|x_ij = train_i-train_j| 
    K = reshape(exp(-delta.^2/theta^2/2), m , m);                % K

    deltas = vecnorm(xTrain-repelem(xTest, 1, m));               %|x_ij = train_i-test_j| euclidean distance
    Ks = exp(-(deltas/theta).^2/2)';                             % K_star             

%     deltass = vecnorm(repmat(xTest, 1, m)-repelem(xTest, 1, m));
%     Kss = reshape((exp(-deltass.^2/theta^2/2)), m , m);
%     Kinv = inv(K);       
%       mu = fTrain(40);
     mu = min(fTrain);                                            % estimated mean of GP
    
%     fTest = min(fTrain) + Ks'*(K\(fTrain'-min(fTrain)));
    fTest = mu + Ks'*(K\(fTrain'-mu));

end

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
    for i = 1:1:length(x)-1
        val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
    end
end