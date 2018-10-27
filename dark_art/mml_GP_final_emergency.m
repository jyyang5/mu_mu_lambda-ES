%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:      Niko's CSA(paper CMA tutorial) and stop in emergency
% Offspring:        Evaluated by GP (ranking)
% Step size:        Adapted by Niko's CSA(paper CMA stop in emergency
% Centroid:         Best so far 
% Version:          Final (fix negative fitness gain/iteration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = mml_GP_final_emergency(fname,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,DECREASE_FACTOR)
% initialization
% f:                  objective function %
%                     1 = linear
%                     2 = quadratic
%                     3 = cubic
% x0:                 mu initial point size [n, mu]
% sigma0:             initial step size
% lambda:             # of offsprings genenerated in each itertaion  
% mu:                 parent size
% NUM_OF_ITERATIONS:  number of maximum iterations
% TRAINING_SIZE:      surrogate training size
% LENGTH_SCALE:       theta in GP

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
% 10. success_rate      success_rate
% 11. normalized delta  (f_centroid(t)-f_centroid(t-1))/factor
%                       linear: factor = R
%                       quadratic: factor = 2*R^2
%                       cubic: factor = 3*R^3
%                       where R=dist(centroid)
% 12. emergency_count   # of emergencies

% OPTIMAL:            global optima
% example input:      f = @(x) x' * x
%                     mml(f,randn(n,mu),1,40,10,1,4000)

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

[n, mu] = size(x0);
% TRAINING_SIZE = 40;
% TRAINING_SIZE = 4*lambda;
xTrain = zeros(n,10000);            % training data for GP size 4*mu
fTrain = zeros(1,10000);


centroid_array = zeros(n,10000);
fcentroid_array = zeros(1,10000);
sigma_array = zeros(1,10000);
sigma_star_array = zeros(1,10000);                                          % store normalized step size 
error_array = zeros(1,10000);                                               % store similar to noise-to-signal ratio
delta_array = zeros(1,10000);                                               % normalized pdf for fitness gain per iteration
y = zeros(n,lambda);                                                        % lambda offspring solution with dim n
z = zeros(n,lambda);                                                        % random directions added for lambda offsprings dim n
fy = zeros(lambda,1);                                                       % objective function value of y     
centroid = mean(x0, 2);                                                     % centroid of parent set, size = [n, 1]
f_centroid = f(centroid);                                                   % fx of centroid
fyep = zeros(lambda,1);                                                     % GP estimate for offsprings
emergency_count = 0;                                                             % emergency counter                                    
convergence_rate = 0;                                                   
sigma = sigma0;

t = 1;       
T = 1;

centroid_array(:,t) = centroid;
fcentroid_array(t) = f_centroid;


% CSA using nico's CMA paper
mu_eff = mu;
c = (mu_eff+2)/(n+mu_eff+5);
D = 1 + 2*max(0,sqrt((mu_eff-1)/(n+1))-1)+c; 
EN = n^0.5*(1-1/(4*n)+1/(21*n^2));
s = 0;


while((T < NUM_OF_ITERATIONS) && f_centroid > 10^(-8))
    % early stopping if diverge
    if(f_centroid > 500)
        % if diverge -> convergence rate = 0 success rate = 0
        success_rate = 0;
        val = {9999,mean(x0, 2),9999,sigma_array, 9999, fcentroid_array,-1,error_array,sigma_star_array,success_rate,delta_array}; 
        return 
    end

    % (mu/mu, lambda)-ES training parese 
    if T <= TRAINING_SIZE
        % offspring genneration 
        for i = 1:1:lambda
            % offspring = mean(parent) + stepsize*z
            z(:,i) = randn(n,1);
            y(:,i) = centroid + sigma*z(:,i);
            fy(i) = f(y(:,i)); 
        end
        xTrain(:,T:T+lambda-1) = y;
        fTrain(T:T+lambda-1) = fy;
        T = T + lambda;
            
    % (mu/mu, lambda)-ES use GP estiate 
    else  
        % update theta  
        theta = sigma*LENGTH_SCALE*sqrt(n);                                 % GP length scale  
        fy_true = zeros(lambda,1);                                          % for model error calculation
        % offspring_generation (lambda offspring) 
        for i = 1:1:lambda
            % offspring = mean(parent) + stepsize*z
            z(:,i) = randn(n,1);
            y(:,i) = centroid + sigma*z(:,i);
            fyep(i) = gp(xTrain(:,T-TRAINING_SIZE:T-1), fTrain(T-TRAINING_SIZE:T-1), y(:,i), theta); 
            fy_true(i) = f(y(:,i));                                         % for calculating the noise 
        end
        % calculate relative model error
        error_array(t) = var(fyep-fy_true)/var(fy_true);
        % for simple calculation 
        fy = fyep;
        
    
    end
  
    % sort fyep (smaller first)
    [index, sorted_order] = sort(fy);
    z = z(:,sorted_order);
    % choose mutation of best mu candidate solutions, take mean for parent
    z = mean(z(:,1:mu),2);
    centroid = centroid + sigma*z;
    f_centroid = f(centroid);
    T = T + 1;
    % update train set
    xTrain(:, T) = centroid;                
    fTrain(T) = f_centroid;
    % current centroid inferior tiger emergency if GP model is built

    if(T > TRAINING_SIZE && t >= 2 && f_centroid > fcentroid_array(t-1))   
        sigma = sigma*DECREASE_FACTOR;          % decrease step size
        centroid = centroid_array(:,t-1);
        f_centroid = fcentroid_array(t-1);
        emergency_count = emergency_count + 1;
    else % GP model not yet built || superior centroid 
        % Use CSA 
        s = (1-c)*s+sqrt(mu*c*(2-c))*z;
        sigma = sigma*exp(c/D*(norm(s)/EN-1));
    end
    % store updates 
    centroid_array(:,t) = centroid;
    fcentroid_array(t) = f_centroid;
    sigma_array(t) = sigma;
    sigma_star_array(t) = sigma*n/norm(centroid);
    
    t = t + 1;
    
    
end
    
    t = t - 1;
    T = T - 1; 
    
    % plot # of emergencey 
    fprintf('Occurance of emergency: %d/%d\n',emergency_count,t);
    
    % convergence rate (overall)
    t_start = ceil(TRAINING_SIZE/lambda);                                   % first iteration when GP is built
    % convergnece rate and normalized fitness gain/iteration
    if(fname==1)
        delta_array(2:t) = -(fcentroid_array(2:t)-fcentroid_array(1:t-1))./vecnorm(centroid_array(:,2:t),2,1);
        convergence_rate = -n*sum(log(fcentroid_array(t_start+2:t)./fcentroid_array(t_start+1:t-1)))/length(fcentroid_array(t_start+1:t-1));
    elseif(fname==2)
        delta_array(2:t) = -(fcentroid_array(2:t)-fcentroid_array(1:t-1))./(vecnorm(centroid_array(:,2:t),2,1)).^2/2;
        convergence_rate = -n/2*sum(log(fcentroid_array(t_start+2:t)./fcentroid_array(t_start+1:t-1)))/length(fcentroid_array(t_start+1:t-1));
    elseif(fname==3)
        delta_array(2:t) = -(fcentroid_array(2:t)-fcentroid_array(1:t-1))./(vecnorm(centroid_array(:,2:t),2,1)).^3/3;        
        convergence_rate = -n/3*sum(log(fcentroid_array(t_start+2:t)./fcentroid_array(t_start+1:t-1)))/length(fcentroid_array(t_start+1:t-1));
    elseif(fname==4)
        delta_array(2:t) = -(fcentroid_array(2:t)-fcentroid_array(1:t-1))./vecnorm(centroid_array(:,2:t),2,1);
        convergence_rate = -n/2*sum(log(fcentroid_array(t_start+2:t)./fcentroid_array(t_start+1:t-1)))/length(fcentroid_array(t_start+1:t-1));
    elseif(fname==5)
    	delta_array(2:t) = -(fcentroid_array(2:t)-fcentroid_array(1:t-1))./vecnorm(centroid_array(:,2:t),2,1);
        convergence_rate = -n/2*sum(log(fcentroid_array(t_start+2:t)./fcentroid_array(t_start+1:t-1)))/length(fcentroid_array(t_start+1:t-1));
    end
    
    % success rate
    success_rate = sum(fcentroid_array(t_start:T-1)>fcentroid_array(t_start+1:T))/length(fcentroid_array(t_start:T-1));
    % calculate the rate of emergency
    emergency_rate = emergency_count/(t-t_start);
    val = {t,centroid,f_centroid,sigma_array, T, fcentroid_array,convergence_rate,error_array,sigma_star_array,success_rate,delta_array,emergency_rate};

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
    mu = min(fTrain);                                            % estimated mean of GP
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
    n = length(x);
    for i = 1:1:n-1
        val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
    end
end