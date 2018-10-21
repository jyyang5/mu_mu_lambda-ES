% Use GP estimate and cummulative step size adaptation
% function evaluation for lambda offsprings with GP estimate 
% In each iteration only centroid is evaluated use true objective Function

function val = mml_GP_CSA(f,x0,sigma0,TRAINING_SIZE,lambda,NUM_OF_ITERATIONS)
% initialization
% f:                  objective function value
% x0:                 mu initial point size [n, mu]
% sigma0:             initial step size
% TRAINING_SIZE:      size of GP training set
% lambda:             # of offsprings genenerated in each itertaion  
% mu:                 parent size
% NUM_OF_ITERATIONS:  number of maximum iterations


% OPTIMAL:            global optima
% example input:      f = @(x) x' * x
%                     mml(f,randn(n,mu),1,40,10,1,4000)
[n, mu] = size(x0);

% TRAINING_SIZE = 4*lambda;
xTrain = zeros(n,10000);            % training data for GP size 4*mu
fTrain = zeros(1,10000);


centroid_array = zeros(n,10000);
fcentroid_array = zeros(1,10000);
sigma_array = zeros(1,10000);
error_array = zeros(1,10000);                                               % store similar to noise-to-signal ratio
s_array = zeros(1,10000);
y = zeros(n,lambda);                                                        % lambda offspring solution with dim n
z = zeros(n,lambda);                                                        % random directions added for lambda offsprings dim n
fy = zeros(lambda,1);                                                       % objective function value of y                                              
centroid = mean(x0, 2);                                                     % centroid of parent set, size = [n, 1]
f_centroid = f(centroid);                                                   % fx of centroid
fyep = zeros(lambda,1);                                                     % GP estimate for offsprings
 
convergence_rate = 0;

t = 1;       
T = 1;

centroid_array(:,t) = centroid;
fcentroid_array(t) = f_centroid;

% parameters for CSA
c = 1/sqrt(n);
D = sqrt(n);
s = 0;
%%%%%%%%%%%%%%%%%%%%%%
% nico's NOT WORK
% mu_eff = mu;
% c = (mu_eff+2)/(n+mu_eff+5);
% d = 1 + 2*max(0,sqrt((mu_eff-1)/(n+1))-1)+c; 
% s = 0;
%%%%%%%%%%%%%%%%%%%%%%
% Dirk small lambda 
% c = 0.63;
% d = 1;
% c = 1/sqrt(n);
% D = sqrt(n);
% s = 0;

length_scale_factor = 8;
sigma = sigma0;

while((T < NUM_OF_ITERATIONS) && f_centroid > 10^(-8))
    % early stopping 
    if(f_centroid > 500)
        % if diverge -> convergence rate = 0
        val = {9999999,mean(x0, 2),99999,sigma_array, 9999999, fcentroid_array,-1,error_array,s_array}; 
        return 
    end
%     dist = norm(centroid);                                                 % distance to optimal
%     /n*dist;                                             % mutation strength/step size(temp)  
    
    % (mu/mu, lambda)-ES 4 times to obtain GP traning set
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
        theta = sigma*length_scale_factor*sqrt(n);                          % length scale for GP
        fy_true = zeros(lambda,1);                                          % for calculating the noise 
        % offspring_generation (lambda offspring) 
        for i = 1:1:lambda
            % offspring = mean(parent) + stepsize*z
            z(:,i) = randn(n,1);
            y(:,i) = centroid + sigma*z(:,i);
            fyep(i) = gp(xTrain(:,T-TRAINING_SIZE:T-1), fTrain(T-TRAINING_SIZE:T-1), y(:,i), theta); 
            fy_true(i) = f(y(:,i));                                         % for calculating the noise 
        end
        % calculate relative error
        error_array(t) = var(fyep-fy_true)/var(fy_true);
        % for simple calculation 
        fy = fyep;
        
    
    end
  
    % sort fyep (smaller first)
    [index, sorted_order] = sort(fy);
    z = z(:,sorted_order);
    % choose the best mu candidate solutions as parent
    z = mean(z(:,1:mu),2);
    centroid = centroid + sigma*z;
    f_centroid = f(centroid);
        
    % CSA
    s = (1-c)*s+sqrt(mu*c*(2-c))*z;
    sigma = sigma*exp((norm(s)^2-n)/2/D/n);
    %%%%%%%%%%%%%%%%%%%%%%
    % nico's NOT WORK
%     s = (1-c)*s + sqrt(c*(2-c)*mu)*z;
%     sigma = sigma*exp(c/d*(norm(s)-1)/sqrt(n));

    s_array(t) = (norm(s)^2-n)/2/D/n;
        
    centroid_array(:,t) = centroid;
    fcentroid_array(t) = f_centroid;
    xTrain(:, T) = centroid;                
    fTrain(T) = f_centroid;
    T = T + 1;
    sigma_array(t) = sigma;
    
    t = t + 1;
    
    
    
    
end

    t = t - 1;
    T = T - 1; 

    % convergence rate (overall)
    t_start = ceil(TRAINING_SIZE/lambda);
    convergence_rate = -n/2*sum(log(fcentroid_array(t_start+2:t)./fcentroid_array(t_start+1:t-1)))/(t-t_start-1);
    
    val = {t,centroid,f_centroid,sigma_array, T, fcentroid_array,convergence_rate,error_array,s_array};

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