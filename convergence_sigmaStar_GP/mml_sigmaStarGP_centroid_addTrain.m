% Add another ADD_TRAIN_POINTS points to the training set per iteration
% TRAING_SIZE = 40 fixed
% Add TRAINING_FACTOR to tune training size
% Use normalized step size for step size update
% function evaluation for lambda offsprings with GP estimate 
% In each iteration only the centroid is evaluated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = mml_sigmaStarGP_centroid_addTrain(f,x0,sigma_star,lambda,NUM_OF_ITERATIONS,ADD_TRAIN_POINTS)
% initialization
% f:                  objective function value
% x0:                 mu initial point size [n, mu]
% sigma0:             initial step size
% sigma_ep_star:      normalized noise-to-signal ratio 
% lambda:             # of offsprings genenerated in each itertaion  
% mu:                 parent size
% sigma0:             initial muttaion strength
% NUM_OF_ITERATIONS:  number of maximum iterations
% TRAINING_FACTOR:    TRAINING_SIZE = TRAINING_FACTOR*lambda

% Return 
% 1.t:                  # of objective function calls                    
% 2.centroid:           last parent x(centroid)
% 3.f_centroid:         last objective function value
% 4.sigma_array:        simage arrary over # of objective function calls  
% 5.T:                  # of objective function calls
% 6.fcentroid_array:    objective function values for parents(centroids)
% 7.convergence_rate:   rate of convergence
% 8.median GP_error:    median of GP error
% 9.GP_error_array:     f_centroid(t)-fep_centroid(t))/(f_centroid(t)-f_centroid(t-1)
%                       effective SIZE = t-6
% 10.fep_centroid
% 11.success_rate       # of centroid superior to its predecessors

% OPTIMAL:            global optima
% example input:      f = @(x) x' * x
%                     mml(f,randn(n,mu),1,0,10,1,4000)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n, mu] = size(x0);
TRAINING_SIZE = 40;
%TRAINING_SIZE = lambda*TRAINING_FACTOR;
xTrain = zeros(n,10000);                    % training data for GP size 4*mu
fTrain = zeros(1,10000);


centroid_array = zeros(n,10000);
fcentroid_array = zeros(1,10000);
sigma_array = zeros(1,10000);
fep_centroid = zeros(1,10000);
theta_array = zeros(1,10000);


y = zeros(n,lambda);                        % lambda offspring solution with dim n
z = zeros(n,lambda);                        % random directions added for lambda offsprings dim n
fy = zeros(lambda,1);                       % objective function value of y                                              
centroid = mean(x0, 2);                     % centroid of parent set, size = [n, 1]
f_centroid = f(centroid);                   % fx of centroid
fyep = zeros(lambda,1);                     % GP estimate for offsprings
GP_error = zeros(1,10000);                  % relative error of GP 


convergence_rate = 0;

t = 1;       
T = 1;


centroid_array(:,t) = centroid;
fcentroid_array(t) = f_centroid;


length_scale_factor = 8;


while((T < NUM_OF_ITERATIONS) && f_centroid > 10^(-8))
    
    if(f_centroid > 10)
        val = {9999999,mean(x0, 2),99999,sigma_array, 9999999, fcentroid_array,0,median(GP_error(1:t-5)),1,fep_centroid, 0,GP_error}; 
        return 
    end
    dist = norm(centroid);                              % distance to optimal
    sigma = sigma_star/n*dist;                          % mutation strength/step size(temp) 
    theta = sigma*length_scale_factor*sqrt(n);          % length scale for GP
    
    % (mu/mu, lambda)-ES 4 times to obtain GP traning set
    if T < TRAINING_SIZE
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
        % offspring_generation (lambda offspring) 
        for i = 1:1:lambda
            % offspring = mean(parent) + stepsize*z
            z(:,i) = randn(n,1);
            y(:,i) = centroid + sigma*z(:,i);
            fyep(i) = gp(xTrain(:,T-TRAINING_SIZE:T-1), fTrain(T-TRAINING_SIZE:T-1), y(:,i), theta);
        end
        % for simple calculation 
        fy = fyep;
    
    end
    
    % sort fyep (smaller first)
    [index, sorted_order] = sort(fy);
    z = z(:,sorted_order);
    if(T >= TRAINING_SIZE && ADD_TRAIN_POINTS~=0)
        % add ADD_TRAIN_POINTS best points to training set
        for i=1:1:ADD_TRAIN_POINTS
            xTrain(:, T) = centroid + sigma*z(:,i);
            fTrain(T) = f(centroid + sigma*z(:,i));
            T = T + 1;
        end
    end
    % choose the best mu candidate solutions as parent
    z = mean(z(:,1:mu),2);
    centroid = centroid + sigma*z;
    f_centroid = f(centroid);
    
    if(T>TRAINING_SIZE)
        fep_centroid(t) = gp(xTrain(:,T-TRAINING_SIZE:T-1), fTrain(T-TRAINING_SIZE:T-1), y(:,i), theta);
    end  
    xTrain(:, T) = centroid;                
    fTrain(T) = f_centroid;
    T = T + 1;
      
    
    centroid_array(:,t) = centroid;
    fcentroid_array(t) = f_centroid;
    
    theta_array(t) = theta;
    sigma_array(t) = sigma;
    t = t + 1;
%     disp(t);
%      disp(T);
%      disp(f_centroid);
    % update normalized step size array 
    
end
    t = t-1;
    T = T-1;
     % # of iterations to build GP model
    temp = ceil(TRAINING_SIZE/lambda);
    % convergence rate (overall)
    convergence_rate = -n/2*sum(log(fcentroid_array(2:t)./fcentroid_array(1:t-1)))/(t-1);
    % relative error for GP |f(y)-fep(y)|/ |f(y)-f(x)|
    GP_error(1:t-temp) = abs(fep_centroid(temp+1:t)-fcentroid_array(temp+1:t))./abs(fcentroid_array(temp:t-1)-fcentroid_array(temp+1:t));
    GP_temp = GP_error(1:t-temp);
    GP_temp = GP_temp(~isnan(GP_temp));
    GP_temp = GP_temp(~isinf(GP_temp));
    % success rate (# of offspring better than parent)/(# total iteartions) 
    success_rate = sum(fcentroid_array(temp:t-1)>fcentroid_array(temp+1:t))/(t-temp);
    val = {t,centroid,f_centroid,sigma_array, T, fcentroid_array,convergence_rate/(TRAINING_FACTOR+1),median(GP_temp),1,fep_centroid, success_rate,GP_error};

end


function fTest = gp(xTrain, fTrain, xTest, theta)
% input: 
%       xTrain(40 training pts)
%       fTrain(true objective function value)
%       xTest(1 test pt)   
%       theta length scale  
% return: the prediction of input test data

    [n, m] = size(xTrain);                                       % m:  # of training data

    delta = vecnorm(repmat(xTrain, 1, m)-repelem(xTrain, 1, m)); %|x_ij = train_i-train_j| 
    K = reshape(exp(-delta.^2/theta^2/2), m , m);                % K

    deltas = vecnorm(xTrain-repelem(xTest, 1, m));               %|x_ij = train_i-test_j| euclidean distance
    Ks = exp(-(deltas/theta).^2/2)';                             % K_star             

    deltass = vecnorm(repmat(xTest, 1, m)-repelem(xTest, 1, m));
    % Kss = reshape((exp(-deltass.^2/theta^2/2)), m , m);
    
    %Kinv = inv(K);       

    mu = min(fTrain);                                            % estimated mean of GP
    fTest = mu + Ks'*(K\(fTrain'-mu));

end