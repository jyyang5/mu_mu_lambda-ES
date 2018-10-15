% refactored 
% Use Gaussian distributed random noise to model GP estimate replace objective 
% function evaluation for lambda offsprings with GP estimate 
% In each iteration only the centroid is evaluated


function val = mml_noise(f,x0,sigma_star,sigma_ep_star,lambda,NUM_OF_ITERATIONS)
% initialization
% f:                  objective function value
% x0:                 mu initial point size [n, mu]
% sigma_star:         normalized step size
% sigma_ep_star:      normalized noise-to-signal ratio 
% lambda:             # of offsprings genenerated in each itertaion  
% mu:                 parent size
% sigma0:             initial muttaion strength
% NUM_OF_ITERATIONS:  number of maximum iterations



% OPTIMAL:            global optima
% example input:      f = @(x) x' * x
%                     mml(f,randn(n,mu),1,0,10,1,4000)
[n, mu] = size(x0);

TRAINING_SIZE = 40;
xtrain = zeros(n,TRAINING_SIZE);            % training data for GP size 4*mu
fTrain = zeros(1,TRAINING_SIZE);


centroid_array = zeros(n,10000);
fcentroid_array = zeros(1,10000);
sigma_array = zeros(1,10000);
s_array = zeros(1,10000);

y = zeros(n,lambda);                        % lambda offspring solution with dim n
z = zeros(n,lambda);                        % random directions added for lambda offsprings dim n
fy = zeros(lambda,1);                       % objective function value of y                                              
centroid = mean(x0, 2);                     % centroid of parent set, size = [n, 1]
f_centroid = f(centroid);                   % fx of centroid
fyep = zeros(lambda,1);                     % GP estimate for offsprings

convergence_rate = 0;

t = 1;       
T = 1;


centroid_array(:,t) = centroid;
fcentroid_array(t) = f_centroid;

% parameters for CSA
c = 1/sqrt(n);
D = sqrt(n);
s = 0;



while((T < NUM_OF_ITERATIONS) && f_centroid > 10^(-8))
    % early stopping 
    if(f_centroid > 50)
        % if diverge -> convergence rate = 0
        val = {9999999,mean(x0, 2),99999,sigma_array, centroid_array, fcentroid_array,-1,s_array}; 
        return 
    end
    dist = norm(centroid);                  % distance to optimal
    sigma = sigma_star/n*dist;              % mutation strength/step size(temp)  
    
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
        % update theta using 
        dist = norm(centroid);
        sigma_ep = sigma_ep_star/n*2*dist^2;      % Gaussian noise 
        
        % offspring_generation (lambda offspring) 
        for i = 1:1:lambda
            % offspring = mean(parent) + stepsize*z
            z(:,i) = randn(n,1);
            y(:,i) = centroid + sigma*z(:,i);
            fyep(i) = f(y(:,i)) + sigma_ep * randn();
        end
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
    s = (1-c)*s + sqrt(mu*c*(2-c))*z;
        
    centroid_array(:,t) = centroid;
    fcentroid_array(t) = f_centroid;
    xTrain(:, T) = centroid;                
    fTrain(T) = f_centroid;
    T = T + 1;
    sigma_array(t) = sigma;
    s_array(t) = norm(s)^2-n;
    
    t = t + 1;
    
    
    
    
end

    t = t - 1;
    T = T - 1; 

    % convergence rate (overall)
    convergence_rate = -n/2*sum(log(fcentroid_array(2:t)./fcentroid_array(1:t-1)))/(t-1);
    
    %plot(1:1:t-1,s_array(1:t-1));
    
    val = {t,centroid,f_centroid,sigma_array, centroid_array, fcentroid_array,convergence_rate,s_array};
%val = {t,centroid,f_centroid,sigma_array, 1, 1,convergence_rate};

end



function fTest = gp(xTrain, fTrain, xTest, theta)
% input: 
%       xTrain(40 training pts)
%       fTrain(true objective function value)
%       xTest(1 test pt)   
%       theta mutation length     
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