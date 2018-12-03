% normal mml-ES
% Evaluate lambda offsprings with true objective function 
% In each iteration only the centroid is added as parent 

function val = mml_normal(f,x0,sigma_star,sigma_ep_star,lambda,NUM_OF_ITERATIONS)
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
fTrain = zeros(TRAINING_SIZE);


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



while((t < NUM_OF_ITERATIONS) && f_centroid > 10^(-8))
    
    dist = norm(centroid);                  % distance to optimal
    sigma = sigma_star/n*dist;              % mutation strength/step size(temp)  
    
    % (mu/mu, lambda)-ES 4 times to obtain GP traning set
    if t <= 4
        % offspring genneration 
        for i = 1:1:lambda
            % offspring = mean(parent) + stepsize*z
            z(:,i) = randn(n,1);
            y(:,i) = centroid + sigma*z(:,i);
            fy(i) = f(y(:,i)); 
        end
        xtrain(:,(t-1)*lambda+1:t*lambda) = y;
        fTrain((t-1)*lambda+1:t*lambda) = fy; 
            
    % (mu/mu, lambda)-ES use GP estiate 
    else  
        xTrain(:, rem(t, 40)+1) = centroid;                
        fTrain(rem(t, 40)+1) = f_centroid;
       
        dist = norm(centroid);
        sigma_ep = sigma_ep_star/n*2*dist^2;      % Gaussian noise 
        
        % offspring_generation (lambda offspring) 
        for i = 1:1:lambda
            % offspring = mean(parent) + stepsize*z
            z(:,i) = randn(n,1);
            y(:,i) = centroid + sigma*z(:,i);
            fyep(i) = f(y(:,i));
        end
        % for simple calculation 
        fy = fyep;
        
    
    end
    
    
    % sort fyep (smaller first)
    [index, sorted_order] = sort(fy);
    z = z(:,sorted_order);
    % choose the best mu candidate solutions as parent
    z = mean(z(:,1:mu),2);
%     z = (mean(y(:,1:mu),2)-centroid)/sigma;           % AVG(step of offspring chosen)
%     centroid = mean(y(:,1:mu), 2);
    centroid = centroid + sigma*z;
    f_centroid = f(centroid);
        
    % CSA
    s = (1-c)*s + sqrt(mu*c*(2-c))*z;
    
%     sigma = sigma*exp((norm(s)^2-n)/(2*D*n));
    
    centroid_array(:,t) = centroid;
    fcentroid_array(t) = f_centroid;
    sigma_array(t) = sigma;
    
    t = t + 1;
    T = T + 1;
    
    
    
end

    for i = 1:1:t-2
        convergence_rate = convergence_rate + (log(fcentroid_array(i+1)/fcentroid_array(i)));
    end
    % convergence rate divided by # of objective function
    % evaluation/iteration
    convergence_rate = -n/2*convergence_rate/(t-2)/(lambda+1);
    
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

%     deltass = vecnorm(repmat(xTest, 1, m)-repelem(xTest, 1, m));
    % Kss = reshape((exp(-deltass.^2/theta^2/2)), m , m);
    %Kinv = inv(K);       

    mu = min(fTrain);                                            % estimated mean of GP
    fTest = mu + Ks'*(K\(fTrain'-mu));

end