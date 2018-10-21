% usual mml-ES using CSA no GP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = mml_CSA(f,x0,sigma0,lambda,NUM_OF_ITERATIONS)
% initialization
% f:                  objective function value
% x0:                 mu initial point size [n, mu]
% sigma0:             initial step size
% sigma_ep_star:      normalized noise-to-signal ratio 
% lambda:             # of offsprings genenerated in each itertaion  
% mu:                 parent size
% sigma0:             initial muttaion strength
% NUM_OF_ITERATIONS:  number of maximum iterations

% Return 
% 1.t:                  # of objective function calls                    
% 2.centroid:           last parent x(centroid)
% 3.f_centroid:         last objective function value
% 4.sigma_array:        simage arrary over # of objective function calls  
% 5.T:                  # of objective function calls
% 6.fcentroid_array:    objective function values for parents(centroids)
% 7.convergence_rate:   rate of convergence
% 8.1:                  no GP estimate
% 9.sigma_star_array:   normalized step size


% OPTIMAL:            global optima
% example input:      f = @(x) x' * x
%                     mml(f,randn(n,mu),1,0,10,1,4000)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n, mu] = size(x0);

TRAINING_SIZE = 4*lambda;
xTrain = zeros(n,10000);                    % training data for GP size 4*mu
fTrain = zeros(1,10000);


centroid_array = zeros(n,10000);
fcentroid_array = zeros(1,10000);
sigma_array = zeros(1,10000);
s_array = zeros(1,10000);
fep_centroid = zeros(1,10000);
theta_array = zeros(1,10000);
sigma_star_array = zeros(1,10000);          % store normalized step size


y = zeros(n,lambda);                        % lambda offspring solution with dim n
z = zeros(n,lambda);                        % random directions added for lambda offsprings dim n
fy = zeros(lambda,1);                       % objective function value of y                                              
centroid = mean(x0, 2);                     % centroid of parent set, size = [n, 1]
f_centroid = f(centroid);                   % fx of centroid
fyep = zeros(lambda,1);                     % GP estimate for offsprings


convergence_rate = 0;
success_rate = 0;

t = 1;       
T = 1;


centroid_array(:,t) = centroid;
fcentroid_array(t) = f_centroid;

% parameters for CSA
c = 1/sqrt(n);
D = sqrt(n);
s = 0;

sigma = sigma0;

while((t < NUM_OF_ITERATIONS) && f_centroid > 10^(-8))
    % early stopping
    if(f_centroid > 500)
        % if diverge -> convergence rate = 0
        val = {999999,mean(x0, 2),99999,sigma_array, 9999999, fcentroid_array,-1,error_array,sigma_star_array,error_array}; 
        return 
    end
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
                    
    % (mu/mu, lambda)-ES use GP estiate 
    else  
       
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
    T = T + lambda;
    
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
    sigma = sigma*exp((norm(s)^2-n)/2/D/n);

%     sigma = sigma*exp((norm(s)^2-n)/(2*D*n));
%     theta = sigma*8*sqrt(n);
    
    centroid_array(:,t) = centroid;
    fcentroid_array(t) = f_centroid;
    xTrain(:, T) = centroid;                
    fTrain(T) = f_centroid;
    
    %fep_centroid(t) = gp(xTrain, fTrain, centroid, theta);
    theta_array(t) = theta;
    
    sigma_array(t) = sigma;
    
    t = t + 1;
    T = T + 1;
    
    dist = norm(centroid);
    sigma_star = sigma*n/dist;
    sigma_star_array(t)=sigma_star;
    
    
    
end
    t = t - 1;
    T = T - 1; 
 
    % convergence rate
    convergence_rate = -n/2*sum(log(fcentroid_array(2:t)./fcentroid_array(1:t-1)))/(t-1);
    % success rate 
    success_rate = sum((fcentroid_array(1:t-1)-fcentroid_array(2:t))>0)/(T-1);
    %plot(1:1:t-1,s_array(1:t-1));
    
    val = {t,centroid,f_centroid,sigma_array, T, fcentroid_array, convergence_rate,-1,sigma_star_array,success_rate};

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