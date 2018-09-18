% Use Gaussian distributed random noise to model GP estimate replace objective 
% function evaluation for lambda offsprings with GP estimate 
% In each iteration only the centroid is evaluated


function val = mml_GP(f,x0,sigma0,lambda,NUM_OF_ITERATIONS)
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
% 5.centroid_array:     parent set(centroids)
% 6.fcentroid_array:    objective function values for parents(centroids)
% 7.convergence_rate:   rate of convergence
% 8.GP_error:           mean[f_centroid(t)-fep_centroid(t))/(f_centroid(t)-f_centroid(t-1))] 
% 9.sigma_star_array:   normalized step size


% OPTIMAL:            global optima
% example input:      f = @(x) x' * x
%                     mml(f,randn(n,mu),1,0,10,1,4000)
[n, mu] = size(x0);

TRAINING_SIZE = 4*lambda;
xTrain = zeros(n,TRAINING_SIZE);            % training data for GP size 4*mu
fTrain = zeros(1,TRAINING_SIZE);


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

t = 1;       
T = 1;


centroid_array(:,t) = centroid;
fcentroid_array(t) = f_centroid;

% parameters for CSA
c = 1/sqrt(n);
D = sqrt(n);
s = 0;

sigma = sigma0;
theta = sigma*8*sqrt(n);

GP_error=0;

while((t < NUM_OF_ITERATIONS) && f_centroid > 10^(-8))
    
    % (mu/mu, lambda)-ES 4 times to obtain GP traning set
    if t <= 4
        % offspring genneration 
        for i = 1:1:lambda
            % offspring = mean(parent) + stepsize*z
            z(:,i) = randn(n,1);
            y(:,i) = centroid + sigma*z(:,i);
            fy(i) = f(y(:,i)); 
        end
        xTrain(:,(t-1)*lambda+1:t*lambda) = y;
        fTrain((t-1)*lambda+1:t*lambda) = fy; 
            
    % (mu/mu, lambda)-ES use GP estiate 
    else  
        xTrain(:, rem(t+35, 40)+1) = centroid;                
        fTrain(rem(t+35, 40)+1) = f_centroid;

        % offspring_generation (lambda offspring) 
        for i = 1:1:lambda
            % offspring = mean(parent) + stepsize*z
            z(:,i) = randn(n,1);
            y(:,i) = centroid + sigma*z(:,i);
%             fyep(i) = f(y(:,i)) + sigma_ep * randn();
            fyep(i) = gp(xTrain, fTrain, y(:,i), theta);
            
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
    if(t>4)
        fep_centroid(t) = gp(xTrain, fTrain, centroid, theta);
    end    
    % CSA
    s = (1-c)*s + sqrt(mu*c*(2-c))*z;
    
    sigma = sigma*exp((norm(s)^2-n)/(2*D*n));
    theta = sigma*8*sqrt(n);
    
    centroid_array(:,t) = centroid;
    fcentroid_array(t) = f_centroid;
    
    
    theta_array(t) = theta;
    
    sigma_array(t) = sigma;
    s_array(t) = norm(s)^2-n;
    
    dist = norm(centroid);
    sigma_star = sigma*n/dist;
    sigma_star_array(t)=sigma_star;

    
    t = t + 1;
    T = T + 1;
    
    % update normalized step size array 
    
   
    
    
end
    
    % convergence rate
    convergence_rate = -n/2*sum(log(fcentroid_array(2:t)./fcentroid_array(1:t-1)))/(t-1);
    % relative error for GP |f(y)-fep(y)|/ |f(y)-f(x)|
    GP_error = mean(abs(fep_centroid(6:t)-fcentroid_array(6:t))./abs(fcentroid_array(5:t-1)-fcentroid_array(6:t)));
    
    val = {t,centroid,f_centroid,sigma_array, centroid_array, fcentroid_array,convergence_rate,GP_error,sigma_star_array};
%val = {t,centroid,f_centroid,sigma_array, 1, 1,convergence_rate};

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