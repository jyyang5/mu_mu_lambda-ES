% CSA uses z = randn(n,1) in offspring genenration 
% two approaches in adapting c
% a) take the mean of mean(AVG(z(1:mu))), c is 1-dim 
% b) use the averaged n-dimensional for c, c is n-dim but later take l2-nurm still work 
% for performance, b) works slightly better
function val = mml_noise_csa(f,x0,sigma0,sigma_GP,sigma_ep_star,lambda,NUM_OF_ITERATIONS)
% initialization
% f:                  objective function value
% x0:                 mu initial point size [n, mu]
% sigma_GP:           a scalar multiplied when compare the fep with fcentroid
% sigma_ep_star:      normalized noise-to-signal ratio 
% lambda:             # of offsprings genenerated in each itertaion  
% mu:                 parent size
% sigma0:             initial muttaion strength
% NUM_OF_ITERATIONS:  number of maximum iterations



% OPTIMAL:            global optima
% example input:      f = @(x) x' * x
%                     mml(f,randn(n,mu),1,0,10,1,4000)
[n, mu] = size(x0);

TRAINING_SIZE = 4*lambda;
xtrain = zeros(n,TRAINING_SIZE);            % training data for GP size 4*mu
fTrain = zeros(TRAINING_SIZE);


centroid_array = zeros(n,10000);
fcentroid_array = zeros(1,10000);
sigma_array = zeros(1,10000);

y = zeros(n,lambda);                        % lambda offspring solution with dim n
fy = zeros(lambda,1);                       % objective function value of y                                              
centroid = mean(x0, 2);                     % centroid of parent set, size = [n, 1]
f_centroid = f(centroid);                   % fx of centroid
fyep = zeros(lambda,1);                     % GP estimate for offsprings

convergence_rate = 0;

t = 1;       
T = 1;
t_gp = 0;


centroid_array(:,t) = centroid;
fcentroid_array(t) = f_centroid;

% parameters for CSA
c = 1/sqrt(n);
D = sqrt(n);
s = 0;

sigma = sigma0;
sigma_array(t) = sigma;


while((t < NUM_OF_ITERATIONS) && f_centroid > 10^(-8))

    % (mu/mu, lambda)-ES 4 times to obtain GP traning set
    if t <= 4
%         dist = norm(centroid);                  % distance to optimal
%         sigma = sigma_star/n*dist;              % mutation strength/step size(temp)  
%         
    
        % offspring genneration 
        for i = 1:1:lambda
            % offspring = mean(parent) + stepsize*z
            y(:,i) = centroid + sigma*randn(n,1);
            fy(i) = f(y(:,i)); 
        end
        xtrain(:,(t-1)*lambda+1:t*lambda) = y;
        fTrain((t-1)*lambda+1:t*lambda) = fy; 
            
    % (mu/mu, lambda)-ES use GP estiate 
    else 

    % offspring_generation (lambda offspring)  
        xTrain(:, rem(t, 40)+1) = centroid;                
        fTrain(rem(t, 40)+1) = f_centroid;
        
%         dist = norm(centroid);                  % distance to optimal
%         sigma = sigma_star/n*dist;              % mutation strength/step size(temp)  
        dist = norm(centroid);
        sigma_ep = sigma_ep_star/n*2*dist^2;      % Gaussian noise 
        
        for i = 1:1:lambda
            % offspring = mean(parent) + stepsize*z
            y_temp = centroid + sigma*randn(n,1);
            fyep_temp = f(y_temp) + sigma_ep * randn();
            % add offspring candidate solution iff. GP estimate superior to centroid
            while(f_centroid*sigma_GP < fyep_temp)
                y_temp = centroid + sigma*randn(n,1);
                fyep_temp = f(y_temp) + sigma_ep * randn();
                %fyep_temp = gp(xTrain, fTrain, y_temp,theta);
                t_gp = t_gp + 1;
            end
            % update  candidate solution
            y(:,i) = y_temp; 
            fyep(i) = fyep_temp;

        end
        % for simple calculation 
        fy = fyep;
        
%         % selection 
%         % sort fyep (smaller first)
%         [index, sorted_order] = sort(fyep);
%         y = y(:,sorted_order);
%         % choose the best mu candidate solutions as parent 
%         % centroid 
%         z = mean(y(:,1:mu),2)-centroid;            % AVG(step of offspring chosen)         
%         centroid = mean(y(:, 1:mu), 2);
%         f_centroid = f(centroid);
%         
%         % CSA
%         s = (1-c)*s + sqrt(mu*c(2-c))*z;
%         sigma = sigma*exp((norm(s)-n)/(2*D*n));
        
        
        
    
    end
    
    disp(t_gp)
    % sort fyep (smaller first)
    [index, sorted_order] = sort(fy);
    y = y(:,sorted_order);
    % choose the best mu candidate solutions as parent
    z = mean(y(:,1:mu),2)-centroid;           % AVG(step of offspring chosen)
    centroid = mean(y(:,1:mu), 2);
    f_centroid = f(centroid);
        
    % CSA
    s = (1-c)*s + sqrt(mu*c*(2-c))*mean(z);
%     s = (1-c)*s + sqrt(mu*c*(2-c))*z; 
    sigma = sigma*exp((norm(s)-n)/(2*D*n));
    
    centroid_array(:,t) = centroid;
    fcentroid_array(t) = f_centroid;
    sigma_array(t+1) = sigma;
    
    t = t + 1;
    T = T + 1;
    
    
    
end

    for i = 1:1:t-2
        convergence_rate = convergence_rate + (log(fcentroid_array(i+1)/fcentroid_array(i)));
    end

    convergence_rate = -n/2*convergence_rate/(t-2);

val = {t,centroid,f_centroid,sigma_array, centroid_array, fcentroid_array,convergence_rate,t_gp};
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
    if(det(K) < 10^(-18))
            %disp(K);
        
    end
    %Kinv = inv(K);       

    mu = min(fTrain);                                            % estimated mean of GP
    fTest = mu + Ks'*(K\(fTrain'-mu));

end