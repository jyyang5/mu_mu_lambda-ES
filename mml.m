function val = mml(f,x0,sigma_star,sigma_ep_star,lambda,sigma0,NUM_OF_ITERATIONS)
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

[n, mu] = size(x0);


x = zeros(n,mu);                           % mu parent solution with dim n 
y = zeros(n,lambda);                       % lambda offspring solution with dim n
fy = zeros(lambda);                        % objective function value of y
sigma = 0;                                 % mutation strength/step size(temp)  
centroid = zeros(n, 1);                    % centroid of parent set 
f_centroid = 100;                          % fx of centroid                     
dist = 0;                                  % distance to optimal
t = 0;       
x = x0;

centroid_array = zeros(n,100);
fcentroid_array = zeros(1,100);
sigma_array = zeros(1,1000);

while((t < NUM_OF_ITERATIONS) && f_centroid > 10^(-8))
    % take centroid of parent 
    dist = norm(x);                           % distance to optimal
    sigma = sigma_star/n*dist;

    centroid = mean(x, 2);
    f_centroid = f(centroid);
    % offspring genneration 
    for i = 1:1:lambda
        % offspring = mean(parent) + stepsize*z
        y(:,i) = centroid + sigma*randn(n,1);
        fy(i) = f(y(:,i)); 
    end

    % sort fyep (smaller first)
    [index, sorted_order] = sort(fy);
    y = y(:,sorted_order);
    % choose the best mu candidate solutions as parent 
    x = y(:,1:mu);
    t = t+1;
    centroid_array(:,t) = centroid;
    fcentroid_array(t) = f_centroid;
    sigma_array(t) = sigma;
end

val = {t,centroid_array,fcentroid_array,sigma_array}; 
end

% TRAINING_SIZE = 4*lambda;

% [mu,n,m] = size(x0);

% x1 = zeros(n,1);                           % parent solution with dim n for (1+1)-ES 
% x = zeros(n,mu);                           % mu parent solution with dim n 
% y1 = zeros(n);                             % candidate solution with dim n for (1+1)-ES
% y = zeros(n,lambda);                       % lambda offspring solution with dim n
% centroid = zeros(n, 1);                    % centroid of parent set 
% f_centroid = 0;                            % fx of centroid                     
% sigma = zeros(1,6000);                     % mutation strength
% xtrain = zeros(n,TRAINING_SIZE);           % training data for GP size 4*mu
% fTrain = zeros(TRAINING_SIZE);

% f_x = zeros(lambda,6000);                  % objectie function value for offspring generated
% fyep = zeros(lambda);                      % GP estimate of the offspring gennerated 

% x_array = zeros(n,6000);                   % store the points evaluated

% D = sqrt(1+n);
% c2 = -0.2;                                 % decrease step size
% c3 = 0.8;                                  % increase step size

% t = 1;                                     % # of iterations
% T = 1;                                     % # of objective function evaluations
% t_gp = 0;


% sigma(1) = sigma0;                         
% x(:, 1) = x0;                               

% while((t < NUM_OF_ITERATIONS) && f(x_array(:,t)) > 10^(-8))

%     % buold the GP model in the first 4*lambda iteraties 
%     if T <= TRAINING_SIZE
%             y1 = x1 +sigma*randn(n,1);
%             fy = f(y1);                       % fitness of offspring(use true objective fn)
%             xTrain(:, T) = y1;                
%             fTrain(T) = fy;
%             dist = norm(y1);
%             if(fy >= fx)                      % bad offspring
%                 sigma = sigma0/n*dist;
%             else
%                 x1 = y;
%                 fx = fy;
%                 sigma = sigma0/n*dist;        % increase step size
%             end
%             % add the point evaluated
%             x_array = x1;
%             % add parent of (mu/mu, lambda)-ES 
%             x(:,rem(T, mu)) = x1;

            
%     % (mu/mu, lambda)-ES 
%     else 

%     % offspring_generation (lambda offspring)
%         centroid = mean(x, 2);
%         f_centroid = f(centroid);
%         xTrain(:, T) = centroid;                
%         fTrain(T) = f_centroid;

%         for i = 1:1:lambda
%             % offspring = mean(parent) + stepsize*z
%             y_temp = centroid + sigma(t)*randn(n,1);
%             fyep_temp = gp(xTrain, fTrain, y_temp,theta);
%             % add offspring candidate solution iff. GP estimate superior to centroid
%             while(f_centroid < fyep_temp)
%                 y_temp = centroid + sigma(t)*randn(n,1);
%                 fyep_temp = gp(xTrain, fTrain, y_temp,theta);
%                 t_gp = t_gp + 1;
%             end
%             % update  candidate solution
%             y(:,i) = y_temp; 
%             fyep(i) = fyep_temp;

%         end

%         % sort fyep (smaller first)
%         [index, sorted_order] = sort(fyep);
%         y = y(:,sorted_order);
%         % choose the best mu candidate solutions as parent 
%         x = y(:, 1:mu);
        
    
%     end

%     t = t + 1;
%     T = T + 1

    
% end 

% x_last = x(:,t);

% val = {t, x_last, f_x, sigma};

% end



% function fTest = gp(xTrain, fTrain, xTest, theta)
% % input: 
% %       xTrain(40 training pts)
% %       fTrain(true objective function value)
% %       xTest(1 test pt)   
% %       theta mutation length     
% % return: the prediction of input test data

%     [n, m] = size(xTrain);                                       % m:  # of training data

%     delta = vecnorm(repmat(xTrain, 1, m)-repelem(xTrain, 1, m)); %|x_ij = train_i-train_j| 
%     K = reshape(exp(-delta.^2/theta^2/2), m , m);                % K

%     deltas = vecnorm(xTrain-repelem(xTest, 1, m));               %|x_ij = train_i-test_j| euclidean distance
%     Ks = exp(-(deltas/theta).^2/2)';                             % K_star             

%     deltass = vecnorm(repmat(xTest, 1, m)-repelem(xTest, 1, m));
%     % Kss = reshape((exp(-deltass.^2/theta^2/2)), m , m);
%     if(det(K) < 10^(-18))
%             %disp(K);
        
%     end
%     %Kinv = inv(K);       

%     mu = min(fTrain);                                            % estimated mean of GP
%     fTest = mu + Ks'*(K\(fTrain'-mu));

% end