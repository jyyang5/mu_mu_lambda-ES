function val = mml_GP(f,x0,sigma_star,sigma_ep_star,lambda,sigma0,NUM_OF_ITERATIONS)
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

TRAINING_SIZE = 4*lambda;
xtrain = zeros(n,TRAINING_SIZE);            % training data for GP size 4*mu
fTrain = zeros(TRAINING_SIZE);


centroid_array = zeros(n,100);
fcentroid_array = zeros(1,100);
sigma_array = zeros(1,1000);

y = zeros(n,lambda);                        % lambda offspring solution with dim n
fy = zeros(lambda,1);                         % objective function value of y                                              
centroid = mean(x0, 2);                     % centroid of parent set, size = [n, 1]
f_centroid = f(centroid);                   % fx of centroid
fyep = zeros(lambda,1);                       % GP estimate for offsprings

t = 1;       
T = 1;
t_gp = 0;


centroid_array(:,t) = centroid;
fcentroid_array(t) = f_centroid;

sigma(1) = sigma0;                         
x(:, 1) = x0;                               

while((t < NUM_OF_ITERATIONS) && f_centroid > 10^(-8))

    % (mu/mu, lambda)-ES 4 times to obtain GP traning set
    if t <= 4
        dist = norm(centroid);                  % distance to optimal
        sigma = sigma_star/n*dist;              % mutation strength/step size(temp)  
    
        sigma_array(t) = sigma;
    
        % offspring genneration 
        for i = 1:1:lambda
            % offspring = mean(parent) + stepsize*z
            y(:,i) = centroid + sigma*randn(n,1);
            fy(i) = f(y(:,i)); 
        end
        xtrain(:,(t-1)*lambda+1:t*lambda) = y;
        fTrain((t-1)*lambda+1:t*lambda) = fy;  
        
        % sort fyep (smaller first)
        [index, sorted_order] = sort(fy);
        y = y(:,sorted_order);
        % choose the best mu candidate solutions as parent 
        centroid = mean(y(:,1:mu), 2);
        f_centroid = f(centroid);
        t = t+1;
    
        centroid_array(:,t) = centroid;
        fcentroid_array(t) = f_centroid;
        

            
    % (mu/mu, lambda)-ES use GP estiate 
    else 

    % offspring_generation (lambda offspring)
        xTrain(:, rem(t, 40)+1) = centroid;                
        fTrain(rem(t, 40)+1) = f_centroid;

        for i = 1:1:lambda
            % offspring = mean(parent) + stepsize*z
            y_temp = centroid + sigma*randn(n,1);
            fyep_temp = gp(xTrain, fTrain, y_temp,theta);
            % add offspring candidate solution iff. GP estimate superior to centroid
            while(f_centroid < fyep_temp)
                y_temp = centroid + sigma(t)*randn(n,1);
                fyep_temp = gp(xTrain, fTrain, y_temp,theta);
                t_gp = t_gp + 1;
            end
            % update  candidate solution
            y(:,i) = y_temp; 
            fyep(i) = fyep_temp;

        end
        % selection 
        % sort fyep (smaller first)
        [index, sorted_order] = sort(fyep);
        y = y(:,sorted_order);
        % choose the best mu candidate solutions as parent 
        % centroid 
        centroid = mean(y(:, 1:mu), 2);
        f_centroid = f(centroid);
        
        
    
    end

    t = t + 1;
    T = T + 1

    
end 

x_last = x(:,t);

val = {t, x_last, f_x, sigma};

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