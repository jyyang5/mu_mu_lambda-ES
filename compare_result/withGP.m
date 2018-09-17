
function val = withGP(f,x0,sigma0,NUM_OF_ITERATIONS)%, OPTIMAL, TARGET_DISTANCE)
% initialization
% f:                  objective function value
% x0:                 initial point
% sigma0:             initial muttaion strength
% NUM_OF_ITERATIONS:  number of maximum iterations

% initialization

[n,m] = size(x0);                         % dim of the data
xTrain = zeros(n, 10000);                  % parent solution with dim n               
sigma = sigma0;                           % mutation strength(temp)

% for graphing
sigma_funEva_array = zeros(1,10000);       % store all sigma over function evaluations
x_array = zeros(n,10000);                  % store all best candiate solution
f_x = zeros(1,10000);                      % store all objective function val of candaite solution
sigma_iterate_array = zeros(1,10000);      % store sigma over iterations   
x_array(:,1) = x0;
f_x(1) = f(x0);
sigma_array(1) = sigma0;
sigma_iterate_array(1) = sigma0;


c1 = 0.05;                                % 1/5 rule
c2 = 0.2;
c3 = 0.6;
D = sqrt(n+1);
theta = sigma*8*sqrt(n);



%f = @(x) x'*x;                           % true objective function
%k = @(xy,theta) exp(-norm(xy)^2/theta/2);% square exponential (SE) with input |x-y|

x = x0;%randn(n, 1);                      % best candiate solution dim:n
fx = f(x);                                % function value
xTrain(:, 1) = x;                         % parent solutions
fTrain(1) = fx;                           % vector: value of parent solutions

t = 1;                                    % # of iterations
T = 1;                                    % # of distinct parent solution 



while(t < NUM_OF_ITERATIONS && f(x_array(:,T))>10^(-8))%(norm(x_array(:,T)-OPTIMAL(:,1)) > TARGET_DISTANCE))

    % offspring_generation
    y = x + sigma*randn(n,1);%temp;

    % update GP iff. there is 40 candiate solutions
    if T > 40
        
        theta = sigma*8*sqrt(n);                                         % NOtE: need to be updated every time
        fy_ep = gp(xTrain(:, T-40:T-1), fTrain(T-40:T-1), y, theta);     % fitness of offspring(use GP)
    end
    
    % update mutation & assign new offspring
    % if GP already built compare
    if(T > 40 && fy_ep >= fx)             % bad offspring
        sigma = sigma * exp(-c1/D);
        
    % GP not built || offspring inferior
    else
        fy = f(y);                        % fitness of offspring(use true objective fn)
        xTrain(:, T) = y;                
        fTrain(T) = fy;
        T = T + 1;
        if(fy >= fx)                      % bad offspring                      
            sigma = sigma * exp(-c2/D);   % reduce step size
        else
            x = y;
            fx = fy;
            sigma = sigma * exp(c3/D);   % increase step size
        end
        x_array(:,T) = x;
        f_x(T) = fx;
        sigma_funEva_array(T) = sigma;
        
    end 
    % new iteration     
    t = t + 1;
    
    sigma_iterate_array(T) = sigma;
    % store vars 
    %x_array(:,t) = x;
    %f_x(t) = fx;
    %sigma_array(t) = sigma;
    
    
end 

val = {T, x_array(:,T), f_x,  sigma_funEva_array};



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
%Kss = reshape((exp(-deltass.^2/theta^2/2)), m , m);

%Kinv = inv(K);       

mu = min(fTrain);                                            % estimated mean of GP
%fTest = mu + Ks'*Kinv*(fTrain - mu)';
fTest = mu + Ks'*(K\(fTrain'-mu));
end
