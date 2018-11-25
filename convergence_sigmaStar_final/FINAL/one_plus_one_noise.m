
function val = one_plus_one_noise(f,x0,sigma_star,sigma_ep_star,NUM_OF_ITERATIONS)%, OPTIMAL, TARGET_DISTANCE)
% initialization
% fname:                objective function value
% x0:                   initial point
% sigma0:               initial muttaion strength
% NUM_OF_ITERATIONS:    number of maximum iterations
% TRAINING_SIZE:        surrogate training size
% LENGTH_SCALE:         theta in GP

% Return 
% 1.T:                  # of objective function calls                    
% 2.x_array(:,T):       last x
% 3.fx:                 last objective function value
% 4.sigma_funEva_array: simage arrary over # of objective function calls  
% 5.T:                  # of objective function calls
% 6.f_x:                objective function values for parents
% 7.convergence_rate:   rate of convergence
% 8.GP error:           relative GP error (after GP is built) SIZE = T-41
% 9.sigma_star_array:   normalized step size
% 10. success_rate      success_rate
% 11. delta_array       normalized delta  (f_centroid(t)-f_centroid(t-1))/factor
%                       linear: factor = R
%                       quadratic: factor = 2*R^2
%                       cubic: factor = 3*R^3
%                       where R=dist(centroid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test functions
% f1 = @(x) (x'*x)^(1/2);  % linear sphere
% f2 = @(x) (x'*x);        % quadratic sphere
% f3 = @(x) (x'*x)^(3/2);  % cubic sphere
% if(fname==1)
%     f=f1;
% elseif(fname==2)
%     f=f2;
% elseif(fname==3)
%     f=f3;
% elseif(fname==4)
%     f=@f4;
% elseif(fname==5)
%     f=@f5;
% end

NUM_OF_ITERATIONS = 20000;
% initialization
[n,m] = size(x0);                         % dim of the data
xTrain = zeros(n, 10000);                 % parent solution with dim n               
% TRAINING_SIZE = 40;

% for graphing
sigma_funEva_array = zeros(1,10000);      % store all sigma over function evaluations
x_array = zeros(n,10000);                 % store all best candiate solution
f_x = zeros(1,10000);                     % store all objective function val of candaite solution
sigma_iterate_array = zeros(1,10000);     % store sigma over iterations   
x_array(:,1) = x0;
f_x(1) = f(x0);
sigma_star_array = zeros(1,10000);        % store normalized step size


c1 = 0.05;                                % 1/5 rule
c2 = 0.2;
c3 = 0.6;
D = sqrt(n+1);


%f = @(x) x'*x;                           % true objective function
%k = @(xy,theta) exp(-norm(xy)^2/theta/2);% square exponential (SE) with input |x-y|

x = x0;%randn(n, 1);                      % best candiate solution dim:n
fx = f(x);                                % function value
xTrain(:, 1) = x;                         % parent solutions
fTrain(1) = fx;                           % vector: value of parent solutions

t = 1;                                    % # of iterations
T = 1;                                    % # of distinct parent solution 



while(t < NUM_OF_ITERATIONS && fx>10^(-8))%(norm(x_array(:,T)-OPTIMAL(:,1)) > TARGET_DISTANCE))
    if(t>100000 || fx>200)
        convergence_rate = 0;
        val = {t, x_array(:,T), fx, sigma_array, T, f_x, convergence_rate};
        return 
    end
    % initial sigma_satr 
    dist = norm(x);                         % distance to optimal
    sigma = sigma_star/n*dist;              % mutation strength/step size(temp)  
    % offspring_generation
    y = x + sigma*randn(n,1);
    % update GP iff. there is 40 candiate solutions
    sigma_ep = sigma_ep_star/n*2*dist^2;      % Gaussian noise 
    fy_ep = f(y)+ sigma_ep * randn();
    if(fy_ep<fx)
        fy = f(y)
        T = T+1;
        if(fy<fx)
            x = y;
            fx = f(x);
        end
        
    end
    x_array(:,T) = x;
    f_x(T) = fx;
    sigma_array(T) = sigma;
    
    % new iteration     
    t = t + 1;
end 
    % convergence rate (overall)
    t_start = 1;
    convergence_rate = -n/2*sum(log(f_x(t_start+1:T)./f_x(t_start:T-1)))/length(f_x(t_start:T-1));
%     % success rate
%     success_rate = sum(f_x(t_start:T-1)>f_x(t_start+1:T))/length(f_x(t_start:T-1));

    val = {t, x_array(:,T), fx, sigma_array, T, f_x, convergence_rate};



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

    % deltass = vecnorm(repmat(xTest, 1, m)-repelem(xTest, 1, m));
    %Kss = reshape((exp(-deltass.^2/theta^2/2)), m , m);

    %Kinv = inv(K);       

    mu = min(fTrain);                                            % estimated mean of GP
    %fTest = mu + Ks'*Kinv*(fTrain - mu)';
    fTest = mu + Ks'*(K\(fTrain'-mu));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schwefel's Problem 1.2
function val = f4(x)
    val = 0;
    for i = 1:1:length(x)
        val = val + sum(x(1:i))^2;
    end
end
% quartic function
function val = f5(x)
    beta = 1;
    val = 0;
    n = length(x);
    for i = 1:1:n-1
        val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
    end
end
