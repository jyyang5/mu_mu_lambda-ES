% [Refactored] training (1+1)-ES + after [bestSoFar] 3 clusters of testFUNs
% step size adaptation: uses 3 cases C1, C2, C3
% function evaluation for lambda offsprings with GP estimate 
% In each iteration only centroid is evaluated use true objective Function

function val = bestSoFar_arashVariant(fname,para,x0,sigma0,lambda,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE,C1,C2,C3)
% initialization
% fname:              an index 
%                       1 for linear
%                       2 for quadratic 
%                       3 for cubic 
%                       4 for schwefel
%                       5 for quartic
%                       6 for generalized sphere functions
%                       7 for quartic function with varying 
%                       8 for ellipsoids function
% para:
% x0:                 mu initial point size [n, mu]
% sigma0:             initial step size
% lambda:             # of offsprings genenerated in each itertaion  
% mu:                 parent size
% NUM_OF_ITERATIONS:  number of maximum iterations
% TRAINING_SIZE:      surrogate training size
% LENGTH_SCALE:       theta in GP
% S:                  success rate
% predicted_bad_ratio:decrease factor when predicted bad 


% Return 
% 1.t:                  # of objective function calls                    
% 2.x_array(:,t):       last x
% 3.fx:                 last objective function value
% 4.sigma:              simage arrary over # of objective function calls  
% 5.T:                  # of objective function calls
% 6.f_x:                objective function values for parents
% 7.convergence_rate:   rate of convergence
% 8.error_array :                 no GP error
% 9.sigma_star_array:   normalized step size
% 10.success_rate       success_rate
% 11.delta_array        normalized delta  (f_centroid(t)-f_centroid(t-1))/factor
%                       linear: factor = R
%                       quadratic: factor = 2*R^2
%                       cubic: factor = 3*R^3
%                       where R=dist(centroid)
% 12. FOUR_COUNT        [TN,FP,FN,TP]
% 13. evaluation ratio  t/iteration 

% OPTIMAL:            global optima
% example input:      f = @(x) x' * x
%                     mml(f,randn(n,mu),1,40,10,1,4000)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test functions
f1 = @(x) (x'*x)^(1/2);  % linear sphere
f2 = @(x) (x'*x);        % quadratic sphere
f3 = @(x) (x'*x)^(3/2);  % cubic sphere
if(fname==1)
    f=f1;
elseif(fname==2)
    f=f2;
elseif(fname==3)
    f=f3;
elseif(fname==4)
    f=@f4;
elseif(fname==5)
    f=@f5;
elseif(fname==6)
    f=@f6;
elseif(fname==7)
    f=@f7;
elseif(fname==8)
    f=@f8;
elseif(fname==9)
    f=@f9;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n, ~] = size(x0);
mu = ceil(lambda/4);
xTrain = zeros(n,50000);            % training data for GP size 4*mu
fTrain = zeros(1,50000);


centroid_array = zeros(n,50000);
fcentroid_array = zeros(1,50000);


sigma_array = zeros(1,50000);
sigma_star_array = zeros(1,50000);                   % store normalized step size 
error_array = zeros(1,50000);                             % store similar to noise-to-signal ratio
delta_array = zeros(1,50000);                             % for pdf of fitness gain in each iteration

y = zeros(n,lambda);                                           % lambda offspring solution with dim n
z = zeros(n,lambda);                                           % random directions added for lambda offsprings dim n
fy = zeros(lambda,1);                                          % objective function value of y     
centroid = x0;                                                     % centroid of parent set, size = [n, 1]
f_centroid = f(centroid);                                      % fx of centroid
fyep = zeros(lambda,1);                                      % GP estimate for offsprings

fy_true = zeros(1,lambda);                                  % for calculating GP error                            
convergence_rate = 0;

t = 1;       
T = 1;
iteration = 1;

centroid_array(:,t) = centroid;
fcentroid_array(t) = f_centroid;
xTrain(:, T) = centroid;                
fTrain(T) = f_centroid;


sigma = sigma0;
D = sqrt(1+n); 
% R = norm(centroid);
% sigma = SIGMA_STAR*R/n;
sigma_array(t) = sigma;
sigma_star_array(t) = n*sigma/norm(centroid);
FOUR_COUNT = zeros(1,4);

% Compute 4 probs
fep_y_array = zeros(1,50000);   % predicted function value of offspring 
f_y_array = zeros(1,50000);     % true function value of offspring 
f_x_array = zeros(1,50000);     % true function value of parent 

while((T < NUM_OF_ITERATIONS) && f_centroid > 10^(-8))
    % early stopping 
    % early stopping 
    if(fname == 6)
        if(f_centroid> 10^30 || sigma <  10^-90)
            % if diverge -> convergence rate = 0 success rate = 0
            success_rate = 0;
            FOUR_COUNT = [0,0,0,0];
            val = {t,centroid,f_centroid,sigma_array, 9999, fcentroid_array,convergence_rate,error_array,sigma_star_array,success_rate,delta_array,FOUR_COUNT,0};
            return 
                        
        end
    elseif(f_centroid > 50000 || iteration >= 500000 || sigma <  10^-45)
        % if diverge -> convergence rate = 0 success rate = 0
        success_rate = 0;
        FOUR_COUNT = [0,0,0,0];
        val = {t,centroid,f_centroid,sigma_array, 9999, fcentroid_array,convergence_rate,error_array,sigma_star_array,success_rate,delta_array,FOUR_COUNT,0};
        return 
        
    end

    
    % (1+1)-ES to obtain GP traning set with 1/5-rule
    if T <= TRAINING_SIZE
        % offspring genneration
        centroid_temp = centroid + sigma*randn(n,1);
        f_centroid_temp = f(centroid_temp);
        T = T + 1;
        xTrain(:, T) = centroid_temp;                
        fTrain(T) = f_centroid_temp;
        if f_centroid_temp < f_centroid
            centroid = centroid_temp;
            f_centroid = f_centroid_temp;
            sigma = sigma*exp(0.8/D);
        else 
            sigma = sigma*exp(-0.2/D);
        end
        t = t + 1;
        % Store data
        fcentroid_array(t) = f_centroid;
        sigma_array(t) = sigma;
        sigma_star_array(t) = n*sigma/norm(centroid);
        
    % (mu/mu, lambda)-ES use GP estiate 
    else  
        % update theta  
        theta = sigma*LENGTH_SCALE*sqrt(n);                          % length scale for GP
        % offspring_generation (lambda offspring) 
        z = randn(n,lambda);
        fyep = gp(xTrain(:,T-TRAINING_SIZE+1:T), fTrain(T-TRAINING_SIZE+1:T), repmat(centroid,1,lambda)+sigma*z, f_centroid, theta);
        for k=1:1:lambda
            fy_true(k) = f(centroid+z(:,k));
        end
        error_array(iteration) = sqrt(var(fyep-fy_true)/var(fy_true-f(centroid)));
        % sort fyep (smaller first)
        [~, sorted_order] = sort(fyep);
        z = z(:,sorted_order);
        centroid_temp = centroid + sigma*mean(z(:,1:mu),2);
        fep_centroid_temp = gp(xTrain(:,T-TRAINING_SIZE+1:T), fTrain(T-TRAINING_SIZE+1:T), centroid_temp, f_centroid, theta);
            % To calculate four probs 
            fep_y_array(iteration) = fep_centroid_temp;
            f_y_array(iteration) = f(centroid_temp);
            f_x_array(iteration) = f_centroid;
        if fep_centroid_temp >= f_centroid      % 1. [estimated inferior] 
            sigma = sigma * exp(-C3/D);
        else
            f_centroid_temp = f(centroid_temp);
            T =  T + 1;
            xTrain(:, T) = centroid_temp;                
            fTrain(T) = f_centroid_temp;
            if f_centroid_temp < f_centroid     % 2.1 [estimated superior & actual superior] 
                centroid = centroid_temp;
                f_centroid = f_centroid_temp;
                sigma = sigma * exp(C1/D);
            else                                % 2.2 [estimated superior & actual inferior] 
                sigma = sigma * exp(-C2/D);
            end
            % Store data
            t = t + 1;
            fcentroid_array(t) = f_centroid;
            sigma_array(t) = sigma;
            sigma_star_array(t) = n*sigma/norm(centroid);
        end
         
    end
        
        
    iteration = iteration + 1;

    
end


    % convergence rate (overall)
    t_start = TRAINING_SIZE;
    if(fname==1)
        convergence_rate = -n*sum(log(fcentroid_array(t_start+2:t)./fcentroid_array(t_start+1:t-1)))/length(fcentroid_array(t_start+1:t-1));
    elseif(fname==2)
        convergence_rate = -n/2*sum(log(fcentroid_array(t_start+2:t)./fcentroid_array(t_start+1:t-1)))/length(fcentroid_array(t_start+1:t-1));
    elseif(fname==3)
        convergence_rate = -n/3*sum(log(fcentroid_array(t_start+2:t)./fcentroid_array(t_start+1:t-1)))/length(fcentroid_array(t_start+1:t-1));
    elseif(fname==4)
        convergence_rate = -n/2*sum(log(fcentroid_array(t_start+2:t)./fcentroid_array(t_start+1:t-1)))/length(fcentroid_array(t_start+1:t-1));
    elseif(fname==5)
    	convergence_rate = -n/2*sum(log(fcentroid_array(t_start+2:t)./fcentroid_array(t_start+1:t-1)))/length(fcentroid_array(t_start+1:t-1));
    end
    % success rate
    success_rate = sum(fcentroid_array(t_start:t-1)>fcentroid_array(t_start+1:t))/length(fcentroid_array(t_start:t-1));
    delta_array = -delta_array;
    % Get four probs 
    %%%%%%%%%%%%%%%%%%%%%%%%  
    range_parent = t_start+1:iteration-1;
    range_offspring = t_start+2:iteration;
    true_superior = f_y_array(range_offspring) < f_x_array(range_parent);
    predicted_superior = fep_y_array(range_offspring) < f_x_array(range_parent);
    tn = sum(~true_superior & ~predicted_superior);
    fp = sum(~true_superior & predicted_superior);
    fn = sum(true_superior & ~predicted_superior);
    tp = sum(true_superior & predicted_superior);
    FOUR_COUNT = [tn,fp,fn,tp];
    val = {iteration,centroid,f_centroid,sigma_array, T, fcentroid_array,convergence_rate,error_array,sigma_star_array,success_rate,delta_array,FOUR_COUNT,t/iteration};

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
    for i = 1:1:length(x)-1
        val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generalized sphere function 
function val = f6(x)
    val = (x'*x)^(para/2);
end
% quartic function with varying beta
function val = f7(x)
%     beta = 100;
    val = 0;
    for i = 1:1:length(x)-1
        val = val + para*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
    end
end

% Ellipsoids function with varying beta
function val = f8(x)
    if length(x) == 1
        val = x'.*para.*x;
    elseif length(x) >= 2
        val = x'*diag([para, ones(1,length(x)-1)])*x;
    end
end
% Schwefel's Problem 1.2
function val = f9(x)
    val = 0;
    for i = 1:1:length(x)
        val = val + sum(x(1:i))^2;
    end
end


end



function fTest = gp(xTrain, fTrain, xTest, mu, theta)
% input: 
%       xTrain(40 training pts)
%       fTrain(true objective function value)
%       xTest(multiple test pt)   
%       theta length scale of GP    
% return: the prediction of input test data

    [n, m] = size(xTrain);                                            % m:  # of training data
    ms = size(xTest, 2);

    if n==1
        delta = repmat(xTrain, 1, m)-repelem(xTrain, 1, m);
    else
        delta = vecnorm(repmat(xTrain, 1, m)-repelem(xTrain, 1, m));  %|x_ij = train_i-train_j| 
    end
    K = reshape(exp(-(delta/theta).^2/2), m , m);                     % K
    if n==1
        deltas = repmat(xTrain, 1, ms)-repelem(xTest, 1, m);
    else
        deltas = vecnorm(repmat(xTrain, 1, ms)-repelem(xTest, 1, m)); %|x_ij = train_i-test_j| euclidean distance
    end
    Ks = reshape(exp(-(deltas/theta).^2/2), m, ms);                        % K_star 
    fTest = (mu + Ks'*(K\(fTrain'-mu)))';
%     [n, m] = size(xTrain);                                       % m:  # of training data
% 
%     delta = vecnorm(repmat(xTrain, 1, m)-repelem(xTrain, 1, m)); %|x_ij = train_i-train_j| 
%     K = reshape(exp(-delta.^2/theta^2/2), m , m);                % K
% 
%     deltas = vecnorm(xTrain-repelem(xTest, 1, m));               %|x_ij = train_i-test_j| euclidean distance
%     Ks = exp(-(deltas/theta).^2/2)';                             % K_star             

%     deltass = vecnorm(repmat(xTest, 1, m)-repelem(xTest, 1, m));
%     Kss = reshape((exp(-deltass.^2/theta^2/2)), m , m);
%     Kinv = inv(K);       
%       mu = fTrain(40);
    
%     fTest = min(fTrain) + Ks'*(K\(fTrain'-min(fTrain)));
%     fTest = mu + Ks'*(K\(fTrain'-mu));


end

