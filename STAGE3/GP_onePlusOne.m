
function val = GP_onePlusOne(fname,para,x0,sigma0,NUM_OF_ITERATIONS,TRAINING_SIZE,LENGTH_SCALE)%, OPTIMAL, TARGET_DISTANCE)
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
% initialization
[n,~] = size(x0);
sigma = sigma0;                           % mutation strength(temp)
% TRAINING_SIZE = 40;

% for graphing
sigma_funEva_array = zeros(1,50000);      % store all sigma over function evaluations
x_array = zeros(n,50000);                 % store all best candiate solution
f_x = zeros(1,50000);                     % store all objective function val of candaite solution
sigma_iterate_array = zeros(1,50000);     % store sigma over iterations   
x_array(:,1) = x0;
f_x(1) = f(x0);
sigma_array(1) = sigma0;
sigma_iterate_array(1) = sigma0;
sigma_star_array = zeros(1,50000);        % store normalized step size
fep_x = zeros(1,50000);                   % store GP estimate of parent
GP_error = zeros(1,50000);                % relative error for GP
GP_error_final = zeros(1);
delta_array = zeros(1,50000);

c1 = 0.05;                                % 1/5 rule
c2 = 0.2;
c3 = 0.6;
D = sqrt(n+1);


%f = @(x) x'*x;                           % true objective function
%k = @(xy,theta) exp(-norm(xy)^2/theta/2);% square exponential (SE) with input |x-y|

x = x0;%randn(n, 1);                      % best candiate solution dim:n
fx = f(x);                                % function value
xTrain = zeros(n,50000);                  % training data for GP size 4*mu
fTrain = zeros(1,50000);
xTrain(:, 1) = x;                         % parent solutions
fTrain(1) = fx;                           % vector: value of parent solutions

iteration = 1;                                    % # of iterations
T = 1;                                    % # of distinct parent solution 


% initial sigma_satr 
dist = norm(x);
sigma_star = sigma*n/dist;
sigma_star_array(T) = sigma_star;

% Compute 4 probs
fep_y_array = zeros(1,50000);   % predicted function value of offspring 
f_y_array = zeros(1,50000);     % true function value of offspring 
f_x_array = zeros(1,50000);     % true function value of parent 


while(T < NUM_OF_ITERATIONS && f(x_array(:,T))>10^(-8))%(norm(x_array(:,T)-OPTIMAL(:,1)) > TARGET_DISTANCE))
    
    % early stopping 
    if(fname == 6)
        if(fx> 1000000 || iteration > 100000)
            success_rate = 0;
            convergence_rate = 0;
            FOUR_COUNT = [0,0,0,0];
            val = {T, x_array(:,T), fx, sigma_funEva_array, 9999, f_x, convergence_rate,GP_error,sigma_star_array,success_rate,delta_array,FOUR_COUNT,0};
            return
        end 
    elseif(fx > 5000 || iteration >= 50000)
        % if diverge -> convergence rate = 0 success rate = 0
        success_rate = 0;
        convergence_rate = 0;
        FOUR_COUNT = [0,0,0,0];
        val = {T, x_array(:,T), fx, sigma_funEva_array, 9999, f_x, convergence_rate,GP_error,sigma_star_array,success_rate,delta_array,FOUR_COUNT,0};

        return 
    end
    
    
    
    % update GP iff. there is 40 candiate solutions
    if T <= TRAINING_SIZE
        % offspring_generation
        y = x + sigma*randn(n,1);
        fy = f(y);
        T = T + 1;
        xTrain(:, T) = y;                
        fTrain(T) = fy;
        if fy < fx
            x = y;
            fx = fy;
            sigma = sigma*exp(0.8/D);
        else 
            sigma = sigma*exp(-0.2/D);
        end
        x_array(:,T) = x;
        f_x(T) = fx;
        sigma_funEva_array(T) = sigma;
        sigma_star_array(T) = sigma*n/norm(x);
    else 
        theta = sigma*LENGTH_SCALE*sqrt(n);                                 % NOtE: need to be updated every time
        fy_ep = gp(xTrain(:, T-40:T-1), fTrain(T-40:T-1), y, fx, theta);    % fitness of offspring(use GP)
        % update GP estimate for parent
        fep_x(T) = gp(xTrain(:, T-40:T-1), fTrain(T-40:T-1), x, fx, theta);
            % To calculate four probs 
            fep_y_array(iteration) = fy_ep;
            f_y_array(iteration) = f(y);
            f_x_array(iteration) = fx;
        % offspring_generation
        y = x + sigma*randn(n,1);
        % bad offspring
        if fy_ep >= fx
            sigma = sigma * exp(-c1/D);
        else
            fy = f(y);                        % fitness of offspring(use true objective fn)
            T = T + 1;
            xTrain(:, T) = y;                
            fTrain(T) = fy;
            GP_error(T) = abs(fy_ep-fy)./abs(fy-f(x));
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
            sigma_star_array(T) = sigma*n/norm(x);
        end
            
    end
    
    % new iteration     
    iteration = iteration + 1;
    
    sigma_iterate_array(T) = sigma;
    
    % store last normalized step size
    sigma_star_array(T) = sigma*n/norm(x);
        
    
    
end 
    % convergence rate (overall)
    t_start = TRAINING_SIZE;
    if(fname==1)
        convergence_rate = -n*sum(log(f_x(t_start+1:T)./f_x(t_start:T-1)))/length(f_x(t_start:T-1));
    elseif(fname==2)
        convergence_rate = -n/2*sum(log(f_x(t_start+1:T)./f_x(t_start:T-1)))/length(f_x(t_start:T-1));
    elseif(fname==3)
        convergence_rate = -n/3*sum(log(f_x(t_start+1:T)./f_x(t_start:T-1)))/length(f_x(t_start:T-1));
    elseif(fname==4)
        convergence_rate = -n/2*sum(log(f_x(t_start+1:T)./f_x(t_start:T-1)))/length(f_x(t_start:T-1));
    else
    	convergence_rate = -n/2*sum(log(f_x(t_start+1:T)./f_x(t_start:T-1)))/length(f_x(t_start:T-1));
    end
    
%     GP_error(1:length(t_start+2:T)) = abs(fep_x(t_start+2:T)-f_x(t_start+2:T))./abs(f_x(t_start+1:T-1)-f_x(t_start+2:T));
    % success rate
    success_rate = sum(f_x(TRAINING_SIZE:T-1)>f_x(TRAINING_SIZE+1:T))/length(f_x(TRAINING_SIZE:T-1));
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
    val = {T, x_array(:,T), fx, sigma_funEva_array, T, f_x, convergence_rate,GP_error,sigma_star_array,success_rate,delta_array,FOUR_COUNT,T/iteration};

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
end


