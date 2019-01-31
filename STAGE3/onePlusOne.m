function val = onePlusOne(fname,para,x0,sigma0,NUM_OF_ITERATIONS)
% initialization
% fname:                objective function value
% x0:                 initial point
% sigma0:             initial muttaion strength
% NUM_OF_ITERATIONS:  number of maximum iterations

% Return 
% 1.t:                  # of objective function calls                    
% 2.x_array(:,t):       last x
% 3.fx:                 last objective function value
% 4.sigma:              simage arrary over # of objective function calls  
% 5.T:                  # of objective function calls
% 6.f_x:                objective function values for parents
% 7.convergence_rate:   rate of convergence
% 8.-1:                 no GP error
% 9.sigma_star_array:   normalized step size
% 10. success_rate      success_rate


% OPTIMAL:            global optima
% TARGET_DISTANCE:    target distance to global optima
% example input:      fun = @(x) x' * x
%                     noGP(fun, randn(10,1),1,1000) 
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

x_array = zeros(n, 50000);                       % parent solution with dim n 
y = zeros(n,1);                            % offspring solution with dim n
sigma = zeros(1,50000);                    % mutation strength
f_x = zeros(1,50000);                      % objectie function value for selected offspring
sigma_star_array = zeros(1,50000);         % store normalized step size

D = sqrt(1+n);
c2 = -0.2;                                 % decrease step size
c3 = 0.8;                                  % increase step size

t = 1;                                     % # of iterations

sigma(1) = sigma0;   
x = x0;
x_array(:,1) = x0;    
% initial sigma_satr 
sigma_star_array(t) = sigma(t)*n/norm(x0);
fx = f(x0);
f_x(t) = fx;

while((t < NUM_OF_ITERATIONS) && fx > 10^(-8))
    if(fname==6)
        if fx > 10^30
            % if diverge -> convergence rate = 0 success rate = 0
            success_rate = 0;
            convergence_rate = 0;
            val = {t, x_array(:,t), fx, sigma, 99999, f_x, convergence_rate, -1, sigma_star_array,success_rate};
            return
        end
        
    elseif(fx > 5000)
        % if diverge -> convergence rate = 0 success rate = 0
        success_rate = 0;
        convergence_rate = 0;
        val = {t, x_array(:,t), fx, sigma, 99999, f_x, convergence_rate, -1, sigma_star_array,success_rate};
        return 
    end
    % offspring_generation
    y = x + sigma(t)*randn(n,1);
    
    % objective function val 
    fy = f(y);
    

    if(fy < fx)                             % good offspring
       x_array(:,t+1) = y;
       fx = fy;
       x = y;
       sigma(t+1) = sigma(t) * exp(0.8/D);   % increase step size
       f_x(t+1) = fy; 
    else                                    % bad offspring
       x_array(:,t+1) = x;  
       sigma(t+1) = sigma(t) * exp(-0.2/D);   % reduce step size
       f_x(t+1) = fx; 
    end
    t = t + 1;
    sigma_star_array(t) = sigma(t)*n/norm(x);
    
    
    % store last normalized step size
    
    
    
end 
convergence_rate = -n/2*sum(log(f_x(2:t)./f_x(1:t-1)))/(t-1);
success_rate = sum(f_x(2:t)<f_x(1:t-1))/(t-1);

val = {t, x_array(:,t), fx, sigma, t, f_x, convergence_rate, -1, sigma_star_array,success_rate};

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
