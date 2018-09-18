function val = noGP(f,x0,sigma0,NUM_OF_ITERATIONS)
% initialization
% f:                  objective function value
% x0:                 initial point
% sigma0:             initial muttaion strength
% NUM_OF_ITERATIONS:  number of maximum iterations

% Return 
% 1.t:                  # of objective function calls                    
% 2.x_array(:,t):       last x
% 3.fx:                 last objective function value
% 4.sigma:              simage arrary over # of objective function calls  
% 5.x:                  parent set for x
% 6.f_x:                objective function values for parents
% 7.convergence_rate:   rate of convergence
% 8.-1:                 no GP error
% 9.sigma_star_array:   normalized step size


% OPTIMAL:            global optima
% TARGET_DISTANCE:    target distance to global optima
% example input:      fun = @(x) x' * x
%                     noGP(fun, randn(10,1),1,1000) 

[n,m] = size(x0);

x = zeros(n, 10000);                       % parent solution with dim n 
y = zeros(n,1);                            % offspring solution with dim n
sigma = zeros(1,10000);                    % mutation strength
f_x = zeros(1,10000);                      % objectie function value for selected offspring
sigma_star_array = zeros(1,10000);         % store normalized step size

D = sqrt(1+n);
c2 = -0.2;                                 % decrease step size
c3 = 0.8;                                  % increase step size

t = 1;                                     % # of iterations

sigma(1) = sigma0;                         
x(:,1) = x0;    
% initial sigma_satr 
dist = norm(x);
sigma_star = sigma(t)*n/dist;
sigma_star_array(t) = sigma_star;


while((t < NUM_OF_ITERATIONS) && f(x(:,t)) > 10^(-8))

    % offspring_generation
    y = x(:,t) + sigma(t)*randn(n,1);
    
    % objective function val 
    fy = f(y);
    fx = f(x(:,t));
    
    if(fy < fx)                             % good offspring
       x(:,t+1) = y;
       sigma(t+1) = sigma(t) * exp(c3/D);   % increase step size
       f_x(t) = fy; 
    else                                    % bad offspring
       x(:,t+1) = x(:,t);  
       sigma(t+1) = sigma(t) * exp(c2/D);   % reduce step size
       f_x(t) = fx; 
    end
    t = t + 1;
    
    % store last normalized step size
    dist = norm(x);
    sigma_star = sigma(t)*n/dist;
    sigma_star_array(t) = sigma_star;
    
    
end 

convergence_rate = -n/2*sum(log(f_x(2:t)./f_x(1:t-1)))/(t-1);

val = {t, x(:,t), fx, sigma, x, f_x, convergence_rate, -1, sigma_star_array};

end
