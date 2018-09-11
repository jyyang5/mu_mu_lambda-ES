function val = mml(f,x0,sigma_star,sigma_ep_star,lambda,sigma0,NUM_OF_ITERATIONS)
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

centroid_array = zeros(n,100);
fcentroid_array = zeros(1,100);
sigma_array = zeros(1,1000);

y = zeros(n,lambda);                        % lambda offspring solution with dim n
fy = zeros(lambda);                         % objective function value of y                                              
t = 1;       
centroid = mean(x0, 2);                     % centroid of parent set, size = [n, 1]
f_centroid = f(centroid);                   % fx of centroid

centroid_array(:,t) = centroid;
fcentroid_array(t) = f_centroid;

while((t < NUM_OF_ITERATIONS) && f_centroid > 10^(-8))
    % take centroid of parent 
    dist = norm(centroid);                  % distance to optimal
    sigma = sigma_star/n*dist;              % mutation strength/step size(temp)  
    
    sigma_array(t) = sigma;
    
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
    centroid = mean(y(:,1:mu), 2);
    f_centroid = f(centroid);
    t = t+1;
    
    centroid_array(:,t) = centroid;
    fcentroid_array(t) = f_centroid;
    
end

val = {t,centroid_array,fcentroid_array,sigma_array}; 
end
