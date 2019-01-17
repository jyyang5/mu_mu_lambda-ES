function [T, rec] = ES_weighted(f, mu, lambda, x, sigma, k, M, c1, c2, c3)

n = size(x, 1);
fx = f(x);
rec = zeros(20000, 5);          % record to store a trace of the run
arx = zeros(n+1, 4000);         % archive of evaluated points
arx(:, 1) = [fx; x];

t = 1;
T = 1;
w = log((lambda+1)/2)-log(1:mu);% weights for recombination
w = w/sum(w);

while t<2000 && fx>1.0e-08 && sigma>1.0e-08

    rec(t, 1:3) = [T, sigma, fx];

    if T<M                      % insufficient training data
        y = x+sigma*randn(n, 1);
        fy = f(y);
        if fy<fx
            x = y;
            fx = fy;
            sigma = sigma*exp(0.8/sqrt(1+n));
        else
            sigma = sigma*exp(-0.2/sqrt(1+n));
        end
        T = T+1;
        arx(:, T) = [fy; y];
        t = t+1;
        continue
    end

    xTrain = arx(2:end, T-M+1:T);
    fTrain = arx(1, T-M+1:T);
    
    z = randn(n, lambda);
    theta = k*sigma*sqrt(n);
    fTest = GP(repmat(x, 1, lambda)+sigma*z, xTrain, fTrain, fx, theta);
    
    [~, idx] = sort(fTest);
    y = x+sigma*mean(repmat(w,n,1).*z(:, idx(1:mu)), 2);

    ey = GP(y, xTrain, fTrain, fx, theta);
    fy = f(y);
    rec(t, 4) = ey;
    rec(t, 5) = fy;

    if ey<fx       
        if fy<fx
            x = y;
            fx = fy;
            sigma = sigma*exp(c1/sqrt(1+n));
        else
            sigma = sigma*exp(-c2/sqrt(1+n));
        end
        T = T+1;
        arx(:, T) = [fy; y];
    else
        sigma = sigma*exp(-c3/sqrt(1+n));
    end

    t = t+1;
end

rec = rec(1:t-1, :);

end