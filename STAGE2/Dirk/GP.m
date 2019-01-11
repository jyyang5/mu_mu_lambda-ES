function fTest = GP(xTest, xTrain, fTrain, mu, theta)
% compute the GP mean biased toward the smallest value in the training set

[n, m] = size(xTrain);
ms = size(xTest, 2);

% compute the covariance matrix
if n==1
    d = repmat(xTrain, 1, m)-repelem(xTrain, 1, m);
else
    d = vecnorm(repmat(xTrain, 1, m)-repelem(xTrain, 1, m));
end
K = reshape(exp(-(d/theta).^2/2), m, m);
%Kinv = inv(K);

if n==1
    d = repmat(xTrain, 1, ms)-repelem(xTest, 1, m);
else
    d = vecnorm(repmat(xTrain, 1, ms)-repelem(xTest, 1, m));
end
Ks = reshape(exp(-(d/theta).^2/2), m, ms);

%fTest = min(fTrain)+Ks'*Kinv*(fTrain'-min(fTrain));
fTest = (mu+Ks'*(K\(fTrain'-mu)))';

end