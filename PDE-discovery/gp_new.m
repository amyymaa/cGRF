% Interface for training a GP using the SE Kernel
% 
% [mu, K] = gp_grad(X, Y, DY)
% 
% Input
% X: n by d matrix representing n training points in d dimensions
% Y: training values corresponding Y = f(X)
% Output
% mu: mean function handle such that calling mu(XX) for some predictive points XX calculates the mean of the GP at XX 
% K: dense kernel matrix

function [mu,mu_1,mu_2,sigma3] = gp_new(X, Y, d1, d2, coeffs)
[ntrain, d] = size(X);

% Initial hyperparameters
ell0 = 0.5*sqrt(d);
s0 = std(Y);
sig0 = 5e-2*s0; 
beta = 1e-6;
% beta = eps;

% Train GP 
cov = @(hyp) cQQ(X, X, hyp, d1, d2, coeffs);
lmlfun = @(x) lml_exact(cov, Y, x, beta);
hyp = struct('cov', log([ell0, s0]), 'lik', log(sig0));
% hyp = struct('cov', log([ell0, s0]));
params = minimize_quiet(hyp, lmlfun, -50);
sigma = sqrt(exp(2*params.lik) + beta);
% sigma = eps;
fprintf('SE with gradients: (ell, s, sigma1) = (%.3f, %.3f, %.3f)\n', exp(params.cov), sigma)

% Calculate interpolation coefficients
sigma2 = sigma^2*ones(1, ntrain);
K = cQQ(X, X, params, d1, d2, coeffs) + diag(sigma2);
% lambda = K\Y;
lambda = K\(Y-cm(X, d1, d2, coeffs));

sigma3 = norm(lambda./diag(inv(K)))^2/ntrain;

% Function handle returning GP mean to be output
mu = @(XX) mean(XX, X, lambda, params, d1, d2, coeffs);
mu_1 = @(XX) mean_1(XX, X, lambda, params, d1, d2, coeffs);
mu_2 = @(XX) mean_2(XX, X, lambda, params, d1, d2, coeffs);

end

function ypred = mean(XX, X, lambda, params, d1, d2, coeffs)

KK = cQQ(X, XX, params, d1, d2, coeffs)';
m = cm(XX, d1, d2, coeffs); 
ypred = m + KK*lambda;

end

function ypred = mean_1(XX, X, lambda, params, d1, d2, coeffs)

KK = cQR(X, XX, params, d1, d2, coeffs)';
m = cdm(XX, d1, d2, coeffs);
ypred = m + KK*lambda;

end

function ypred = mean_2(XX, X, lambda, params, d1, d2, coeffs)

KK = cQS(X, XX, params, d1, d2, coeffs)';
m = cd2m(XX, d1, d2, coeffs);
ypred = m + KK*lambda;

end