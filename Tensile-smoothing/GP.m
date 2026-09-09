function [f2_mean, f2_sigma, evidence] = GP(y1, m1, m2, K11, K12, K22, noise)
n = size(K11, 1);
L = chol(K11 + eye(n)*noise, "lower");
alpha = L' \ (L\(y1-m1));
f2_mean = m2 + K12' * alpha;
v = L \ K12;
f2_sigma = K22 - v'*v;
evidence = - ((y1-m1)'*alpha)/2 - sum(log(diag(L))) - n*log(2*pi)/2;
end