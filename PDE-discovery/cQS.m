function K = cQS(u1, v1, params, a1, a2, coeffs)

l = exp(params.cov(1)); 
s = exp(params.cov(2));

u1 = reshape(u1, length(u1), 1);
v1 = reshape(v1, length(v1), 1);
[u1w1, u1w2] = W(u1, a1, a2, coeffs);
[v1w1, v1w2] = W(v1, a1, a2, coeffs);
[du1w1, du1w2] = dW(u1, a1, a2, coeffs);
[dv1w1, dv1w2] = dW(v1, a1, a2, coeffs);
[d2u1w1, d2u1w2] = d2W(u1, a1, a2, coeffs);
[d2v1w1, d2v1w2] = d2W(v1, a1, a2, coeffs);

c11 = coeffs(1); c12 = coeffs(2); c21 = coeffs(3); c22 = coeffs(4);

K = s^2 * ( QS(u1, v1, l) - ...
    repmat(d2v1w1', length(u1), 1) .* kL(u1, zeros(size(v1))+a1, l, c11, c12) - ...
    repmat(d2v1w2', length(u1), 1) .* kL(u1, zeros(size(v1))+a2, l, c21, c22) - ...
    repmat(u1w1, 1, length(v1)) .* d2kL(v1, zeros(size(u1))+a1, l, c11, c12)' + ...
    repmat(u1w1, 1, length(v1)) .* repmat(d2v1w1', length(u1), 1) .* LkL(zeros(size(u1))+a1, zeros(size(v1))+a1, l, c11, c12, c11, c12) + ...
    repmat(u1w1, 1, length(v1)) .* repmat(d2v1w2', length(u1), 1) .* LkL(zeros(size(u1))+a1, zeros(size(v1))+a2, l, c11, c12, c21, c22) - ...
    repmat(u1w2, 1, length(v1)) .* d2kL(v1, zeros(size(u1))+a2, l, c21, c22)' + ...
    repmat(u1w2, 1, length(v1)) .* repmat(d2v1w1', length(u1), 1) .* LkL(zeros(size(u1))+a2, zeros(size(v1))+a1, l, c21, c22, c11, c12) + ...
    repmat(u1w2, 1, length(v1)) .* repmat(d2v1w2', length(u1), 1) .* LkL(zeros(size(u1))+a2, zeros(size(v1))+a2, l, c21, c22, c21, c22) );

end

%% functions

function kl = kL(u1, v1, l, c11, c12)
kl = c11.*QR(u1, v1, l) + c12.*QQ(u1, v1, l);
end

function dkl = dkL(u1, v1, l, c11, c12)
dkl = c11.*RR(u1, v1, l) + c12.*QR(v1, u1, l)';
end

function d2kl = d2kL(u1, v1, l, c11, c12)
d2kl = c11.*RS(v1, u1, l)' + c12.*QS(v1, u1, l)';
end

function lkl = LkL(u1, v1, l, c11, c12, c21, c22)
lkl = c11.*c21.*RR(u1, v1, l) + c12.*c22.*QQ(u1, v1, l) +...
      c21.*c12.*QR(u1, v1, l) + c11.*c22.*QR(v1, u1, l)';
end
