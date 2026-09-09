function [K, dKhyp] = cQQ(u1, v1, params, a1, a2, coeffs)

l = exp(params.cov(1)); 
s = exp(params.cov(2));

u1 = reshape(u1, length(u1), 1);
v1 = reshape(v1, length(v1), 1);
[u1w1, u1w2] = W(u1, a1, a2, coeffs); % remove l and alpha for simplicity
[v1w1, v1w2] = W(v1, a1, a2, coeffs);
c11 = coeffs(1); c12 = coeffs(2); c21 = coeffs(3); c22 = coeffs(4);

K = s^2 * ( QQ(u1, v1, l) - ...
    repmat(v1w1', length(u1), 1) .* kL(u1, zeros(size(v1))+a1, l, c11, c12) - ...
    repmat(v1w2', length(u1), 1) .* kL(u1, zeros(size(v1))+a2, l, c21, c22) - ...
    repmat(u1w1, 1, length(v1)) .* kL(v1, zeros(size(u1))+a1, l, c11, c12)' + ...
    repmat(u1w1, 1, length(v1)) .* repmat(v1w1', length(u1), 1) .* LkL(zeros(size(u1))+a1, zeros(size(v1))+a1, l, c11, c12, c11, c12) + ...
    repmat(u1w1, 1, length(v1)) .* repmat(v1w2', length(u1), 1) .* LkL(zeros(size(u1))+a1, zeros(size(v1))+a2, l, c11, c12, c21, c22) - ...
    repmat(u1w2, 1, length(v1)) .* kL(v1, zeros(size(u1))+a2, l, c21, c22)' + ...
    repmat(u1w2, 1, length(v1)) .* repmat(v1w1', length(u1), 1) .* LkL(zeros(size(u1))+a2, zeros(size(v1))+a1, l, c21, c22, c11, c12) + ...
    repmat(u1w2, 1, length(v1)) .* repmat(v1w2', length(u1), 1) .* LkL(zeros(size(u1))+a2, zeros(size(v1))+a2, l, c21, c22, c21, c22) );    
dKl = s^2 * ( dQQl(u1, v1, l) - ...
      repmat(v1w1', length(u1), 1) .* dkLl(u1, zeros(size(v1))+a1, l, c11, c12) - ...
      repmat(v1w2', length(u1), 1) .* dkLl(u1, zeros(size(v1))+a2, l, c21, c22) - ...
      repmat(u1w1, 1, length(v1)) .* dkLl(v1, zeros(size(u1))+a1, l, c11, c12)' + ...
      repmat(u1w1, 1, length(v1)) .* repmat(v1w1', length(u1), 1) .* dLkLl(zeros(size(u1))+a1, zeros(size(v1))+a1, l, c11, c12, c11, c12) + ...
      repmat(u1w1, 1, length(v1)) .* repmat(v1w2', length(u1), 1) .* dLkLl(zeros(size(u1))+a1, zeros(size(v1))+a2, l, c11, c12, c21, c22) - ...
      repmat(u1w2, 1, length(v1)) .* dkLl(v1, zeros(size(u1))+a2, l, c21, c22)' + ...
      repmat(u1w2, 1, length(v1)) .* repmat(v1w1', length(u1), 1) .* dLkLl(zeros(size(u1))+a2, zeros(size(v1))+a1, l, c21, c22, c11, c12) + ...
      repmat(u1w2, 1, length(v1)) .* repmat(v1w2', length(u1), 1) .* dLkLl(zeros(size(u1))+a2, zeros(size(v1))+a2, l, c21, c22, c21, c22) );    

if nargout == 2
    dKhyp = {dKl, 2*K/s};
end

end

%% functions

function kl = kL(u1, v1, l, c11, c12)
kl = c11.*QR(u1, v1, l) + c12.*QQ(u1, v1, l);
end

function dkll = dkLl(u1, v1, l, c11, c12)
dkll = c11.*dQRl(u1, v1, l) + c12.*dQQl(u1, v1, l);
end

function lkl = LkL(u1, v1, l, c11, c12, c21, c22)
lkl = c11.*c21.*RR(u1, v1, l) + c12.*c22.*QQ(u1, v1, l) + ...
      c21.*c12.*QR(u1, v1, l) + c11.*c22.*QR(v1, u1, l)';
end

function dlkll = dLkLl(u1, v1, l, c11, c12, c21, c22)
dlkll = c11.*c21.*dRRl(u1, v1, l) + c12.*c22.*dQQl(u1, v1, l) + ...
        c21.*c12.*dQRl(u1, v1, l) + c11.*c22.*dQRl(v1, u1, l)';
end

function dqql = dQQl(u1, v1, l)
u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);
dqql = exp(-(u-v).^2./(2*l.^2)) .* ((u-v).^2./(l.^3));
end

function dqrl = dQRl(u1, v1, l)
u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);
dqrl = exp(-(u - v).^2 ./ (2 * l.^2)) .* (-2 * (u - v) ./ l.^3 + (u - v).^3 ./ l.^5);
end

function drrl = dRRl(u1, v1, l)
u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);
drrl = exp(-(u - v).^2 ./ (2 * l.^2)) .* (-2 ./ l.^3 + 4 * (u - v).^2 ./ l.^5 + (1 ./ l.^2 - (u - v).^2 ./ l.^4) .* (u - v).^2 ./ l.^3);
end
