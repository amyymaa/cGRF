
function cqq = cQQ(u1, v1, l, a1, a2, coeffs)
% u1, v1: input locations (x and x' in writing)
% l: lengthscale, prior precision
% a1, a2: boundary locations
% coeffs: [c11, c12, c21, c22], coefficients for mixed boundary constraints

u1 = reshape(u1, length(u1), 1);
v1 = reshape(v1, length(v1), 1);
[u1w1, u1w2] = W(u1, l, a1, a2, coeffs);
[v1w1, v1w2] = W(v1, l, a1, a2, coeffs);
c11 = coeffs(1); c12 = coeffs(2); c21 = coeffs(3); c22 = coeffs(4);

cqq = QQ(u1, v1, l) - ...
    repmat(v1w1', length(u1), 1) .* kL(u1, zeros(size(v1))+a1, l, c11, c12) - ...
    repmat(v1w2', length(u1), 1) .* kL(u1, zeros(size(v1))+a2, l, c21, c22) - ...
    repmat(u1w1, 1, length(v1)) .* kL(v1, zeros(size(u1))+a1, l, c11, c12)' + ...
    repmat(u1w1, 1, length(v1)) .* repmat(v1w1', length(u1), 1) .* LkL(zeros(size(u1))+a1, zeros(size(v1))+a1, l, c11, c12, c11, c12) + ...
    repmat(u1w1, 1, length(v1)) .* repmat(v1w2', length(u1), 1) .* LkL(zeros(size(u1))+a1, zeros(size(v1))+a2, l, c11, c12, c21, c22) - ...
    repmat(u1w2, 1, length(v1)) .* kL(v1, zeros(size(u1))+a2, l, c21, c22)' + ...
    repmat(u1w2, 1, length(v1)) .* repmat(v1w1', length(u1), 1) .* LkL(zeros(size(u1))+a2, zeros(size(v1))+a1, l, c21, c22, c11, c12) + ...
    repmat(u1w2, 1, length(v1)) .* repmat(v1w2', length(u1), 1) .* LkL(zeros(size(u1))+a2, zeros(size(v1))+a2, l, c21, c22, c21, c22);

end

%% functions

function kl = kL(u1, v1, l, c11, c12)
kl = c11.*QR(u1, v1, l) + c12.*QQ(u1, v1, l);
end

% function dkl = dkL(u1, v1, l, c11, c12)
% dkl = c11.*RR(u1, v1, l) + c12.*QR(v1, u1, l)';
% end

% function d2kl = d2kL(u1, v1, l, c11, c12)
% d2kl = c11.*RS(v1, u1, l)' + c12.*QS(v1, u1, l)';
% end

function lkl = LkL(u1, v1, l, c11, c12, c21, c22)
lkl = c11.*c21.*RR(u1, v1, l) + c12.*c22.*QQ(u1, v1, l) +...
      c21.*c12.*QR(u1, v1, l) + c11.*c22.*QR(v1, u1, l)';
end
