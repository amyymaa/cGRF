function cm = cm(u1, a1, a2, coeffs)

u1 = reshape(u1, length(u1), 1);
[u1w1, u1w2] = W(u1, a1, a2, coeffs);
g1 = coeffs(5); g2 = coeffs(6);

cm = u1w1 .* g1 + u1w2 .* g2;

end

% %% functions
% 
% function kl = kL(u1, v1, l, alpha, c11, c12)
% kl = c11.*QR(u1, v1, l, alpha) + c12.*QQ(u1, v1, l, alpha);
% end
% 
% % function dkl = dkL(u1, v1, l, alpha, c11, c12)
% % dkl = c11.*RR(u1, v1, l, alpha) + c12.*QR(v1, u1, l, alpha)';
% % end
% 
% % function d2kl = d2kL(u1, v1, l, alpha, c11, c12)
% % d2kl = c11.*RS(v1, u1, l, alpha)' + c12.*QS(v1, u1, l, alpha)';
% % end
% 
% function lkl = LkL(u1, v1, l, alpha, c11, c12, c21, c22)
% lkl = c11.*c21.*RR(u1, v1, l, alpha) + c12.*c22.*QQ(u1, v1, l, alpha) +...
%       c21.*c12.*QR(u1, v1, l, alpha) + c11.*c22.*QR(v1, u1, l, alpha)';
% end
