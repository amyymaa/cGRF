function [dw1, dw2] = dW(u1, l, a1, a2, coeffs)

c11 = coeffs(1); c12 = coeffs(2); c21 = coeffs(3); c22 = coeffs(4);

if c21 == 0 && c22 == 0
    C = LkL(a1, a1, l, c11, c12, c11, c12);
    dw1 = dkL(u1, a1, l, c11, c12) ./ C; 
    dw2 = 0;
elseif c11 == 0 && c12 == 0
    C = LkL(a2, a2, l, c21, c22, c21, c22);
    dw1 = 0; 
    dw2 = dkL(u1, a2, l, c21, c22) ./ C;
else
    C = LkL(a1, a1, l, c11, c12, c11, c12).*LkL(a2, a2, l, c21, c22, c21, c22) - ...
        LkL(a1, a2, l, c11, c12, c21, c22).*LkL(a2, a1, l, c21, c22, c11, c12);
    dw1 = ( dkL(u1, a1, l, c11, c12).*LkL(a2, a2, l, c21, c22, c21, c22) - ...
          dkL(u1, a2, l, c21, c22).*LkL(a2, a1, l, c21, c22, c11, c12) ) ./ C;
    dw2 = ( dkL(u1, a2, l, c21, c22).*LkL(a1, a1, l, c11, c12, c11, c12) - ...
          dkL(u1, a1, l, c11, c12).*LkL(a1, a2, l, c11, c12, c21, c22) ) ./ C;
end

end

%%

function kl = kL(x, y, l, c11, c12)
kl = c11.*(x-y)./(l.^2).*exp(-(x-y).^2./(2*l.^2)) + ...
     c12.*exp(-(x-y).^2./(2*l.^2));
end

% function dkl = dkL(u1, v1, l, c11, c12)
% dkl = c11.*RR(u1, v1, l) + c12.*QR(v1, u1, l)';
% end
function dkl = dkL(x, y, l, c11, c12)
dkl = c11.*(1./(l.^2)-(x-y).^2/(l.^4)).*exp(-(x-y).^2./(2*l.^2)) + ...
      c12.*(y-x)./(l.^2).*exp(-(y-x).^2./(2*l.^2));
end

function lkl = LkL(x, y, l, c11, c12, c21, c22)
lkl = c11.*c21.*(1./(l.^2)-(x-y).^2/(l.^4)).*exp(-(x-y).^2./(2*l.^2)) + ...
      c12.*c22.*exp(-(x-y).^2./(2*l.^2)) + ...
      c21.*c12.*(x-y)./(l.^2).*exp(-(x-y).^2./(2*l.^2)) + ...
      c11.*c22.*(y-x)./(l.^2).*exp(-(y-x).^2./(2*l.^2));
end
