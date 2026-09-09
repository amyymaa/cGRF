function [d2w1, d2w2] = d2W(u1, a1, a2, coeffs)

c11 = coeffs(1); c12 = coeffs(2); c21 = coeffs(3); c22 = coeffs(4);

C = LkL(a1, a1, c11, c12, c11, c12).*LkL(a2, a2, c21, c22, c21, c22) - ...
    LkL(a1, a2, c11, c12, c21, c22).*LkL(a2, a1, c21, c22, c11, c12);
d2w1 = ( d2kL(u1, a1, c11, c12).*LkL(a2, a2, c21, c22, c21, c22) - ...
     d2kL(u1, a2, c21, c22).*LkL(a2, a1, c21, c22, c11, c12) ) ./ C;
d2w2 = ( d2kL(u1, a2, c21, c22).*LkL(a1, a1, c11, c12, c11, c12) - ...
     d2kL(u1, a1, c11, c12).*LkL(a1, a2, c11, c12, c21, c22) ) ./ C;

end

%%

function kl = kL(x, y, c11, c12)
kl = c11.*(x-y).*exp(-(x-y).^2./2) + ...
     c12.*exp(-(x-y).^2./2);
end

function dkl = dkL(x, y, c11, c12)
dkl = c11.*(1-(x-y).^2).*exp(-(x-y).^2./2) + ...
      c12.*(y-x).*exp(-(y-x).^2./2);
end

function d2kl = d2kL(x, y, c11, c12)
d2kl = c11.*(3.*(y-x)-(y-x).^3).*exp(-(y-x).^2./2) + ...
       c12.*((-1)+(y-x).^2).*exp(-(y-x).^2./2);
end

function lkl = LkL(x, y, c11, c12, c21, c22)
lkl = c11.*c21.*(1-(x-y).^2).*exp(-(x-y).^2./2) + ...
      c12.*c22.*exp(-(x-y).^2./2) + ...
      c21.*c12.*(x-y).*exp(-(x-y).^2./2) + ...
      c11.*c22.*(y-x).*exp(-(y-x).^2./2);
end
