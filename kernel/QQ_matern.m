function qqt = QQ_matern(u1, v1, nu)

u1 = reshape(u1, length(u1), 1);
v1 = reshape(v1, length(v1), 1);
u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);

d = abs(u - v);
C1 = 1 / (gamma(nu) * (2^(nu-1)));
C2 = sqrt(2 * nu);
qqt = C1 * (C2 * d).^nu .* besselk(nu, C2 * d);
qqt(d==0) = 1; % use Taylor series to derive d==0

end