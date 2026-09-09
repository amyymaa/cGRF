function qqt = QQ_periodic(u1, v1, p, l)

u1 = reshape(u1, length(u1), 1);
v1 = reshape(v1, length(v1), 1);
u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);

qqt = exp(-2.*sin(pi.*(u-v)./p).^2./(l.^2));

end