function rst = RS(u1, v1, l)

u1 = reshape(u1, length(u1), 1);
v1 = reshape(v1, length(v1), 1);
u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);

rst = (3.*(u-v)./(l.^4)-(u-v).^3/(l.^6)).*exp(-(u-v).^2./(2*l.^2));

end