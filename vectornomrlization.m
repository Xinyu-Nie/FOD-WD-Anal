function vectorn = vectornomrlization(v)
vq=sqrt(sum(v.^2));
if vq>0
    vectorn = v/vq;
else
    vectorn=v;
end
end