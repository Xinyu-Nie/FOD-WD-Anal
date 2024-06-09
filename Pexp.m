function vt=Pexp(v1,v)       

normv=sqrt(sum(v.^2));
if normv>0.0001
vt=cos(normv)*v1+sin(normv)*v/normv;
else
    vt=v1;
end

end