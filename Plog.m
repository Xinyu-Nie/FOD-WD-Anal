function vt=Plog(v,v1)       

inv=sum(v.*v1);
if inv<(1-0.0001)
vt=(v1-inv*v)*abs(acos(inv))/sqrt(1-inv^2);
else
    vt=zeros(size(v));
end

end