function vout=SphereKarcherMean(PX,weight,tol) 
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

[L,N]=size(PX);
vmean=mean(PX,2);
q=sqrt(sum(vmean.^2));
if q>0
   vbar=vmean/q;
else
   vbar=PX(:,1);
end
T=1000;
lambda=1;

for t=1:T
    phi=zeros(size(vbar));
    for n=1:N
        phi=phi+weight(n)*Plog(vbar,PX(:,n));
    end
    vbar=Pexp(vbar,lambda*phi);
    phinorm=sqrt(sum(phi.^2));
    if phinorm<tol
        break;
    end
end

vout=vbar;
end