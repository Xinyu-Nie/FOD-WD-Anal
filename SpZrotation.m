function U = SpZrotation(maxOrder,alfa)
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

N=(maxOrder+1)*(maxOrder+2)/2;

U=zeros(N,N);
U(1,1)=1;
t=1;
for order = 2:2:maxOrder
    len=2*order+1;
    A=zeros(len,len);
    A(order+1,order+1)=1;
    for m=1:order
         A(order+1+m,order+1+m)=exp(i*m*alfa);
         A(order+1-m,order+1-m)=exp(-i*m*alfa);
    end
    U(t+1:t+len,t+1:t+len)=A;
    t=t+len;
end

end
