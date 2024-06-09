function dist=WassdistSph(x,P1,y,P2,p)
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu
%  Computaion the Wasserstein distance between FOD functions by solving a
%  linear programming.
M=size(x,1);
N=size(y,1);
Dist=zeros(M,N);
for m=1:M
    for k=1:N
        inv=abs(sum(x(m,:).*y(k,:))/(norm(x(m,:),2)*norm(y(k,:),2)));
        Dist(m,k)=abs(acos(inv));
    end
end
f=zeros(M*N,1);
for m=1:M
    f(((m-1)*N+1):m*N)=(Dist(m,:).^p)';
end


P=P1;
P(end+1:end+N)=P2;
A1=sparse(M,M*N);
for m=1:M
    A1(m,((m-1)*N+1):m*N)=1;
end
A2=sparse(N,M*N);
for m=1:M
    A2(:,((m-1)*N+1):m*N)=eye(N);
end
Aeq=sparse(M+N,M*N);
Aeq(1:M,:)=A1;
Aeq((M+1):end,:)=A2;
A=-speye(M*N);
gamma = linprog(f,A,zeros(M*N,1),Aeq,P);
if length(gamma)>0
   dsum=sum(gamma.*f);
else
    dsum=1;
end
dist=dsum.^(1/p);

end