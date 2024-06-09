function Ybar=SphWasscenterlinprog(Delta,X,BS1,BS2,weight,p,Ymax,epsilon,epsilon1)
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu
%  Compution of the Wasserstein Barycenter based on linear programming
K=size(X,1);
Q=length(Delta);
Dist=cell(Q,1);
M=zeros(Q,1);
F=cell(Q,1);
Lenf=zeros(Q,1);
Gamax=zeros(Q,1);
for q=1:Q
delta=Delta{q};
x=delta(:,1:3);
Mq=size(delta,1);
dist=zeros(Mq,K);
for m=1:Mq
    for k=1:K
        inv=abs(sum(x(m,:).*X(k,:))/(norm(x(m,:),2)*norm(X(k,:),2)));
        dist(m,k)=abs(acos(inv));
    end
end
Dist{q} =dist;
MqK=K*Mq;
M(q)=Mq;
Lenf(q)=MqK;
fq=zeros(MqK,1);
for m=1:Mq
    fq(((m-1)*K+1):m*K)=repmat(weight(q),K,1).*(dist(m,:).^p)';
end
F{q}=fq;
Gamax(q)=max(delta(:,4));
if q==1
   Pconstraint=delta(:,4);
else
   Pconstraint(end+1:end+Mq,1)=delta(:,4);
end
end

LENF=sum(Lenf);
f=zeros(LENF,1);
for q=1:Q
    if q==1
        f(1:Lenf(1))=F{1};
    else
        Lenq=sum(Lenf(1:(q-1)));
        f(1+Lenq:Lenq+Lenf(q))=F{q};
    end
end
Lenp=length(Pconstraint);
Aeq1=sparse(Lenp,LENF);
for q=1:Q
    if q==1
    Aeq1(1:M(1),1:(M(1)*K)) = kron(eye(M(1)),ones(1,K));
    else
    Lenq=sum(M(1:(q-1)));
    Aeq1(Lenq+1:Lenq+M(q),(Lenq*K+1):(Lenq*K+K*M(q))) = kron(eye(M(q)),ones(1,K));
    end
end

Aeq2=sparse(K*Q,LENF);
for q=1:Q
    if q==1
        Len1=0;
    else
        Len1=sum(M(1:(q-1)));
    end
    A21 = kron(ones(1,M(q)),eye(K));
    Aeq2((q-1)*K+1:q*K,Len1*K+1:Len1*K+M(q)*K)=A21;
end

A3=[-BS1 eye(K)];
Aeq3=repmat(A3,Q,1);

Ly=size(BS1,2);
Lcontrol=size(BS2,1);
Aeq=sparse(Lenp+Q*K+Lcontrol,LENF+Ly+K+Lcontrol);
Aeq(1:Lenp,1:LENF)=Aeq1;
Aeq(Lenp+1:Lenp+Q*K,1:LENF)=Aeq2;
Aeq(Lenp+1:Lenp+Q*K,LENF+1:LENF+Ly+K)=Aeq3;

Aeq4=[-BS2 zeros(Lcontrol,K) eye(Lcontrol)];
Aeq(Lenp+Q*K+1:Lenp+Q*K+Lcontrol,LENF+1:LENF+Ly+K+Lcontrol)=Aeq4;

P=zeros(Lenp+Q*K+Lcontrol,1);
P(1:Lenp,1)=Pconstraint;
fcost=zeros(LENF+Ly+K+Lcontrol,1);
fcost(1:LENF)=f;
A=sparse(LENF+Ly+K+Lcontrol,LENF+Ly+K+Lcontrol);
A(1:LENF,1:LENF)=-speye(LENF);
lb=zeros(LENF+Ly+K+Lcontrol,1);
ub=lb;
Gmax=max(Gamax);
lb(1:LENF)=-Gmax;
ub(1:LENF)=Gmax;
lb(LENF+1:LENF+Ly)=-1.5*Ymax;
ub(LENF+1:LENF+Ly)=1.5*Ymax;
lb(LENF+Ly+1:LENF+Ly+K)=-epsilon;
ub(LENF+Ly+1:LENF+Ly+K)=epsilon;
lb(LENF+Ly+K+1:LENF+Ly+K+Lcontrol)=-epsilon1;
ub(LENF+Ly+K+1:LENF+Ly+K+Lcontrol)=epsilon1;

gamma = linprog(fcost,A,zeros(LENF+Ly+K+Lcontrol,1),Aeq,P,lb,ub);
if length(gamma)>1
Ybar=gamma(LENF+1:LENF+Ly);
else
    Ybar=zeros(Ly,1);
end
%pbar=kron(ones(1,M(1)),eye(K))*gamma(1:M(1)*K,1);


end