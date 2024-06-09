function Ybar=FODwassbary(coord,Y,weight,p,THD,epsilon,epsilon1)
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu
%Wasserstein Barycenter Computations
[Sorder,Q]=size(Y);
SPHARM_Order = (sqrt(8*Sorder+1)-3)/2;
Lv=size(coord,1);
Ymax=max(max(abs(Y)));
BS0 = matrixB(coord,SPHARM_Order);
p1=max(BS0*Y(:,1),0);
peakind=find(p1==max(p1));
peakcoord=coord(peakind(1),:);
Prodist=sum(repmat(peakcoord,Lv,1).*coord,2);
uppspindex=find(Prodist>=0);
upsphcoord=coord(uppspindex,:);
BS =BS0(uppspindex,:);

Delta=cell(Q,1);
scal=zeros(1,Q);
flagq=zeros(Q,1);
K=size(upsphcoord,1);
label=zeros(K,1);
for q=1:Q
px=BS*Y(:,q);
px=max(px,0);
xind=find(px>THD);
if length(xind)>0
flagq(q)=1;
pxth=px(xind);
%scal(q)=sum(abs(BS0*Y(:,q)));
scal(q)=sqrt(sum(Y(:,q).^2));
pxth=pxth/sum(pxth);
Delta{q}=double([upsphcoord(xind,:) pxth]);
x=upsphcoord(xind,:);
dist=zeros(size(x,1),K);
for m=1:size(x,1)
    for k=1:K
        inv=abs(sum(x(m,:).*upsphcoord(k,:))/(norm(x(m,:),2)*norm(upsphcoord(k,:),2)));
        dist(m,k)=180*abs(acos(inv))/pi;
    end
end
d=min(dist);
ind=find(d<20);
label(ind)=1;
end
end
indf=find(flagq>0);
Delta=Delta(indf);
weightk=weight(indf);
if length(indf)>1
index=find(label>0);
rescale=sum(scal(indf).*weightk);
BS1 = matrixB(upsphcoord(index,:),SPHARM_Order);
Ld=1:K;
Indexc=setdiff(Ld',index);
BS2 = matrixB(upsphcoord(Indexc,:),SPHARM_Order);

Ybar=SphWasscenterlinprog(Delta,double(upsphcoord(index,:)),BS1,BS2,weightk,p,Ymax,epsilon,epsilon1);
%ybar=BS0*Ybar;
%nybar=sum(abs(ybar));
nybar=sqrt(sum(Ybar.^2));
if nybar > 0
   Ybar=(Ybar*rescale)/nybar;
end

elseif length(indf)>0
Ybar=weight(indf)*Y(:,indf);
else
    Ybar=zeros(Sorder,1);
end

end