function SPHARcom = foddecomcstr (SPHARMs,VertNbr,BS,Pks,Pkres,SPHARM_Order,coord)
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright:The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

lambda1=10;
lambda2=10;
Numpeaks=numel(find(Pks>0));

VertNbronering=VertNbr{1};
VertNbrtworing=VertNbr{2};
%VertNbrthreering=VertNbr{3};

if SPHARM_Order>10
  % VertNbrring = VertNbr{1};   %adaptive choose of suppress constraint
   VertNbrringSup = VertNbr{3};
else
   %VertNbrring = VertNbr{2};  
   VertNbrringSup = VertNbr{4};
end
peakNbr=cell(Numpeaks,1);
CM=cell(Numpeaks,1);
BM=CM;
Lu=size(BS,2);
LenC=zeros(Numpeaks,1);
LenVoid=zeros(Numpeaks,1);
bm=cell(Numpeaks,1);
for k=1:Numpeaks
    if Pks(k)>0
        onering=VertNbronering{Pks(k)};
        tworing=VertNbrtworing{Pks(k)};
       % threering=VertNbrthreering{Pks(k)};
        %threering=setdiff(threering,tworing);
        tworing=setdiff(tworing,onering);
        tworing=setdiff(tworing,Pks(k));
        vonering=coord(onering',:);
        vtworing=coord(tworing',:);
       % vthreering=coord(threering',:);
        p0v=matrixB(coord(Pks(k),:),SPHARM_Order)*SPHARMs;
        p1v=matrixB(vonering,SPHARM_Order)*SPHARMs;
        p2v=matrixB(vtworing,SPHARM_Order)*SPHARMs;
        %p3v=matrixB(vthreering,SPHARM_Order)*SPHARMs;
        if p0v>0
            p1v=p1v/p0v;
            p2v=p2v/p0v;
        end
        Ind1=onering(p1v'>0.6);
        Ind2=tworing(p2v'>0.6);
        Indep=[Pks(k) Ind1 Ind2];

        onering=VertNbronering{Pkres(k)};
        tworing=VertNbrtworing{Pkres(k)};
        %threering=VertNbrthreering{Pkres(k)};
       % threering=setdiff(threering,tworing);
        tworing=setdiff(tworing,onering);
        tworing=setdiff(tworing,Pkres(k));
        vonering=coord(onering',:);
        vtworing=coord(tworing',:);
        %vthreering=coord(threering',:);
        p0v=matrixB(coord(Pkres(k),:),SPHARM_Order)*SPHARMs;
        p1v=matrixB(vonering,SPHARM_Order)*SPHARMs;
        p2v=matrixB(vtworing,SPHARM_Order)*SPHARMs;
        %p3v=matrixB(vthreering,SPHARM_Order)*SPHARMs;
        if p0v>0
            p1v=p1v/p0v;
            p2v=p2v/p0v;
        end
        Ind1=onering(p1v'>0.6);
        Ind2=tworing(p2v'>0.6);
        Inder=[Pkres(k) Ind1 Ind2];
        %ring=VertNbrring{Pks(k)};
        %ringre=VertNbrring{Pkres(k)};
        %ring=[ring ringre];

        ring=[Indep Inder];   %Peak lobes constraint
        peakNbr{k}=ring';
        CM{k}=BS(ring',:);
        ringsup=VertNbrringSup{Pks(k)};
        ringresup=VertNbrringSup{Pkres(k)};
        ringsup=[ringsup ringresup];
        BSV=BS;
        BSV(ringsup',:)=[];
        BSV=BSV(1:5:end,:);
        BM{k}=BSV;
        LenC(k)=length(ring);
        LenVoid(k)=size(BSV,1);
        bm{k}= CM{k}*SPHARMs;
    end
end
A0=zeros(sum(LenC)+sum(LenVoid)+Lu,Lu*Numpeaks);
b0=zeros(sum(LenC)+sum(LenVoid)+Lu,1);
b0(end-Lu+1:end)=SPHARMs;

Lksum=0;
for k=1:Numpeaks
    Ak=zeros(LenC(k)+LenVoid(k),Lu);
    bk=zeros(LenC(k)+LenVoid(k),1);
    Ak(1:LenC(k),1:Lu)=lambda1*CM{k};
    bk(1:LenC(k))=lambda1*bm{k};
    Ak((LenC(k)+1):(LenC(k)+LenVoid(k)),1:Lu)=lambda2*BM{k};
    Lk=LenC(k)+LenVoid(k);
   
  A0((1+Lksum):(Lksum+Lk),(1+(k-1)*Lu):(k*Lu))=Ak;
  b0((1+Lksum):(Lksum+Lk))=bk;
  A0(end-Lu+1:end,(1+(k-1)*Lu):(k*Lu))=eye(Lu);
  Lksum=Lksum+Lk;
end

[SPHARcom,flag]=lsqr(A0,b0);
end