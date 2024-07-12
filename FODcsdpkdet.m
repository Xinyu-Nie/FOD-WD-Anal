function [DirVol,FracVol,peaks,peakres] = FODcsdpkdet(Y0,BS,THD,NumPeaks,coord,VertNbronering)
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright: The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

%Inputs:
%Y0,
%BS,
%THD,
%NumPeaks,
%coord,
%VertNbronering.

DirVol = zeros(3,NumPeaks);
FracVol = zeros(1,NumPeaks);
peaks=zeros(1,NumPeaks);
peakres=peaks;

SPHARM_Order=16;
BS1=matrixB(coord,SPHARM_Order);
magnitude=500;
p=100;
Lorder=0:2:SPHARM_Order;
Lenord=2*Lorder+1;
Sorder=sum(Lenord);
stepsize=0.01;
x0=-1:stepsize:1;
x=x0(1:end-1)+stepsize/2;
Pl=zeros(length(Lorder),length(x));
for k=1:length(Lenord)
Pl(k,:)=2*pi*legendreP(Lorder(k),x);
end
y=((1+abs(x))/2).^p;
y=magnitude*y;
yy=repmat(y,length(Lorder),1);
Gl=Pl.*yy;
InGl=sum(Gl*stepsize,2);
Gy=zeros(1,Sorder);
for k=1:length(Lorder)
    if k==1
        Gy(1)=InGl(1);
    else
        Ink=sum(Lenord(1:(k-1)));
        Gy(1,Ink+1:Ink+Lenord(k))=InGl(k);
    end
end
A=zeros(size(BS1));
Lenc=size(coord,1);
for k=1:length(Lorder)
     if k==1
        A(:,1)=InGl(1)*BS1(:,1);
    else
        Ink=sum(Lenord(1:(k-1)));
        A(:,Ink+1:Ink+Lenord(k))=BS1(:,Ink+1:Ink+Lenord(k)).*repmat(Gy(1,Ink+1:Ink+Lenord(k)),Lenc,1);
    end
end
             
               y=BS*Y0;
               s=y;
               s(y<THD)=0;
               if max(s)>THD
               [Ydecon,flag1]=lsqr(A,s);
               sy=BS1*Ydecon;
               sy(y<THD)=0;
               [Dir1,Frac1,peak1,peakre1] = Peakdetect(sy,y,THD,NumPeaks,coord,VertNbronering);
               for num=1:NumPeaks
                   DirVol(:,num)=Dir1(:,num);
                   FracVol(num)=Frac1(num);
                   peaks(num)=peak1(num);
                   peakres(num)=peakre1(num);
               end
               end
            


end
