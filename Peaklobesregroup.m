function [Clabel,Cweights,Cmatrix] = Peaklobesregroup(DecomFOD,FracVol,Numpeak,THD,Sorder,coordM)
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright: The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

%'Peaklobesregroup' is the function that regroups the peak lobes of the FOD functions
%Inputs:
%'DecomFOD' are the spherical harmonics coefficients of the peak lobes
%'FracVol' represents the magnitudes of the peak lobes
%'Numpeak' is the maximum number of peak lobes for any FOD function
%'THD' is the cutoff threshold for FPD functions
%'Sorder' is the number of spherical harmonics coefficients 
%'coordM' comes from the coarser spherical mesh

%Outputs:
%'Clabel' and 'Cweights' are the labels and weights for all peak-lobe fields.

Ththeta=44; %Threshold for clustering
Numpoints=size(FracVol,1);
peaklobes=zeros(Sorder,Numpoints*Numpeak);
LenPk=zeros(Numpoints,1);
Lenv=0;
for k=1:Numpoints
    lent=0;
    for num=1:Numpeak
        if FracVol(k,num)>THD
           lent=lent+1;
           Lenv=Lenv+1;
           peaklobes(:,Lenv)=DecomFOD(((num-1)*Sorder+1):(num*Sorder),k);
         end
     end
     LenPk(k)=lent;
end

Cmatrix=zeros(Numpoints,Numpeak);
if Lenv>0
    peaklobes=peaklobes(:,1:Lenv);
    Simv=zeros(Lenv,Lenv);
    for k=1:Lenv
        for h=(k+1):Lenv
             [dist,scaldist]=FODdistfast(peaklobes(:,k),peaklobes(:,h),2,THD,coordM);
             Simv(k,h)=180*dist/pi;
        end
    end
    Simv=Simv+Simv';
    Y=squareform(Simv);
    Z = linkage(Y,'complete');
    %figure
    %dendrogram(Z,'ColorThreshold',25);
    idx = cluster(Z,'cutoff',Ththeta,'Criterion','distance'); 
end

if Lenv>0
    for k=1:Numpoints
        if k==1
        L=0;
        else
        L=sum(LenPk(1:(k-1)));
        end
        if LenPk(k)>0
        Cmatrix(k,1:LenPk(k))=idx((L+1):(L+LenPk(k)))';
        end
    end
end

Pklbfields=zeros(Numpoints,Lenv);
    lent=0;
    for k=1:Numpoints     
        for num=1:Numpeak
            cnum=Cmatrix(k,num);
            if cnum>0
                lent=lent+1;
                Pklbfields(k,lent)=num;
                peaklobe0=DecomFOD(((num-1)*Sorder+1):(num*Sorder),k);
                for t=1:Numpoints
                    if t~=k
                       Ind= find(Cmatrix(t,:)==cnum);
                       if length(Ind)>0
                           n0=Ind(1);
                           d0=1000;
                           for h=1:length(Ind)
                               dist=FODdistfast(peaklobe0,DecomFOD(((Ind(h)-1)*Sorder+1):(Ind(h)*Sorder),t),2,THD,coordM);
                               dist=180*dist/pi;
                               if dist<d0
                                   n0=Ind(h);
                                   d0=dist;
                               end
                           end
                           Pklbfields(t,lent)=n0;
                       end
                    end
                end
            end
        end
    end

   Clabel=unique(Pklbfields','rows','stable');
   Clabel=Clabel';
   NumIntp=size(Clabel,2);
   Cweights=zeros(size(Clabel));
   for k=1:Numpoints
       for h=1:NumIntp
           if Clabel(k,h)>0
           indh=numel(find(Clabel(k,:)==Clabel(k,h)));
           if indh>0
               Cweights(k,h)=1/indh;
           end
           end
       end
   end

end
