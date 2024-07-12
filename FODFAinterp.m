function Sphout = FODFAinterp(weight,DecomFOD,peaks,coord,Sorder,SPHARM_Order,R2C1,C2R1,UY90p,UY90n,Clusters,Cweights)
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright: The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu

%'FODFAinterp' is the Fast Approximation for FOD interpolation
%Inputs:
%'DecomFOD' are the spherical harmonics coefficients of the peak lobes
%'weight' represents the weights for interpolation
%'peaks' are the peaks indices of the peak lobes on the spherical mesh
%'coord' comes from the dense spherical mesh
%'Sorder' is the number of spherical harmonics coefficients 
%'SPHARM Order' is the highest order of spherical harmonics
%'R2C1','C2R1','UY90p','UY90n' are matrices used to rotate spherical harmonics
%'Clusters' and 'Cweights' are the labels and weights for all peak-lobe fields

%Outputs:
%'Sphout' are spherical harmonics coefficients of the interpolated FOD function


NumPeaks=size(Clusters,2);
Sph=zeros(Sorder,NumPeaks);

Lv=length(weight);
for num = 1:NumPeaks
    Plist=Clusters(:,num);
    Pind=find(Plist>0);
    Ymean=zeros(Sorder,1);
    if length(Pind)>0
       weightk=weight(Pind');
       weightk=weightk/sum(weightk);
       Vectorint=zeros(3,length(Pind));
       for m=1:length(Pind)
           Pm=peaks(Pind(m),:);
           vm=coord(Pm(Plist(Pind(m))),:)';
           if m==1
               Vectorint(:,m)=vm;
           else
               inv=sum(Vectorint(:,1).*vm);
               if inv>0
                   Vectorint(:,m)=vm;
               else
                   Vectorint(:,m)=-vm;
               end
           end
       end
       vout=SphereKarcherMean(Vectorint,weightk,0.01);
    for k=1:Lv
    if  Plist(k)>0 && weight(k)>0
        Pk=peaks(k,:);
        nk=Plist(k);
        SpPk=DecomFOD(:,k);
        Spnumk=SpPk(((nk-1)*Sorder+1):(nk*Sorder));
        vpk=coord(Pk(nk),:);
        Y0=SPHrotation(vout',vpk,Spnumk,SPHARM_Order,R2C1,C2R1,UY90p,UY90n);
        Ymean=Ymean+Y0*weight(k)*Cweights(k,num);
    end
    end
        Sph(:,num)=Ymean;
    end

end
Sphout=sum(Sph,2);
end
