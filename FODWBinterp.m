function Sphout = FODWBinterp(weight,DecomFOD,coord,Sorder,Clusters,Cweights)
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright: The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu


%'FODFAinterp' is Wasserstein Barycenter for FOD function interpolations
%Inputs:
%'weight' represents the weights for interpolation
%'DecomFOD' are the spherical harmonics coefficients of the peak lobes
%'coord' comes from the dense spherical mesh
%'Sorder' is the number of spherical harmonics coefficients 
%'Clusters' and 'Cweights' are the labels and weights for all peak-lobe fields

%Outputs:
%'Sphout' are spherical harmonics coefficients of the interpolated FOD function

NumPeaks=size(Clusters,2);
Sph=zeros(Sorder,NumPeaks);

for num = 1:NumPeaks
    Plist=Clusters(:,num);
    Pind=find(Plist>0);
    if length(Pind)>0
       weightk=weight(Pind');
       if sum(weightk)>0
       weightk=weightk/sum(weightk);
       Yint=zeros(Sorder,length(Pind));
       for k=1:length(Pind)
           nk=Plist(Pind(k));
           SpPk=DecomFOD(:,Pind(k));
           Spnumk=SpPk(((nk-1)*Sorder+1):(nk*Sorder));
           Yint(:,k)=Cweights(Pind(k),nk)*Spnumk;
       end
       Ymean=FODwassbary(coord,Yint,weightk,2,0.04,0.001,0.001);
       Sph(:,num)=Ymean; 
       end
    end
end

Sphout=sum(Sph,2);
end
