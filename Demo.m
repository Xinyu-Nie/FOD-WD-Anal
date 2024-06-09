%This is a demo to demonstrate the computation of FOD interpolation.
%All required data is in DATA.mat.
%The FOD functions are processed with peak detection and decomposition.
%The first figure is the interpolation by the Fast Approximation.
%The second figure is the Wasserstein space-based interpolation.
%Variables:
%Numpeaks is the maximum number of peak lobes for FOD functions.
%'Sorder' is the number of coefficients for Spherical Harmonics.
%'SPHARM_Order' is the highest order of the shperical harmonics. 
% 'R2C1','C2R1','UY90p','UY90n' are matrices used to rotate shperical harmonics.
%'coord' and 'trg' come from a dense spherical mesh, and 'coordM' and 'trgM'
% represent a sparse mesh for fast computation.
%'VertNbr' and 'VertNbronering' are neighboring structures of the dense%mesh.
%'FODrecon' is the coefficients of the FOD functions


load('DATA.mat')
THD=0.01;
Numpeaks=3;

decomfod=zeros(Sorder*Numpeaks,8);
onecube=zeros(8,3);
pks=zeros(8,Numpeaks);
pkres=zeros(8,Numpeaks);
dirvol=zeros(3,8,Numpeaks);
fracvol=zeros(8,Numpeaks);
t=0;
for a=0:1:1
        for b=0:1:1
            for c=0:1:1
                t=t+1;
                [dirvol(:,t,:),fracvol(t,:),pks(t,:),pkres(t,:)] = FODcsdpkdet(FODrecon(:,t),BS,THD,Numpeaks,coord,VertNbronering); %FOD peak dection based on Spherical Deconvolution
                decomfod(:,t) = foddecomcstr (FODrecon(:,t),VertNbr,BS,pks(t,:),pkres(t,:),SPHARM_Order,coord);   %FOD peaklobes decomposition
                onecube(t,:)=[a b c];
            end
        end
end


figure
for t=1:8
     y=BS*FODrecon(:,t);
     soln = (coord.*repmat(y,1,3))*0.2;
     color0 = abs(coord)+0.01;
     coord0 = soln+repmat([onecube(t,1) onecube(t,2) onecube(t,3)],size(coord,1),1);
     p = ViewMesh(coord0,trg);
     view(90,0)
     set(p,'FaceVertexCData',color0(:,[2 1 3]))
     set(p,'facecolor','interp') 
     alpha(0.3)
     hold on 
end

for a=0:0.5:1
    for b=0:0.5:1
        for c=0:0.5:1
points=[a b c];
P=round(points);
if ~(points(1)==P(1) && points(2)==P(2) && points(3)==P(3))
weights=LinearWeights(points);
[Clabel,Cweights] = Peaklobesregroup(decomfod,fracvol,Numpeaks,THD,Sorder,coordM); %Regroup of the FOD peak lobes
Sph = FODFAinterp(weights,decomfod,pks,coord,Sorder,SPHARM_Order,R2C1,C2R1,UY90p,UY90n,Clabel,Cweights); %Fast Approximation Interpolation
y=BS*Sph;
y = max(0,y);
soln = (coord.*repmat(y,1,3))*0.2;
color0 = abs(coord)+0.1;
coord0 = soln+repmat(points,size(coord0,1),1);
p = ViewMesh(coord0,trg);
view(90,0)
set(p,'FaceVertexCData',color0(:,[2 1 3]))
set(p,'facecolor','interp')
alpha(0.3)
hold on
end
        end
    end
end
title('Fast Approximation Interpolation')

figure
for t=1:8
     y=BS*FODrecon(:,t);
     soln = (coord.*repmat(y,1,3))*0.2;
     color0 = abs(coord)+0.01;
     coord0 = soln+repmat([onecube(t,1) onecube(t,2) onecube(t,3)],size(coord0,1),1);
     p = ViewMesh(coord0,trg);
     view(90,0)
     set(p,'FaceVertexCData',color0(:,[2 1 3]))
     set(p,'facecolor','interp') 
     alpha(0.3)
     hold on 
end

for a=0:0.5:1
    for b=0:0.5:1
        for c=0:0.5:1
points=[a b c];
P=round(points);
if ~(points(1)==P(1) && points(2)==P(2) && points(3)==P(3))
weights=LinearWeights(points);
[Clabel,Cweights] = Peaklobesregroup(decomfod,fracvol,Numpeaks,THD,Sorder,coordM);
Sph = FODWBinterp(weights,decomfod,coordM,Sorder,Clabel,Cweights);  %Wasserstein Barycenter Interpolation
y=BS*Sph;
y = max(0,y);
soln = (coord.*repmat(y,1,3))*0.2;
color0 = abs(coord)+0.1;
coord0 = soln+repmat(points,size(coord0,1),1);
p = ViewMesh(coord0,trg);
view(90,0)
set(p,'FaceVertexCData',color0(:,[2 1 3]))
set(p,'facecolor','interp')
alpha(0.3)
hold on
end
        end
    end
end
title('Wasserstein Barycenter Interpolation')
