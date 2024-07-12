function [dist,scaldist]=FODdistfast(Y1,Y2,p,THD,coord)
% Contribution
%  Author: Xinyu Nie
%  Created: 2024/6/1
%  Copyright: The Neuro Image Computing Research (NICR) group at the Mark and Mary Stevens Neuroimaging 
%  and Informatics Institute of USC Laboratory of NeuroImaging 
%  USC Stevens Neuroimaging and Informatics Institute
%  email: xnie@usc.edu
% Wasserstein Distance bewteen FOD functions

Sorder=length(Y1);
SPHARM_Order = (sqrt(8*Sorder+1)-3)/2;
Lv=size(coord,1);
Meantriarea=4*pi/Lv;

uppspindex=find(coord(:,3)>=0);
upsphcoord=coord(uppspindex,:);
BS = matrixB(upsphcoord,SPHARM_Order);
px=BS*Y1;
py=BS*Y2;
px=max(px,0);
py=max(py,0);
xind=find(px>THD);
yind=find(py>THD);
pxth=px(xind);
pyth=py(yind);
scaldist=abs(2*sum(px*Meantriarea)-2*sum(py*Meantriarea));
pxth=pxth/sum(pxth);
pyth=pyth/sum(pyth);

%figure
%p = patch('Vertices',coord, 'Faces',trg,'FaceVertexCData',px+py,'facecolor','flat');
%set(p,'edgecolor','none');

if length(xind)>0 && length(yind)>0
dist=WassdistSph(upsphcoord(xind,:),double(pxth),upsphcoord(yind,:),double(pyth),p);
elseif length(xind)>0 || length(yind)>0
    dist=1;
else 
    dist=0;  
end
end
